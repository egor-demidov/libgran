//
// Created by egor on 2/2/24.
// Edited by gurdeep on 12/30/25
//

#ifndef LIBGRAN_ALT_SINTER_BRIDGE_H
#define LIBGRAN_ALT_SINTER_BRIDGE_H

#include <tuple>
#include <vector>
#include <cstddef>
#include <forward_list>
#include <algorithm>
#include <numeric>
#include <random>

#include "../contact_force/contact_force_var_size.h"

template <typename field_value_t, typename real_t, typename matrix_t, typename BoxType, typename bond_t>
struct alt_sinter_functor {
    alt_sinter_functor(size_t n_part,            // Number of particles in the system
                          std::vector<field_value_t> x0,         // Initial positions
                          std::vector<bond_t> bond_types,// possible types of bonds
                          real_t dt,                // Time step for spring update (same as integration time step for 1st order schemes)
                          field_value_t field_zero, // Zero-valued field_value_t
                          real_t real_zero,         // Zero-valued real_t
                          std::vector<real_t> & r,  // vector of particle radii
                          std::vector<real_t> & m,  // vector of particle masses
                          BoxType & box,                // Simulation Periodic Box
                          real_t critical_separation, // Critical separation between particles to make them necked
                          real_t clay_ratio,           // ratio of clay bonds to cemented bonds
                          real_t bond_percentage,      // percent of bonds between grains
                          real_t rand_seed,             // randoms seed for random generator
                          contact_force_functor_var_size<field_value_t, real_t, matrix_t> contact_force) : // Instance of contact force functor that handles non-bonded contacts
        n_part(n_part),
        dt(dt),
        real_zero(real_zero),
        field_zero(field_zero),
        contact_force(std::move(contact_force)),
        vertex_subsets(x0.size())
    {
        // generating bonds
        contact_springs.resize(n_part * n_part);
        std::fill(contact_springs.begin(), contact_springs.end(), std::make_tuple(field_zero, field_zero, field_zero));

        bonded_contacts.resize(n_part * n_part);
        std::fill(bonded_contacts.begin(), bonded_contacts.end(), false);

        particle_to_bond_map.resize(n_part);

        // Initialize vertex_subsets for use with the cycle prevention algorithm
        std::iota(vertex_subsets.begin(), vertex_subsets.end(), 0);

        initial_normal_dist.resize(n_part * n_part);
        bond_params.resize(n_part * n_part);

        std::mt19937 gen(rand_seed);
        std::uniform_real_distribution<> dist(0.0, 1.0);

        for (size_t i = 0; i < n_part - 1; i ++) {
            for (size_t j = i+1; j < n_part; j ++) {
                // minimum image convention
                field_value_t d = x0[i] - x0[j];
                d = box.minimumImage(d);

                if (abs((d).norm() - (r[i] + r[j])) < critical_separation) {
                    // cycle pervention
                    // if (vertex_subsets[i] == vertex_subsets[j]) {
                    //     std::cout << "Warning: preventing neck insertion to avoid a cycle" << std::endl;
                    //     continue;
                    // }

                    bond_params[i*n_part+j] = bond_types[0];
                    bond_params[j*n_part+i] = bond_types[0];

                    bonded_contacts[i*n_part + j] = true;
                    bonded_contacts[j*n_part + i] = true;
                    particle_to_bond_map[i].emplace_front(j);
                    particle_to_bond_map[j].emplace_front(i);
                    
                    initial_normal_dist[i*n_part + j] = d;
                    initial_normal_dist[j*n_part + i] = -d;

                    // Update the cycle detection data structures
                    undirected_graph_edges.emplace_back(i, j);
                    // Merge the subsets
                    size_t subset_j = vertex_subsets[j];
                    for (auto & subset : vertex_subsets) {
                        if (subset == subset_j) subset = vertex_subsets[i];
                    }
                }
            }
        }

        // removing bonds so only bond percentage of bonds remain
        std::vector<std::pair<size_t, size_t>> bond_list;
        for (size_t i = 0; i < n_part - 1; i++) {
            for (size_t j = i + 1; j < n_part; j++) {
                if (bonded_contacts[i*n_part + j]) {
                    bond_list.emplace_back(i, j);
                }
            }
        }
        std::shuffle(bond_list.begin(), bond_list.end(), gen);

        real_t bond_removal_percentage = 1-bond_percentage;
        int num_bonds_remove = bond_removal_percentage * bond_list.size();

        for (int k = 0; k < num_bonds_remove; k++) {
            size_t i = bond_list[k].first;
            size_t j = bond_list[k].second;

            // Remove bond (both directions)
            bonded_contacts[i*n_part + j] = false;
            bonded_contacts[j*n_part + i] = false;

            // clear bond params
            bond_params[i*n_part + j] = bond_t();
            bond_params[j*n_part + i] = bond_t();

            // remove from adjacency list
            particle_to_bond_map[i].remove(j);
            particle_to_bond_map[j].remove(i);
        }

        // converting bonds to clay with clay ratio amount
        std::vector<std::pair<size_t, size_t>> remaining_bonds;

        for (size_t i = 0; i < n_part - 1; i++) {
            for (size_t j = i + 1; j < n_part; j++) {
                if (bonded_contacts[i*n_part + j]) {
                    remaining_bonds.emplace_back(i, j);
                }
            }
        }
        std::shuffle(remaining_bonds.begin(), remaining_bonds.end(), gen);

        int num_clay_bonds = static_cast<int>(
            std::round(clay_ratio * remaining_bonds.size())
        );

        for (int k = 0; k < num_clay_bonds; k++) {
            size_t i = remaining_bonds[k].first;
            size_t j = remaining_bonds[k].second;

            bond_params[i*n_part + j] = bond_types[1]; // clay
            bond_params[j*n_part + i] = bond_types[1];
        }
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i,
                                                         size_t j,
                                                         std::vector<field_value_t> const & x,
                                                         std::vector<field_value_t> const & v,
                                                         std::vector<field_value_t> const & theta [[maybe_unused]],
                                                         std::vector<field_value_t> const & omega,
                                                         std::vector<real_t> const & r,
                                                         std::vector<real_t> const & m,
                                                         std::vector<std::array<double, 9>> & p,
                                                         BoxType const & box,
                                                         real_t t [[maybe_unused]]) {

        if (!bonded_contacts[i*n_part + j]) [[likely]]
            return contact_force(i, j, x, v, theta, omega, r, m, p, box, t);

        bond_t current_bond = bond_params[i*n_part+j];

        // minimum image convention
        field_value_t d_raw = x[i] - x[j];
        field_value_t d = box.minimumImage(d_raw);

        field_value_t d0 = initial_normal_dist[i*n_part + j];

        field_value_t n = (d).normalized();
        //real_t overlap = (r[i] + r[j]) - (d).dot(n);

        // only add equilbirium dist for clay bonds
        real_t overlap;
        if(current_bond.type == "clay"){
            overlap = d0.norm() - d.norm() + current_bond.equilibrium_dist;
        }else{
            overlap = d0.norm() - d.norm();
        }

        // real_t r_part_prime = r[i] - 1/2 * overlap;
        real_t r_i_prime = r[i] - 1/2 * overlap;
        real_t r_j_prime = r[j] - 1/2 * overlap;
        real_t r_ij_prime = r_i_prime * r_j_prime / (r_i_prime + r_j_prime);

        const matrix_t& D = box.get_deformation_rate();
        // Subtract affine streaming velocity
        field_value_t uij = (v[i] - v[j]) - D * d_raw;
        field_value_t velocity_jump = D * (d_raw - d);

        real_t v_n = -(v[i] - v[j]).dot(n); // Normal relative velocity
        field_value_t uij_tangential = (v[i] - v[j]) - velocity_jump;

        real_t f_n = current_bond.k_n_bond * overlap // Elastic contribution
                + current_bond.gamma_n_bond * v_n; // Viscous contribution

        // Add rotational contributions
        field_value_t v_ij = uij_tangential + r_i_prime * n.cross(omega[i]) + r_j_prime * n.cross(omega[j]);

        field_value_t v_t = v_ij - v_ij.dot(n) * n; // Tangential relative velocity
        field_value_t v_r = r_ij_prime * (-n.cross(omega[i]) + n.cross(omega[j])); // Rolling velocity
        field_value_t v_o = r_ij_prime * (n.dot(omega[i]) - n.dot(omega[j])) * n; // Spin velocity

        field_value_t f_t = compute_shear_contribution<0>(i, j, n, current_bond.k_t_bond, current_bond.gamma_t_bond, v_t); // Sliding/sticking
        field_value_t f_r = compute_shear_contribution<1>(i, j, n, current_bond.k_r_bond, current_bond.gamma_r_bond, v_r); // Rolling
        field_value_t f_o = compute_shear_contribution<2>(i, j, n, current_bond.k_o_bond, current_bond.gamma_o_bond, v_o); // Torsion

        // Compute the torques associated with all the shear contributions
        field_value_t tau_t = r_i_prime * n.cross(f_t);
        field_value_t tau_r = r[i] * n.cross(f_r);
        field_value_t tau_o = r[i] * f_o;

        real_t inertia = 2.0 / 5.0 * m[i] * pow(r[i], 2.0);
        field_value_t F = f_n * n + f_t;

        // updating particle pressures for barostat
        update_particle_pressures(p, F, d, i);

        return std::make_pair((F) / m[i], (-tau_t + tau_r + tau_o) / inertia);
    }

    void reset_springs(size_t i, size_t j) {
        contact_springs[i * n_part + j] = std::make_tuple(field_zero, field_zero, field_zero);
    }

    // Computes either sliding/sticking, rolling, or torsion contribution
    // Use model_num 0 for sliding/sticking, 1 for rolling, 2 for torsion
    template<size_t model_num>
    field_value_t compute_shear_contribution(size_t i, size_t j, field_value_t const & n,
                                                   real_t stiffness, real_t damping,
                                                   field_value_t const & relative_velocity) {


        field_value_t & xi = std::get<model_num>(contact_springs[i * n_part + j]); // Access the respective spring
        field_value_t xi_new = xi - xi.dot(n) * n; // rotate the tangential spring
        // Rescale the spring to preserve its magnitude after rotation
        if (xi_new.norm() > 0.0) [[likely]] {
            xi_new.normalize();
            xi_new *= xi.norm();
        }
        // Update the spring in the spring buffer
        xi = xi_new;
        field_value_t f_0 = -stiffness * xi - damping * relative_velocity; // Compute the test force
        xi += relative_velocity * dt;

        return f_0;
    }

    [[nodiscard]]
    std::vector<std::tuple<field_value_t, field_value_t, field_value_t>> const & get_contact_springs() const {
        return contact_springs;
    }

    std::vector<bond_t> const & get_bond_params(){
        return bond_params;
    }

    std::vector<field_value_t> const & get_bond_initial_dist(){
        return initial_normal_dist;
    }

    std::array<real_t, 3> get_normal_force(){
        return contact_force.get_normal_force();
    }

    std::vector<bool> bonded_contacts;

private:
    const size_t n_part;
    const real_t dt, real_zero;
    const field_value_t field_zero;
    std::vector<std::forward_list<size_t>> particle_to_bond_map;
    std::vector<std::tuple<field_value_t, field_value_t, field_value_t>> contact_springs;
    contact_force_functor_var_size<field_value_t, real_t, matrix_t> contact_force;

    std::vector<field_value_t> initial_normal_dist;
    std::vector<bond_t> bond_params;

    // Data structures for the cycle prevention algorithm
    std::vector<std::pair<size_t, size_t>> undirected_graph_edges;
    std::vector<size_t> vertex_subsets;
};

#endif //LIBGRAN_ALT_SINTER_BRIDGE_H
