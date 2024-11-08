//
// Created by egor on 2/2/24.
//

#ifndef LIBGRAN_ALT_SINTER_BRIDGE_H
#define LIBGRAN_ALT_SINTER_BRIDGE_H

#include <tuple>
#include <vector>
#include <cstddef>
#include <forward_list>
#include <algorithm>
#include <numeric>

#include "../contact_force/contact_force.h"

template <typename field_value_t, typename real_t>
struct alt_sinter_functor {
    alt_sinter_functor(size_t n_part,            // Number of particles in the system
                          std::vector<field_value_t> x0,         // Initial positions
                          real_t k,                 // Normal stiffness coefficient
                          real_t gamma_n,           // Normal damping coefficient
                          real_t k_t,               // Stiffness coefficient for sticking/sliding
                          real_t gamma_t,           // Damping coefficient for sticking/sliding
                          real_t k_r,               // Stiffness coefficient for rolling
                          real_t gamma_r,           // Damping coefficient for rolling
                          real_t k_o,               // Stiffness for torsion
                          real_t gamma_o,           // Damping coefficient for torsion
                          real_t r_part,            // Radius of a particle
                          real_t mass,              // Mass
                          real_t inertia,           // Moment of inertia
                          real_t dt,                // Time step for spring update (same as integration time step for 1st order schemes)
                          field_value_t field_zero, // Zero-valued field_value_t
                          real_t real_zero,         // Zero-valued real_t
                          real_t critical_separation, // Critical separation between particles to make them necked
                          contact_force_functor<field_value_t, real_t> contact_force) : // Instance of contact force functor that handles non-bonded contacts
        n_part(n_part),
        k(k),
        gamma_n(gamma_n),
        k_t(k_t),
        gamma_t(gamma_t),
        k_r(k_r),
        gamma_r(gamma_r),
        k_o(k_o),
        gamma_o(gamma_o),
        r_part(r_part),
        mass(mass),
        inertia(inertia),
        dt(dt),
        real_zero(real_zero),
        field_zero(field_zero),
        contact_force(std::move(contact_force)),
        vertex_subsets(x0.size())
    {

        contact_springs.resize(n_part * n_part);
        std::fill(contact_springs.begin(), contact_springs.end(), std::make_tuple(field_zero, field_zero, field_zero));

        bonded_contacts.resize(n_part * n_part);
        std::fill(bonded_contacts.begin(), bonded_contacts.end(), false);

        particle_to_bond_map.resize(n_part);

        // Initialize vertex_subsets for use with the cycle prevention algorithm
        std::iota(vertex_subsets.begin(), vertex_subsets.end(), 0);

        for (size_t i = 0; i < n_part - 1; i ++) {
            for (size_t j = i+1; j < n_part; j ++) {
                if (abs((x0[i] - x0[j]).norm() - 2.0*r_part) < critical_separation) {
                    if (vertex_subsets[i] == vertex_subsets[j]) {
                        std::cout << "Warning: preventing neck insertion to avoid a cycle" << std::endl;
                        continue;
                    }

                    bonded_contacts[i*n_part + j] = true;
                    bonded_contacts[j*n_part + i] = true;
                    particle_to_bond_map[i].emplace_front(j);
                    particle_to_bond_map[j].emplace_front(i);

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
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i,
                                                         size_t j,
                                                         std::vector<field_value_t> const & x,
                                                         std::vector<field_value_t> const & v,
                                                         std::vector<field_value_t> const & theta [[maybe_unused]],
                                                         std::vector<field_value_t> const & omega,
                                                         real_t t [[maybe_unused]]) {

        if (!bonded_contacts[i*n_part + j]) [[likely]]
            return contact_force(i, j, x, v, theta, omega, t);

        field_value_t n = (x[i] - x[j]).normalized();
        real_t overlap = 2.0 * r_part - (x[i] - x[j]).dot(n);

        real_t r_part_prime = r_part - overlap / 2.0;

        real_t v_n = -(v[i] - v[j]).dot(n); // Normal relative velocity

        real_t f_n = k * overlap // Elastic contribution
                + gamma_n * v_n; // Viscous contribution

        field_value_t v_ij = v[i] - v[j] + r_part_prime * n.cross(omega[i]) + r_part_prime * n.cross(omega[j]);

        field_value_t v_t = v_ij - v_ij.dot(n) * n; // Tangential relative velocity
        field_value_t v_r = -r_part_prime * (n.cross(omega[i]) - n.cross(omega[j])); // Rolling velocity
        field_value_t v_o = r_part * (n.dot(omega[i]) - n.dot(omega[j])) * n; // Spin velocity

        field_value_t f_t = compute_shear_contribution<0>(i, j, n, k_t, gamma_t, v_t); // Sliding/sticking
        field_value_t f_r = compute_shear_contribution<1>(i, j, n, k_r, gamma_r, v_r); // Rolling
        field_value_t f_o = compute_shear_contribution<2>(i, j, n, k_o, gamma_o, v_o); // Torsion

        // Compute the torques associated with all the shear contributions
        field_value_t tau_t = r_part_prime * n.cross(f_t);
        field_value_t tau_r = r_part * n.cross(f_r);
        field_value_t tau_o = r_part * f_o;

        return std::make_pair((f_n * n + f_t) / mass, (-tau_t + tau_r + tau_o) / inertia);
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

    std::vector<bool> bonded_contacts;

private:
    const size_t n_part;
    const real_t k, gamma_n,
        k_t, gamma_t,
        k_r, gamma_r,
        k_o, gamma_o,
        r_part, mass, inertia, dt, real_zero;
    const field_value_t field_zero;
    std::vector<std::forward_list<size_t>> particle_to_bond_map;
    std::vector<std::tuple<field_value_t, field_value_t, field_value_t>> contact_springs;
    contact_force_functor<field_value_t, real_t> contact_force;

    // Data structures for the cycle prevention algorithm
    std::vector<std::pair<size_t, size_t>> undirected_graph_edges;
    std::vector<size_t> vertex_subsets;
};

#endif //LIBGRAN_ALT_SINTER_BRIDGE_H
