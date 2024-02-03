//
// Created by egor on 2/2/24.
//

#ifndef LIBGRAN_ALT_SINTER_BRIDGE_H
#define LIBGRAN_ALT_SINTER_BRIDGE_H

#include <tuple>
#include <vector>
#include <cstddef>
#include <forward_list>

//#include <libtimestep/rotational_step_handler/rotational_step_handler.h>
//
//template <typename field_container_t, typename field_value_t, typename real_t>
//struct sinter_step_handler {
//    sinter_step_handler(std::vector<std::tuple<field_container_t, field_container_t, field_container_t>> & contact_springs,
//                        std::vector<bool> const & bonded_contacts,
//                        std::vector<std::forward_list<size_t>> const & particle_to_bond_map,
//                        real_t r_part, size_t n_part) :
//            contact_springs(contact_springs),
//            bonded_contacts(bonded_contacts),
//            particle_to_bond_map(particle_to_bond_map),
//            r_part(r_part), n_part(n_part) {}
//
//    void increment_x(size_t n,
//                     field_value_t const & dx,
//                     typename field_container_t::iterator x_begin_itr,
//                     typename field_container_t::iterator v_begin_itr,
//                     typename field_container_t::iterator a_begin_itr,
//                     typename field_container_t::iterator theta_begin_itr,
//                     typename field_container_t::iterator omega_begin_itr,
//                     typename field_container_t::iterator alpha_begin_itr) {
//
//        for (auto const & bond : particle_to_bond_map[n]) {
//            size_t i = n;
//            size_t j = bond;
//
//            // Relative tangential displacement
//            field_value_t d_ij = *(v_begin_itr + i) - *(v_begin_itr + j);
//            field_value_t d_t = d_ij - d_ij.dot(*())
//        }
//        base_step_handler.increment_x(n, dx, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
//    }
//
//    void increment_v(size_t n,
//                     field_value_t const & dv,
//                     typename field_container_t::iterator x_begin_itr,
//                     typename field_container_t::iterator v_begin_itr,
//                     typename field_container_t::iterator a_begin_itr,
//                     typename field_container_t::iterator theta_begin_itr,
//                     typename field_container_t::iterator omega_begin_itr,
//                     typename field_container_t::iterator alpha_begin_itr) {
//
//        base_step_handler.increment_v(n, dv, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
//    }
//
//    void increment_theta(size_t n,
//                     field_value_t const & dtheta,
//                     typename field_container_t::iterator x_begin_itr,
//                     typename field_container_t::iterator v_begin_itr,
//                     typename field_container_t::iterator a_begin_itr,
//                     typename field_container_t::iterator theta_begin_itr,
//                     typename field_container_t::iterator omega_begin_itr,
//                     typename field_container_t::iterator alpha_begin_itr) {
//
//        base_step_handler.increment_theta(n, dtheta, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
//    }
//
//    void increment_omega(size_t n,
//                         field_value_t const & domega,
//                         typename field_container_t::iterator x_begin_itr,
//                         typename field_container_t::iterator v_begin_itr,
//                         typename field_container_t::iterator a_begin_itr,
//                         typename field_container_t::iterator theta_begin_itr,
//                         typename field_container_t::iterator omega_begin_itr,
//                         typename field_container_t::iterator alpha_begin_itr) {
//
//        base_step_handler.increment_omega(n, domega, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
//    }
//
//private:
//    const rotational_step_handler<field_container_t, field_value_t> base_step_handler;
//    std::vector<std::forward_list<size_t>> const & particle_to_bond_map;
//    std::vector<bool> const & bonded_contacts;
//    std::vector<std::tuple<field_container_t, field_container_t, field_container_t>> & contact_springs;
//    const size_t n_part;
//    const real_t r_part;
//};

template <typename field_value_t, typename real_t>
struct alt_sinter_functor {
    alt_sinter_functor(size_t n_part,            // Number of particles in the system
                          std::vector<field_value_t> x0,         // Initial positions
                          real_t k,                 // Normal stiffness coefficient
                          real_t gamma_n,           // Normal damping coefficient
                          real_t k_t,               // Stiffness coefficient for sticking/sliding
                          real_t gamma_t,           // Damping coefficient for sticking/sliding
                          real_t mu_s,              // Static friction coefficient for sticking/sliding
                          real_t phi_d,             // Coulomb coefficient for sticking/sliding
                          real_t k_r,               // Stiffness coefficient for rolling
                          real_t gamma_r,           // Damping coefficient for rolling
                          real_t mu_r,              // Static friction coefficient for rolling
                          real_t phi_r,             // Coulomb coefficient for rolling
                          real_t k_o,               // Stiffness for torsion
                          real_t gamma_o,           // Damping coefficient for torsion
                          real_t mu_o,              // Static friction coefficient for torsion
                          real_t phi_o,             // Coulomb coefficient for torsion
                          real_t r_part,            // Radius of a particle
                          real_t mass,              // Mass
                          real_t inertia,           // Moment of inertia
                          real_t dt,                // Time step for spring update (same as integration time step for 1st order schemes)
                          field_value_t field_zero, // Zero-valued field_value_t
                          real_t real_zero) :       // Zero-valued real_t
        n_part(n_part),
        k(k),
        gamma_n(gamma_n),
        k_t(k_t),
        gamma_t(gamma_t),
        mu_s(mu_s),
        phi_d(phi_d),
        k_r(k_r),
        gamma_r(gamma_r),
        mu_r(mu_r),
        phi_r(phi_r),
        k_o(k_o),
        gamma_o(gamma_o),
        mu_o(mu_o),
        phi_o(phi_o),
        r_part(r_part),
        mass(mass),
        inertia(inertia),
        dt(dt),
        real_zero(real_zero),
        field_zero(field_zero) {

        contact_springs.resize(n_part * n_part);
        std::fill(contact_springs.begin(), contact_springs.end(), std::make_tuple(field_zero, field_zero, field_zero));

        bonded_contacts.resize(n_part * n_part);
        std::fill(bonded_contacts.begin(), bonded_contacts.end(), false);

        particle_to_bond_map.resize(n_part);

        for (size_t i = 0; i < n_part - 1; i ++) {
            for (size_t j = i+1; j < n_part; j ++) {
                if (abs((x0[i] - x0[j]).norm() - 2.0*r_part) < 1.0e-9) {
                    bonded_contacts[i*n_part + j] = true;
                    bonded_contacts[j*n_part + i] = true;
                    particle_to_bond_map[i].emplace_front(j);
                    particle_to_bond_map[j].emplace_front(i);
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

        if (!bonded_contacts[i*n_part + j])
            return std::make_pair(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());

        field_value_t n = (x[i] - x[j]).normalized();
        real_t overlap = 2.0 * r_part - (x[i] - x[j]).dot(n);

//        if (overlap <= 0) {
//            reset_springs(i, j); // Reset the accumulated tangential springs
//            return std::make_pair(field_zero, field_zero); // Return zeros - there is no force or torque for interparticle contact
//        }

        real_t r_part_prime = r_part - overlap / 2.0;

        real_t v_n = -(v[i] - v[j]).dot(n); // Normal relative velocity

        real_t f_n = k * overlap // Elastic contribution
                + gamma_n * v_n; // Viscous contribution

        field_value_t v_ij = v[i] - v[j] + r_part_prime * n.cross(omega[i]) + r_part_prime * n.cross(omega[j]);

        field_value_t v_t = v_ij - v_ij.dot(n) * n; // Tangential relative velocity
        field_value_t v_r = -r_part_prime / 2.0 * (n.cross(omega[i]) - n.cross(omega[j])); // Rolling velocity
        field_value_t v_o = r_part / 2.0 * (n.dot(omega[i]) - n.dot(omega[j])) * n; // Spin velocity

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
        if (xi_new.norm() > 0.0) {
            xi_new.normalize();
            xi_new *= xi.norm();
        }
        // Update the spring in the spring buffer
        xi = xi_new;
        field_value_t f_0 = -stiffness * xi - damping * relative_velocity; // Compute the test force
        xi += relative_velocity * dt;

        return f_0;
    }

    std::vector<bool> bonded_contacts;

private:
    const size_t n_part;
    const real_t k, gamma_n,
        k_t, gamma_t, mu_s, phi_d,
        k_r, gamma_r, mu_r, phi_r,
        k_o, gamma_o, mu_o, phi_o,
        r_part, mass, inertia, dt, real_zero;
    const field_value_t field_zero;
    std::vector<std::forward_list<size_t>> particle_to_bond_map;
    std::vector<std::tuple<field_value_t, field_value_t, field_value_t>> contact_springs;
};

#endif //LIBGRAN_ALT_SINTER_BRIDGE_H
