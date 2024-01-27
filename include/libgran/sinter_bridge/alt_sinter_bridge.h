//
// Created by egor on 1/26/24.
//

#ifndef LIBGRAN_ALT_SINTER_BRIDGE_H
#define LIBGRAN_ALT_SINTER_BRIDGE_H

#include <forward_list>
#include <functional>

#include "../contact_force/contact_force.h"

template <typename field_container_t, typename field_value_t, typename real_t>
struct alt_sinter_step_handler {

    // Do not initialize directly, use get_step_handler() of sinter_functor instead
    alt_sinter_step_handler(typename std::vector<std::array<field_value_t, 5>>::iterator begin_spring_connectors_itr,
                            std::vector<bool>::const_iterator begin_bonded_contacts,
                            std::vector<std::forward_list<size_t>>::const_iterator begin_particle_to_contact_map,
                            real_t r_part, size_t n_part) :
            begin_spring_connectors_itr(begin_spring_connectors_itr), begin_bonded_contacts(begin_bonded_contacts),
            begin_particle_to_contact_map(begin_particle_to_contact_map), r_part(r_part), n_part(n_part) {}

    void increment_x(size_t n,
                     field_value_t const & dx,
                     typename field_container_t::iterator x_begin_itr,
                     typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator a_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator omega_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator alpha_begin_itr [[maybe_unused]]) const {
        *(x_begin_itr + n) += dx;

        std::forward_list<size_t> const & contact_index = *(begin_particle_to_contact_map + n);
        // Increment spring connectors for bonded contacts
        for (size_t index : contact_index) {
            // This is a bonded contact. Do increments
            for (size_t k = 0; k < 5; k ++) {
                (*(begin_spring_connectors_itr + index))[k] += dx;
            }
        }
    }

    void increment_v(size_t n,
                     field_value_t const & dv,
                     typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator v_begin_itr,
                     typename field_container_t::iterator a_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator omega_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator alpha_begin_itr [[maybe_unused]]) const {
        *(v_begin_itr + n) += dv;
    }

    void increment_theta(size_t n,
                         field_value_t const & dtheta,
                         typename field_container_t::iterator x_begin_itr,
                         typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator a_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator theta_begin_itr,
                         typename field_container_t::iterator omega_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator alpha_begin_itr [[maybe_unused]]) const {
        *(theta_begin_itr + n) += dtheta;

        std::forward_list<size_t> const & contact_index = *(begin_particle_to_contact_map + n);
        // Increment spring connectors for bonded contacts
        for (size_t index : contact_index) {
            // This is a bonded contact. Do increments
            for (size_t k = 0; k < 5; k ++) {
                // Compute the spring vector
                field_value_t spring_vector = (*(begin_spring_connectors_itr + index))[k] - *(x_begin_itr + n);
                (*(begin_spring_connectors_itr + index))[k] += dtheta.cross(spring_vector);
            }
        }
    }

    void increment_omega(size_t n,
                         field_value_t const & domega,
                         typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator a_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator omega_begin_itr,
                         typename field_container_t::iterator alpha_begin_itr [[maybe_unused]]) const {
        *(omega_begin_itr + n) += domega;
    }

private:
    typename std::vector<std::array<field_value_t, 5>>::iterator begin_spring_connectors_itr;
    std::vector<bool>::const_iterator begin_bonded_contacts;
    std::vector<std::forward_list<size_t>>::const_iterator begin_particle_to_contact_map;
    const real_t r_part;
    const size_t n_part;
};

// Alternative sinter bridge implementation
// This force model is combined with contact_force
// No separate contact model is needed
template <typename field_value_t, typename real_t>
struct alt_sinter_functor {
    alt_sinter_functor(size_t n_part,
                       real_t k_d,                  // Stiffness of the sinter springs
                       real_t gamma_d,              // Damping coefficient of the sinter springs
                       real_t k,
                       real_t gamma_n,
                       real_t k_t,
                       real_t gamma_t,
                       real_t mu_s,
                       real_t phi_d,
                       real_t k_r,
                       real_t gamma_r,
                       real_t mu_r,
                       real_t phi_r,
                       real_t k_o,
                       real_t gamma_o,
                       real_t mu_o,
                       real_t phi_o,
                       real_t r_part,
                       real_t mass,
                       real_t inertia,
                       real_t dt,
                       field_value_t field_zero,
                       real_t real_zero,
                       typename std::vector<field_value_t>::const_iterator x0_begin_itr,
                       std::function<field_value_t()> random_unit_vector_generator,
                       real_t critical_separation) :
       n_part(n_part),
       k_d(k_d),
       gamma_d(gamma_d),
       r_part(r_part),
       mass(mass),
       inertia(inertia),
       real_zero(real_zero),
       field_zero(field_zero),
       contact_force(n_part, k, gamma_n, k_t, gamma_t, mu_s, phi_d, k_r, gamma_r, mu_r, phi_r,
                     k_o, gamma_o, mu_o, phi_o, r_part, mass, inertia, dt, field_zero, real_zero) {

        particle_to_connector_map.resize(n_part);
        spring_connectors.resize(n_part * n_part);
        spring_lengths.resize(n_part * n_part);
        bonded_contacts.resize(n_part * n_part);

        // All contacts are non-bonded by default
        std::fill(bonded_contacts.begin(), bonded_contacts.end(), false);

        // By default, sinter bridges will be inserted between all neighboring particles
        // Iterate over all particles
        for (size_t i = 0; i < n_part - 1; i ++) {
            for (size_t j = i + 1; j < n_part; j ++) {
                // Create references for convenience
                field_value_t const & part_i = *(x0_begin_itr + i);
                field_value_t const & part_j = *(x0_begin_itr + j);

                // Compute the separation between particles
                real_t separation = (part_j - part_i).norm() - 2.0 * r_part;
                if (abs(separation) > critical_separation)
                    continue; // Particles are too far apart to be considered neighbors

                // Mark this pair as necked
                bonded_contacts[i * n_part + j] = true;
                bonded_contacts[j * n_part + i] = true;

                // Normal unit vectors
                field_value_t n_i = (part_j - part_i).normalized();
                field_value_t n_j = - n_i;

                field_value_t v; // test vector
                field_value_t db; // binormal vector

                // Generator a test unit vector that is not collinear with n_i/n_j
                do {
                    v = random_unit_vector_generator();
                    // db is a vector in the contact plane orthogonal to (future) vector ds
                    db = v.cross(n_j);
                } while (db.norm() == 0.0);

                // ds is the projection of the test vector onto the contact plane
                field_value_t ds = v - v.dot(n_i) * n_i;

                // Normalize the ds and db vectors
                db.normalize();
                ds.normalize();

                // Spring connectors attached to particle i
                auto & connector_part_i = spring_connectors[i * n_part + j];
                connector_part_i[0] = part_i + ds * r_part;
                connector_part_i[1] = part_i - ds * r_part;
                connector_part_i[2] = part_i + db * r_part;
                connector_part_i[3] = part_i - db * r_part;
                connector_part_i[4] = part_i + n_i * r_part;

                // Spring connectors attached to particle j
                auto & connector_part_j = spring_connectors[j * n_part + i];
                connector_part_j[0] = part_j + ds * r_part;
                connector_part_j[1] = part_j - ds * r_part;
                connector_part_j[2] = part_j + db * r_part;
                connector_part_j[3] = part_j - db * r_part;
                connector_part_j[4] = part_j + n_j * r_part;

                // Set spring lengths
                for (size_t n = 0; n < 4; n ++) {
                    spring_lengths[i * n_part + j][n] = 2.0 * r_part;
                    spring_lengths[j * n_part + i][n] = 2.0 * r_part;
                }
                spring_lengths[i * n_part + j][4] = real_zero;
                spring_lengths[j * n_part + i][4] = real_zero;

                particle_to_connector_map[i].push_front(i * n_part + j);
                particle_to_connector_map[j].push_front(j * n_part + i);
            }
        }
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i,
                                                         size_t j,
                                                         std::vector<field_value_t> const & x,
                                                         std::vector<field_value_t> const & v,
                                                         std::vector<field_value_t> const & theta,
                                                         std::vector<field_value_t> const & omega,
                                                         real_t t) {

        if (bonded_contacts[i * n_part + j])
            return bonded_contact_accelerations(i, j, x);
        return contact_force(i, j, x, v, theta, omega, t);
    }

    // Create an instance of sinter_step_handler to be used with granular_system
    template <typename field_container_t>
    alt_sinter_step_handler<field_container_t, field_value_t, real_t> get_step_handler() {
        return alt_sinter_step_handler<field_container_t, field_value_t, real_t> (spring_connectors.begin(),
                          bonded_contacts.begin(), particle_to_connector_map.begin(), r_part, n_part);
    }

private:
    std::pair<field_value_t, field_value_t> bonded_contact_accelerations(size_t i, size_t j, std::vector<field_value_t> const & x) {
        field_value_t force = field_zero, torque = field_zero;

        auto const & part_i_connectors = spring_connectors[i * n_part + j];
        auto const & part_j_connectors = spring_connectors[j * n_part + i];

        // Iterate over the attached springs
        for (size_t n = 0; n < 5; n ++) {
            field_value_t spring = part_j_connectors[n] - part_i_connectors[n];
            real_t spring_length = spring.norm();
            real_t equilibrium_length = spring_lengths[i * n_part + j][n];
            field_value_t dforce = -(equilibrium_length - spring_length) * k_d * spring.normalized();
            force += dforce;

            field_value_t arm = part_i_connectors[n] - x[i];
            torque += arm.cross(dforce);
        }

        return std::make_pair(force / mass, torque / inertia);
    }

    std::vector<std::forward_list<size_t>> particle_to_connector_map; // not strictly needed, but will speed up the step handler
    std::vector<std::array<field_value_t, 5>> spring_connectors;
    std::vector<std::array<real_t, 5>> spring_lengths;
    std::vector<bool> bonded_contacts;

    const size_t n_part;
    const real_t k_d, gamma_d, r_part, mass, inertia, real_zero;
    const field_value_t field_zero;
    contact_force_functor<field_value_t, real_t> contact_force;

    friend int main(); // TODO: remove after debugging
};

#endif //LIBGRAN_ALT_SINTER_BRIDGE_H
