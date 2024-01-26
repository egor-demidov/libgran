//
// Created by egor on 1/24/24.
//

#ifndef LIBGRAN_SINTER_BRIDGE_H
#define LIBGRAN_SINTER_BRIDGE_H

#include <functional>

template <typename field_container_t, typename field_value_t, typename real_t>
struct sinter_step_handler {

    // Do not initialize directly, use get_step_handler() of sinter_functor instead
    sinter_step_handler(typename std::vector<std::array<field_value_t, 5>>::iterator begin_spring_connectors_itr,
                        std::vector<bool>::const_iterator begin_enabled_contacts, real_t r_part, size_t n_part) :
            begin_spring_connectors_itr(begin_spring_connectors_itr), begin_enabled_contacts(begin_enabled_contacts), r_part(r_part), n_part(n_part) {}

    void increment_x(size_t n,
                     field_value_t const & dx,
                     typename field_container_t::iterator x_begin_itr,
                     typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator omega_begin_itr [[maybe_unused]]) const {
        *(x_begin_itr + n) += dx;
        // TODO: create an index with spring connectors attached to each particle to improve performance
        // Increment spring connectors for bonded contacts
        for (size_t m = n * n_part; m < (n + 1) * n_part; m ++) {
            if (!*(begin_enabled_contacts + m)) {
                // This is a bonded contact. Do increments
                for (size_t k = 0; k < 5; k ++) {
                    (*(begin_spring_connectors_itr + m))[k] += dx;
                }
            }
        }
    }

    void increment_v(size_t n,
                     field_value_t const & dv,
                     typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator v_begin_itr,
                     typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator omega_begin_itr [[maybe_unused]]) const {
        *(v_begin_itr + n) += dv;
    }

    void increment_theta(size_t n,
                         field_value_t const & dtheta,
                         typename field_container_t::iterator x_begin_itr,
                         typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator theta_begin_itr,
                         typename field_container_t::iterator omega_begin_itr [[maybe_unused]]) const {
        *(theta_begin_itr + n) += dtheta;

        // Increment spring connectors for bonded contacts
        // TODO: create an index with spring connectors attached to each particle to improve performance
        for (size_t m = n * n_part; m < (n + 1) * n_part; m ++) {
            if (!*(begin_enabled_contacts + m)) {
                // This is a bonded contact. Do increments
                for (size_t k = 0; k < 5; k ++) {
                    // Compute the spring vector
                    field_value_t spring_vector = (*(begin_spring_connectors_itr + m))[k] - *(x_begin_itr + n);
                    (*(begin_spring_connectors_itr + m))[k] += dtheta.cross(spring_vector);
                }
            }
        }
    }

    void increment_omega(size_t n,
                         field_value_t const & domega,
                         typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator omega_begin_itr) const {
        *(omega_begin_itr + n) += domega;
    }

private:
    typename std::vector<std::array<field_value_t, 5>>::iterator begin_spring_connectors_itr;
    std::vector<bool>::const_iterator begin_enabled_contacts;
    const real_t r_part;
    const size_t n_part;
};

template <typename field_value_t, typename real_t>
struct sinter_functor {
    sinter_functor(size_t n_part,
                   real_t k,
                   real_t gamma_d,
                   real_t r_part,
                   real_t mass,
                   real_t inertia,
                   field_value_t field_zero,
                   real_t real_zero,
                   real_t critical_separation,
                   typename std::vector<field_value_t>::const_iterator x_begin_itr,
                   typename std::vector<field_value_t>::const_iterator x_end_itr,
                   std::vector<bool> & enabled_contacts,
                   std::function<field_value_t()> random_unit_vector_generator) :
                       n_part(n_part),
                       k(k), gamma_d(gamma_d), r_part(r_part), mass(mass),
                       inertia(inertia), real_zero(real_zero), field_zero(field_zero),
                       enabled_contacts(enabled_contacts), random_unit_vector_generator(random_unit_vector_generator) {

        spring_connectors.resize(n_part * n_part);
        spring_lengths.resize(n_part * n_part);

        // By default, sinter bridges will be inserted between all neighboring particles
        // Iterate over all particles
        for (size_t i = 0; i < n_part - 1; i ++) {
            for (size_t j = i + 1; j < n_part; j ++) {
                // Create references for convenience
                field_value_t const & part_i = *(x_begin_itr + i);
                field_value_t const & part_j = *(x_begin_itr + j);

                // Compute the separation between particles
                real_t separation = (part_j - part_i).norm() - 2.0 * r_part;
                if (abs(separation) > critical_separation)
                    continue; // Particles are too far apart to be considered neighbors

                // Mark this pair as necked
                enabled_contacts[i * n_part + j] = false;
                enabled_contacts[j * n_part + i] = false;

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
            }
        }
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i,
                                                         size_t j,
                                                         std::vector<field_value_t> const & x,
                                                         std::vector<field_value_t> const & v [[maybe_unused]],
                                                         std::vector<field_value_t> const & theta [[maybe_unused]],
                                                         std::vector<field_value_t> const & omega [[maybe_unused]],
                                                         real_t t [[maybe_unused]]) {

        if (enabled_contacts[i * n_part + j])
            return std::make_pair(field_zero, field_zero); // No acceleration if contact is enabled

        field_value_t force = field_zero, torque = field_zero;

        auto const & part_i_connectors = spring_connectors[i * n_part + j];
        auto const & part_j_connectors = spring_connectors[j * n_part + i];

        // Iterate over the attached springs
        for (size_t n = 0; n < 5; n ++) {
            field_value_t spring = part_j_connectors[n] - part_i_connectors[n];
            real_t spring_length = spring.norm();
            real_t equilibrium_length = spring_lengths[i * n_part + j][n];
            field_value_t dforce = -(equilibrium_length - spring_length) * k * spring.normalized();
            force += dforce;

            field_value_t arm = part_i_connectors[n] - x[i];
            torque += arm.cross(dforce);
        }

//        if (i == 0 && force.norm() > 0) {
//            std::cout << force.normalized() << std::endl;
//            exit(0);
//        }

        return std::make_pair(force / mass, torque / inertia);
    }

    // Create an instance of sinter_step_handler to be used with granular_system
    template <typename field_container_t>
    sinter_step_handler<field_container_t, field_value_t, real_t> get_step_handler() {
        return sinter_step_handler<field_container_t, field_value_t, real_t> (spring_connectors.begin(), enabled_contacts.begin(), r_part, n_part);
    }

private:


    std::vector<std::array<field_value_t, 5>> spring_connectors; // Positions of spring connectors
    std::vector<std::array<real_t, 5>> spring_lengths; // Equilibrium lengths of springs

    const size_t n_part;
    const real_t k, gamma_d, r_part, mass, inertia, real_zero;
    const field_value_t field_zero;
    std::vector<bool> & enabled_contacts; // Reference to enabled_contacts from contact_force
    std::function<field_value_t()> random_unit_vector_generator;

    friend int main();
};

#endif //LIBGRAN_SINTER_BRIDGE_H
