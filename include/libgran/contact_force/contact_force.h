//
// Created by egor on 1/23/24.
//

#ifndef LIBGRAN_CONTACT_FORCE_H
#define LIBGRAN_CONTACT_FORCE_H

template <typename field_value_t, typename real_t>
struct contact_force_functor {
    contact_force_functor(size_t n_part,            // Number of particles in the system
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
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i,
                                                         size_t j,
                                                         std::vector<field_value_t> const & x,
                                                         std::vector<field_value_t> const & v,
                                                         std::vector<field_value_t> const & theta [[maybe_unused]],
                                                         std::vector<field_value_t> const & omega [[maybe_unused]],
                                                         real_t t [[maybe_unused]]) {

        field_value_t n = (x[i] - x[j]).normalized();
        real_t overlap = 2.0 * r_part - (x[i] - x[j]).dot(n);

        if (overlap <= 0) {
            reset_springs(i, j);
            return std::make_pair(field_zero, field_zero);
        }

        real_t r_part_prime = r_part - overlap / 2.0;

        field_value_t v_ij = v[i] - v[j] + r_part_prime * n.cross(omega[i]) + r_part_prime * n.cross(omega[j]);

        real_t v_n = -v_ij.dot(n); // Normal relative velocity

        real_t f_n = k * overlap // Elastic contribution
                + gamma_n * v_n; // Viscous contribution

        field_value_t v_t = v_ij - v_n * n; // Tangential relative velocity
        field_value_t v_r = -r_part_prime / 2.0 * (n.cross(omega[i]) - n.cross(omega[j])); // Rolling velocity
        field_value_t v_o = r_part / 2.0 * (n.dot(omega[i]) - n.dot(omega[j])) * n; // Spin velocity

        field_value_t f_t = compute_shear_contribution<0>(i, j, k_t, gamma_t, f_n, mu_s, phi_d, v_t); // Sliding/sticking
        field_value_t f_r = compute_shear_contribution<1>(i, j, k_r, gamma_r, f_n, mu_r, phi_r, v_r); // Rolling
        field_value_t f_o = compute_shear_contribution<2>(i, j, k_o, gamma_o, f_n, mu_o, phi_o, v_o); // Torsion

        // Compute the torques associated with all the shear contributions
        field_value_t tau_t = r_part_prime * n.cross(f_t);
        field_value_t tau_r = r_part * n.cross(f_r);
        field_value_t tau_o = r_part * f_o;

        return std::make_pair((f_n * n + f_t) / mass, (tau_t + tau_r + tau_o) / inertia);
    }

    void reset_springs(size_t i, size_t j) {
        contact_springs[i * n_part + j] = std::make_tuple(field_zero, field_zero, field_zero);
    }

    // Computes either sliding/sticking, rolling, or torsion contribution
    // Use model_num 0 for sliding/sticking, 1 for rolling, 2 for torsion
    template<size_t model_num>
    field_value_t compute_shear_contribution(size_t i, size_t j,
                                                   real_t stiffness, real_t damping,
                                                   real_t normal_force,
                                                   real_t mu_static,
                                                   real_t phi_dynamic,
                                                   field_value_t const & relative_velocity) {

        real_t mu_dynamic = mu_static * phi_dynamic; // Compute the dynamic friction coefficient
        field_value_t & xi = std::get<model_num>(contact_springs[i * n_part + j]); // Access the respective spring
        field_value_t f_0 = -stiffness * xi - damping * relative_velocity; // Compute the test force
        real_t static_friction = mu_static * normal_force; // Compute the static friction force
        real_t dynamic_friction = mu_dynamic * normal_force; // Compute the dynamic friction force

        // Safely compute the unit tangent vector
        field_value_t t;
        if (f_0.norm() > 0.0)
            t = f_0.normalized();
        else
            t = field_zero;

        // Select whether static of dynamic friction should be used based on the test force
        real_t f_selected;
        if (f_0.norm() < static_friction) {
            // This is static friction
            xi += relative_velocity * dt;
            f_selected = static_friction;
        } else {
            // This is dynamic friction
            xi = -1.0 / stiffness * (dynamic_friction * t + damping * relative_velocity);
            f_selected = dynamic_friction;
        }

        return std::min(f_selected, f_0.norm()) * t;
    }

private:
    const size_t n_part;
    const real_t k, gamma_n,
        k_t, gamma_t, mu_s, phi_d,
        k_r, gamma_r, mu_r, phi_r,
        k_o, gamma_o, mu_o, phi_o,
        r_part, mass, inertia, dt, real_zero;
    const field_value_t field_zero;
    std::vector<std::tuple<field_value_t, field_value_t, field_value_t>> contact_springs;
};

#endif //LIBGRAN_CONTACT_FORCE_H
