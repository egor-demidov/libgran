//
// Created by egor on 1/23/24.
// Editted by gurdeep on 11/29/25.
//

#ifndef LIBGRAN_CONTACT_FORCE_H
#define LIBGRAN_CONTACT_FORCE_H

template <typename field_value_t>
void update_particle_pressures(std::vector<std::array<double, 9>> & p, field_value_t force, field_value_t rIJ, int i){
    for(int k = 0; k < 3; k++){
        p[i][k] += force[k] * rIJ[k];
    }
    p[i][3] += force[1] * rIJ[0];
    p[i][4] += force[2] * rIJ[0];
    p[i][5] += force[2] * rIJ[1];

    p[i][6] += force[0] * rIJ[1];
    p[i][7] += force[0] * rIJ[2];
    p[i][8] += force[1] * rIJ[2];
}

template <typename field_value_t, typename real_t, typename matrix_t>
struct contact_force_functor_var_size {
    contact_force_functor_var_size(size_t n_part,   // Number of particles in the system
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
                          real_t RH,                // relative humidity
                          real_t dt,                // Time step for spring update (same as integration time step for 1st order schemes)
                          field_value_t field_zero, // Zero-valued field_value_t
                          real_t real_zero          // Zero-valued real_t
                          ) :       
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
        RH(RH),
        dt(dt),
        real_zero(real_zero),
        field_zero(field_zero) {

        contact_springs.resize(n_part * n_part);
        std::fill(contact_springs.begin(), contact_springs.end(), std::make_tuple(field_zero, field_zero, field_zero));

        contact_active.resize(n_part * n_part, false);

        a_prev.resize(n_part * n_part);
        std::fill(a_prev.begin(), a_prev.end(), 0.0);       
    }

    template<typename BoxType>
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
        size_t idx = i * n_part + j;

        // minimum image convention
        field_value_t d_raw = x[i] - x[j];
        field_value_t d = box.minimumImage(d_raw);

        field_value_t n = d.normalized();
        real_t overlap = (r[i] + r[j]) - d.dot(n);

        // if (overlap <= 0) [[likely]] {
        //     reset_springs(i, j); // Reset the accumulated tangential springs
        //     // return std::make_pair(field_zero, field_zero); // Return zeros - there is no force or torque for interparticle contact
        // }

        // real_t r_part_prime = (r[i] + r[j])/2.0 - overlap/2.0;
        real_t r_i_prime = r[i] - 1.0/2.0 * overlap;
        real_t r_j_prime = r[j] - 1.0/2.0 * overlap;
        real_t r_ij_prime = r_i_prime * r_j_prime / (r_i_prime + r_j_prime);

        const matrix_t& D = box.get_deformation_rate();
        // Subtract affine streaming velocity
        field_value_t uij = (v[i] - v[j]) - D * d_raw;

        real_t v_n = -uij.dot(n); // Normal relative velocity

        // JKR model
        real_t Reff = (r[i] * r[j]) / (r[i] + r[j]);

        real_t E_grain = 88.7e9; // Pa
        real_t nu = 0.166;
        real_t E_star = E_grain / (2.0 * (1.0 - nu * nu));

        real_t gamma_rh = calculate_gamma(this->RH);

        real_t delta_to = std::pow( (3.0 * M_PI * M_PI * gamma_rh * gamma_rh * Reff) / (16.0 * E_star * E_star), 1.0/3.0 );

        if (!contact_active[idx]) {
            if (overlap > 0) {
                contact_active[idx] = true; 
                a_prev[idx] = std::pow((2.0 * M_PI * gamma_rh * Reff * Reff / E_star), 1.0/3.0);
            } else { 
                contact_radius = 0;
                normal_force = 0;
                normalized_overlap = 0;
                return std::make_pair(field_zero, field_zero);
            }
        } else {
            if (overlap < -delta_to) {
                contact_active[idx] = false; 
                reset_springs(i, j);
                a_prev[idx] = 0.0;
                return std::make_pair(field_zero, field_zero);
            }
        }

        real_t a = a_prev[idx];
        real_t a0 = std::pow((4.5 * M_PI * gamma_rh * Reff * Reff) / E_star, 1.0/3.0);
        if (a < 1e-15) a = a0;

        for (int iter = 0; iter < 10; ++iter) {
            real_t term = std::sqrt(2.0 * M_PI * gamma_rh * a / E_star); 
            real_t f = (a * a / Reff) - term - overlap;
            real_t df = (2.0 * a / Reff) - (0.5 * term / a);

            if (std::abs(f) < 1e-12) break;
            
            a -= f / df;
            if (a < 0.3 * a0) { a = 0.3 * a0; break; }
        }
        a_prev[idx] = a;

        real_t f_n_elastic = (4.0 * E_star * std::pow(a, 3.0) / (3.0 * Reff)) 
                   - std::sqrt(8.0 * M_PI * gamma_rh * E_star * std::pow(a, 3.0));

        // remove damping for no hysterisis
        // real_t f_n = f_n_elastic;// + gamma_n * v_n;
        real_t f_n = f_n_elastic + gamma_n * v_n;
        if (i == 0 && j == 1) {
            real_t F_po = 1.5 * M_PI * gamma_rh * Reff;
            this->normal_force = f_n/F_po;
            this->contact_radius = a/a0;
            this->normalized_overlap = overlap/delta_to;
        }

        // Add rotational contributions
        field_value_t v_ij = uij + r_i_prime * n.cross(omega[i]) + r_j_prime * n.cross(omega[j]);

        field_value_t v_t = v_ij - v_ij.dot(n) * n; // Tangential relative velocity
        field_value_t v_r = r_ij_prime * (-n.cross(omega[i]) + n.cross(omega[j])); // Rolling velocity
        field_value_t v_o = r_ij_prime * (n.dot(omega[i]) - n.dot(omega[j])) * n; // Spin velocity

        field_value_t f_t = compute_shear_contribution<0>(i, j, n, k_t, gamma_t, f_n, mu_s, phi_d, v_t); // Sliding/sticking
        field_value_t f_r = compute_shear_contribution<1>(i, j, n, k_r, gamma_r, f_n, mu_r, phi_r, v_r); // Rolling
        field_value_t f_o = compute_shear_contribution<2>(i, j, n, k_o, gamma_o, f_n, mu_o, phi_o, v_o); // Torsion

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

    std::array<real_t, 3> get_normal_force(){
        return {normal_force, contact_radius, normalized_overlap};
    }

private:
    void reset_springs(size_t i, size_t j) {
        contact_springs[i * n_part + j] = std::make_tuple(field_zero, field_zero, field_zero);
    }

    // Computes either sliding/sticking, rolling, or torsion contribution
    // Use model_num 0 for sliding/sticking, 1 for rolling, 2 for torsion
    template<size_t model_num>
    field_value_t compute_shear_contribution(size_t i, size_t j, field_value_t const & n,
                                                   real_t stiffness, real_t damping,
                                                   real_t normal_force,
                                                   real_t mu_static,
                                                   real_t phi_dynamic,
                                                   field_value_t const & relative_velocity) {

        real_t mu_dynamic = mu_static * phi_dynamic; // Compute the dynamic friction coefficient
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
        real_t static_friction = mu_static * normal_force; // Compute the static friction force
        real_t dynamic_friction = mu_dynamic * normal_force; // Compute the dynamic friction force

        // Safely compute the unit tangent vector
        field_value_t t;
        if (f_0.norm() > 0.0) [[likely]]
            t = f_0.normalized();
        else [[unlikely]]
            t = field_zero;

        // Select whether static of dynamic friction should be used based on the test force
        real_t f_selected;
        if (f_0.norm() <= static_friction) [[likely]] {
            // This is static friction
            xi += relative_velocity * dt;
            f_selected = static_friction;
        } else [[unlikely]] {
            // This is dynamic friction
            xi = -1.0 / stiffness * (dynamic_friction * t + damping * relative_velocity);
            f_selected = dynamic_friction;
        }

        return std::min(f_selected, f_0.norm()) * t;
    }

    real_t calculate_gamma(real_t current_RH) {
        const real_t gamma_s = 0.45;       // J/m^2
        const real_t RH_c = 0.70;          // 70%
        const real_t Rg = 8.314;           
        const real_t T = 298.0;            
        const real_t Gamma_m = 3.7e-5;     // mol/m^2 (from 0.037 mmol/m^2)
        const real_t C = 21.3;

        if (current_RH <= 0.0) return gamma_s;
        
        // Clamp RH for the plateau after capillary condensation
        real_t effective_RH = std::min(current_RH, RH_c);
        
        // Equation 5: Adsorption-induced surface energy reduction
        real_t reduction = Rg * T * Gamma_m * (std::log(1.0 + (C - 1.0) * effective_RH) - std::log(1.0 - effective_RH));
        return gamma_s - reduction;
    }

    const size_t n_part;
    const real_t k, gamma_n,
        k_t, gamma_t, mu_s, phi_d,
        k_r, gamma_r, mu_r, phi_r,
        k_o, gamma_o, mu_o, phi_o, RH,
        dt, real_zero;
    const field_value_t field_zero;
    std::vector<std::tuple<field_value_t, field_value_t, field_value_t>> contact_springs;
    std::vector<bool> contact_active;
    real_t normal_force, contact_radius, normalized_overlap;
    std::vector<real_t> a_prev;
};

#endif //LIBGRAN_CONTACT_FORCE_H
