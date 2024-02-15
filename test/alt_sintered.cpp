//
// Created by egor on 2/2/24.
//

#include <iostream>
#include <vector>
#include <chrono>

#ifdef _GNU_SOURCE
#include <cfenv>
#define enable_fp_exceptions() feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO)
#elif defined(_MSC_VER)
#define enable_fp_exceptions()
#elif defined(__APPLE__)
#define enable_fp_exceptions()
#else
#error "Unsupported system"
#endif

#include <Eigen/Eigen>

#include <libgran/sinter_bridge/alt_sinter_bridge.h>
#include <libgran/granular_system/granular_system.h>

#include "mass_distribution.h"
#include "../writer.h"

using contact_functor_t = contact_force_functor<Eigen::Vector3d, double>;
using alt_sinter_functor_t = alt_sinter_functor<Eigen::Vector3d, double>;
using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, alt_sinter_functor_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;
using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half, rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main () {
    enable_fp_exceptions();

    // General simulation parameters
    const double dt = 5e-14;
    const double t_tot = 1.0e-7;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 300;
    const size_t dump_period = n_steps / n_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

    // Parameters for the non-bonded contact model
    const double k = 10000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the bonded contact model
    const double k_bond = 1000.0;
    const double gamma_n_bond = 2.0*sqrt(2.0*mass*k_bond);
    const double gamma_t_bond = 0.2 * gamma_n_bond;
    const double gamma_r_bond = 0.05 * gamma_n_bond;
    const double gamma_o_bond = 0.05 * gamma_n_bond;
    const double d_crit = 1.0e-9; // Critical separation

    // Initialize the particles
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
    for (size_t i = 0; i < 3; i ++) {
        auto y = -1.0 * r_part + double(i) * 2.0 * r_part;
        for (size_t j = 0; j < 3; j ++) {
            auto x = -1.0 * r_part + double(j) * 2.0 * r_part;
            x0.emplace_back(x, y, 0.0);
            x0.emplace_back(x, y, -2.0*r_part);
            x0.emplace_back(x, y, -4.0*r_part);

            // Particles on the x=x_min edge will have an initial velocity in the +z direction
            // Particle on the x=x_max edge will have an initial velocity in the +y direction
            Eigen::Vector3d init_vel;
            if (j == 0)
                init_vel = {0.0, 0.0, 5.0};
            else if (j == 2)
                init_vel = {0.0, 5.0, 0.0};
            else
                init_vel = Eigen::Vector3d::Zero();
            for (size_t n = 0; n < 3; n ++)
                v0.emplace_back(init_vel);
        }
    }

    // Initialize the remaining buffers
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    // Create an instance of step_handler
    // Using field type Eigen::Vector3d with container std::vector
    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;

    contact_functor_t  contact_force(x0.size(),
                                     k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                     r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of contact force model
    // Using field type Eigen::Vector3d with real type double
    alt_sinter_functor_t sinter_model(x0.size(), x0,
                               k_bond, gamma_n_bond, k_bond, gamma_t_bond, k_bond, gamma_r_bond, k_bond, gamma_o_bond,
                               r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0, d_crit, contact_force);

    binary_force_container_t
        binary_force_functors{sinter_model};

    unary_force_container_t
            unary_force_functors;

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system_t system(x0,
                   v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0, step_handler_instance, binary_force_functors, unary_force_functors);

    Eigen::Vector3d center_0 = center_of_mass(system.get_x());

    auto start_time = std::chrono::high_resolution_clock::now();
    for (size_t n = 0; n < n_steps; n ++) {
        system.do_step(dt);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "Elapsed time: " << duration << " sec" << std::endl;

    Eigen::Vector3d center_f = center_of_mass(system.get_x());

    double rg = radius_of_gyration(system.get_x());
    double center_dist = (center_f - center_0).norm();

    std::cout << rg << " " << center_dist << std::endl;

    double target_rg = 3.9598e-08;
    double target_center_dist = 2.35702e-07;

    double rg_tolerance = 1.0; // Percent
    double center_dist_tolerance = 1.0; // Percent

    if (abs(center_dist - target_center_dist) / target_center_dist > center_dist_tolerance / 100.0 ||
        abs(rg - target_rg) / target_rg > rg_tolerance / 100.0)
        return EXIT_FAILURE;

    return 0;
}
