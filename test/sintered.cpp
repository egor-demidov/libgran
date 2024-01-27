//
// Created by egor on 1/25/24.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <random>

#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/granular_system/granular_system.h>
#include <libgran/sinter_bridge/sinter_bridge.h>

#include "mass_distribution.h"
#include "../writer.h"

std::mt19937_64 mt(0);
std::uniform_real_distribution<double> dist(-1.0, 1.0);

Eigen::Vector3d generate_random_unit_vector() {
    Eigen::Vector3d vec;
    do {
        vec = {dist(mt), dist(mt), dist(mt)};
    } while (vec.norm() == 0);
    return vec.normalized();
}

// We will be using a custom step handler in this simulation - sinter_step_handler
// which has three template arguments: field_container_t, field_value_t, and real_t
// but granular_system expects the step handler to have only two template arguments:
// field_container_t and field_value_t
// Therefore, we need to create an alias where real_t is specialized
template <typename field_container_t, typename field_value_t>
using sinter_step_handler_double = sinter_step_handler<field_container_t, field_value_t, double>;

int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 1.0e-7;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 300;
    const size_t dump_period = n_steps / n_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

    // Parameters for the contact model
    const double k = 10000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

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

    sinter_functor<Eigen::Vector3d, double> sinter_model(x0.size(), k, gamma_n, k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                                         r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0, x0.begin(), generate_random_unit_vector, 1.0e-9);

    // We need to use a custom step handler with the sintering model
    // Get the custom step handler instance from the sinter model
    auto step_handler_instance = sinter_model.get_step_handler<std::vector<Eigen::Vector3d>>();

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        sinter_step_handler_double, sinter_functor<Eigen::Vector3d, double>> system(x0,
                                                                                    v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0, step_handler_instance,
                                                                                    sinter_model);

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

    double target_rg = 3.95981e-08;
    double target_center_dist = 2.35702e-07;

    double rg_tolerance = 1.0; // Percent
    double center_dist_tolerance = 1.0; // Percent

    if (abs(center_dist - target_center_dist) / target_center_dist > center_dist_tolerance / 100.0 ||
        abs(rg - target_rg) / target_rg > rg_tolerance / 100.0)
        return EXIT_FAILURE;

    return 0;
}