//
// Created by egor on 1/25/24.
//

#include <iostream>
#include <vector>
#include <chrono>

#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

#include "../writer.h"

// The driver program is responsible for:
// 1) Instantiating the force models
// 2) Instantiating the step handler
// 3) Initializing the initial positions and velocities
// 4) Initializing sinter bridges and sinter-contact map
// 4) Instantiating granular system
// 5) Performing the time steps

int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 3.0e-7;
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

    // Parameters for the Van der Waals model
    const double A = 1.0e-20;
    const double h0 = 1.0e-9;

    // Initialize the particles
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
    // Two cubes
    size_t cube_size = 10;
    for (size_t i = 0; i < cube_size; i ++) {
        auto y = -2.0 * double(cube_size) * r_part + double(i) * 2.0 * r_part;
        for (size_t j = 0; j < cube_size; j ++) {
            auto x = -2.0 * double(cube_size) * r_part + double(j) * 2.0 * r_part;

            // First cube
            for (size_t n = 0; n < cube_size; n ++) {
                auto z = -2.0 * double(n) * r_part;
                x0.emplace_back(x, y, z);
                v0.emplace_back(0.0, 0.0, 5.0); // Cube moving up
            }

            // Second cube
            for (size_t n = 0; n < cube_size; n ++) {
                auto z = 8.0 * r_part + 2.0 * double(n) * r_part;
                x0.emplace_back(x + r_part * 0.5 * double(cube_size), y + r_part * 0.5 * double(cube_size), z);
                v0.emplace_back(0.0, 0.0, -5.0); // Cube moving down
            }
        }
    }

    std::cout << "Number of particles: " << x0.size() << std::endl;

    // Initialize the remaining buffers
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    // Set the initial velocity of the above-plane particle

    // Create an instance of step_handler
    // Using field type Eigen::Vector3d with container std::vector
    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;

    // Create an instance of contact force model
    // Using field type Eigen::Vector3d with real type double
    contact_force_functor<Eigen::Vector3d, double> contact_force_model(x0.size(),
                                                                       k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                                                       r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of Hamaker model
    hamaker_functor<Eigen::Vector3d, double> hamaker_model(A, h0,
                                                           r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
            rotational_step_handler, contact_force_functor<Eigen::Vector3d, double>,
            hamaker_functor<Eigen::Vector3d, double>> system(x0,
                                                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                                                             step_handler_instance, contact_force_model, hamaker_model);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump #" << n / dump_period << std::endl;
            write_particles("run", system.get_x(), system.get_theta(), r_part);
        }
        system.do_step(dt);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "Elapsed time: " << duration << " sec" << std::endl;

    return 0;
}