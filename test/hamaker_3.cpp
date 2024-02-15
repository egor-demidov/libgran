//
// Created by egor on 1/25/24.
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

#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

#include "../writer.h"

int main() {
    enable_fp_exceptions();

    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 3.0e-7 / 3.0;
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
    // A pentahedron
    for (size_t i = 0; i < 6; i ++) {
        double delta_z = 2.0 * r_part * 1.0 / sqrt(2.0);
        auto z = double(i) * delta_z;
        for (size_t j = 0; j < i + 1; j ++) {
            auto x = -4.0 * r_part -1.0 * r_part * double(i+1) + 2.0 * r_part * double (j);
            for (size_t m = 0; m < i + 1; m ++) {
                auto y = -1.0 * r_part * double(i+1) + 2.0 * r_part * double (m);
                x0.emplace_back(x, y, z);
                v0.emplace_back(3.0, 0.0, -5.0); // The pentahedron is moving in the -z direction
            }
        }
    }

    // An octahedron
    size_t width = 6, height = 6; // Must be multiples of two
    size_t depth = 3;
    for (size_t i = 0; i < depth; i ++) {
        auto z = -4.0 * r_part - 2.0 * r_part * double(i);
        for (size_t j = 0; j < width; j ++) {
            auto x = -double(width) * r_part + double(j) * 2.0 * r_part;
            for (size_t m = 0; m < height; m ++) {
                auto y = -double(height) * r_part + double(m) * 2.0 * r_part;
                x0.emplace_back(x, y, z);
                v0.emplace_back(Eigen::Vector3d::Zero()); // The octahedron is static
            }
        }
    }

    std::cout << x0.size() << " particles total" << std::endl;

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

    binary_force_functor_container<Eigen::Vector3d, double,
        contact_force_functor<Eigen::Vector3d, double>,
        hamaker_functor<Eigen::Vector3d, double>>
            binary_force_functors{contact_force_model, hamaker_model};

    unary_force_functor_container<Eigen::Vector3d, double>
            unary_force_functors;

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        rotational_step_handler,
        binary_force_functor_container<Eigen::Vector3d, double,
                contact_force_functor<Eigen::Vector3d, double>,
                hamaker_functor<Eigen::Vector3d, double>>,
        unary_force_functor_container<Eigen::Vector3d, double>> system(x0,
            v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
            step_handler_instance, binary_force_functors, unary_force_functors);

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
}