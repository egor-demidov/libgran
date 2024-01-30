//
// Created by egor on 1/30/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/granular_system/granular_system.h>
#include <libgran/sinter_bridge/sinter_bridge.h>

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

// "Assemble" the force models used in this simulation
using sinter_functor_t = sinter_functor<Eigen::Vector3d, double>; // Contact force
using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, sinter_functor_t>; // Binary force container

using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>; // Unary force container (empty)

// We will be using a custom step handler in this simulation - sinter_step_handler
// which has three template arguments: field_container_t, field_value_t, and real_t
// but granular_system expects the step handler to have only two template arguments:
// field_container_t and field_value_t
// Therefore, we need to create an alias where real_t is specialized
template <typename field_container_t, typename field_value_t>
using sinter_step_handler_double = sinter_step_handler<field_container_t, field_value_t, double>;

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        sinter_step_handler_double, binary_force_container_t, unary_force_container_t>; // Granular system representation

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
    x0.emplace_back(0, 0, 0);
    x0.emplace_back(0, 2.0*r_part, 0);
    x0.emplace_back(2.0*r_part, 2.0*r_part, 0);
    v0.emplace_back(1, 0, 0);
    v0.emplace_back(0, 0, 0);
    v0.emplace_back(-1, 0, 0);

    theta0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    omega0.resize(x0.size());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    sinter_functor_t sinter_model(x0.size(), k, 0, k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                                         r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0, x0.begin(), generate_random_unit_vector, 1.0e-9);

    binary_force_container_t
            binary_force_functors{sinter_model};

    unary_force_container_t
            unary_force_functors;

    // We need to use a custom step handler with the sintering model
    // Get the custom step handler instance from the sinter model
    auto step_handler_instance = sinter_model.get_step_handler<std::vector<Eigen::Vector3d>>();

    granular_system_t system(x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(),
                             0.0, step_handler_instance,
                             binary_force_functors, unary_force_functors);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump # " << n / dump_period << std::endl;
            write_particles("run", system.get_x(), system.get_omega(), r_part);
        }
        system.do_step(dt);
    }

    return 0;
}
