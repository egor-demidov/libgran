//
// Created by egor on 2/2/24.
//

#include <iostream>
#include <vector>

#ifdef _GNU_SOURCE
#include <cfenv>
#define enable_fp_exceptions() feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO)
#elif defined(_MSC_VER)
#define enable_fp_exceptions()
#endif

#include <Eigen/Eigen>

#include <libgran/sinter_bridge/alt_sinter_bridge.h>
#include <libgran/granular_system/granular_system.h>

#include "../writer.h"

using sinter_functor_t = alt_sinter_functor<Eigen::Vector3d, double>;
using alt_sinter_functor_t = alt_sinter_functor<Eigen::Vector3d, double>;
using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, alt_sinter_functor_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;
using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half, rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main () {
    enable_fp_exceptions();

    // General simulation parameters
    const double dt = 1e-14;
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

    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0.emplace_back(0, 0, 0);
    x0.emplace_back(0, 2.0*r_part, 0);
    x0.emplace_back(2.0*r_part, 2.0*r_part, 0);
    x0.emplace_back(2.0*r_part, 4.0*r_part, 0);
    x0.emplace_back(2.0*r_part, 6.0*r_part, 0);
    x0.emplace_back(2.0*r_part, 8.0*r_part, 0);

    v0.emplace_back(1, 0, 0);
    v0.emplace_back(0, 0, 0);
    v0.emplace_back(0, 0, 0);
    v0.emplace_back(0, 0, 0);
    v0.emplace_back(0, 0, 0);
    v0.emplace_back(0, 0, 0);

//    omega0.emplace_back(0, 100000000, 0);
    omega0.emplace_back(0, 0, 0);
    omega0.emplace_back(0, 0, 0);
    omega0.emplace_back(0, 0, 0);
    omega0.emplace_back(0, 0, 0);
    omega0.emplace_back(0, 0, 0);
    omega0.emplace_back(0, 0, 0);

    theta0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());

    // Create an instance of step_handler
    // Using field type Eigen::Vector3d with container std::vector
    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;

    // Create an instance of contact force model
    // Using field type Eigen::Vector3d with real type double
    sinter_functor_t sinter_model(x0.size(), x0,
                               k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                               r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    binary_force_container_t
        binary_force_functors{sinter_model};

    unary_force_container_t
            unary_force_functors;

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system_t system(x0,
                   v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0, step_handler_instance, binary_force_functors, unary_force_functors);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            write_particles("run", system.get_x(), system.get_theta(), r_part);
            std::cout << "Dump #" << n / dump_period << std::endl;
        }
        system.do_step(dt);
    }

    return 0;
}
