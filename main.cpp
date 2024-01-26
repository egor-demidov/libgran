#include <iostream>
#include <vector>
#include <random>

#include <Eigen/Eigen>

#include <boost/lambda/construct.hpp>

#include <libgran/granular_system/granular_system.h>
#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/hamaker_force/simplified_hamaker_force.h>
#include <libgran/sinter_bridge/sinter_bridge.h>

#include "writer.h"

// The driver program is responsible for:
// 1) Instantiating the force models
// 2) Instantiating the step handler
// 3) Initializing the initial positions and velocities
// 4) Initializing sinter bridges and sinter-contact map
// 4) Instantiating granular system
// 5) Performing the time steps

template <typename field_container_t, typename field_value_t>
using sinter_step_handler_double = sinter_step_handler<field_container_t, field_value_t, double>;

std::mt19937_64 mt(0);
std::uniform_real_distribution<double> dist(-1.0, 1.0);

Eigen::Vector3d generate_random_unit_vector() {
    Eigen::Vector3d vec;
    do {
        vec = {dist(mt), dist(mt), dist(mt)};
    } while (vec.norm() == 0);
    return vec.normalized();
}

int main() {
    const double dt = 1e-16;
    const double t_tot = 3.0e-8;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 100;
    const size_t dump_period = n_steps / n_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * cube(r_part) * rho;
    const double inertia = 2.0 / 5.0 * mass * square(r_part);

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
//    const double A = 0.0;
    const double h0 = 1.0e-9;

    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0.emplace_back(-1.00 * r_part, 0.0, 0.0);
    x0.emplace_back(1.00 * r_part, 0.0, 0.0);
    x0.emplace_back(1.00 * r_part, 0.0, 2.00 * r_part);

    v0.emplace_back(1.0, 0.0, 0.0);
    v0.emplace_back(0.0, 0.0, 0.0);
    v0.emplace_back(0.0, 0.0, 0.0);

//    v0.resize(x0.size());
//    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());

    theta0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());

//    omega0.emplace_back(0.0, 100000000.0, 0.0);
    omega0.emplace_back(0.0, 0.0, 0.0);
    omega0.emplace_back(0.0, 0.0, 0.0);
    omega0.emplace_back(0.0, 0.0, 0.0);

    contact_force_functor<Eigen::Vector3d, double> contact_force(x0.size(), k, gamma_n, k, gamma_t, mu, phi, k, gamma_r,
                                                                 0.0, phi, k, gamma_o, mu_o, phi, r_part, mass, inertia,
                                                                 dt, Eigen::Vector3d::Zero(), 0.0);

    hamaker_functor<Eigen::Vector3d, double> hamaker(A, h0, r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    sinter_functor<Eigen::Vector3d, double> sinter_bridge(x0.size(), k, gamma_n, r_part, mass, inertia,
                                                          Eigen::Vector3d::Zero(), 0.0, 1.0e-9, x0.begin(),
                                                          x0.end(), contact_force.enabled_contacts, generate_random_unit_vector);

    auto step_handler = sinter_bridge.get_step_handler<std::vector<Eigen::Vector3d>>();

//    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler;

    granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half, sinter_step_handler_double,
        contact_force_functor<Eigen::Vector3d, double>,
        hamaker_functor<Eigen::Vector3d, double>,
        sinter_functor<Eigen::Vector3d, double>> gran_system(x0, v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                                                                   step_handler, contact_force, hamaker, sinter_bridge);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << gran_system.get_v()[0].norm() << std::endl;
            write_particles("run", gran_system.get_x(), gran_system.get_theta(), r_part);
            write_spring_connectors("run", r_part, sinter_bridge.spring_connectors, contact_force.enabled_contacts);
        }

        gran_system.do_step(dt);
    }

    return 0;
}
