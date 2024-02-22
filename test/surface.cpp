//
// Created by egor on 2/22/24.
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
#include <libgran/contact_force/surface_contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/hamaker_force/surface_hamaker_force.h>
#include <libgran/granular_system/granular_system.h>
#include <libgran/surface_force/triangular_facet.h>

#include "../writer.h"
#include "compute_energy.h"

using facet_t = std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>;

using contact_force_functor_t = contact_force_functor<Eigen::Vector3d, double>;
using binary_force_container_t =
    binary_force_functor_container<Eigen::Vector3d, double,
    contact_force_functor_t>;

using surface_contact_force_functor_t = surface_contact_force_functor<Eigen::Vector3d, double>;
using surface_hamaker_functor_t = surface_hamaker_functor<Eigen::Vector3d, double>;
using triangular_facet_t = triangular_facet<Eigen::Vector3d, double, surface_contact_force_functor_t, surface_hamaker_functor_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double, triangular_facet_t>;

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 5.0e-7;
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
    const double mu = 0.5;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the Van der Waals model
    const double A = 0.0;
    const double h0 = 1.0e-9;

    const double a = 50.0 * r_part;
    const facet_t facet {
        {-a, -a, 0},
        {-a, a, 0},
        {a, -a, 0}
    };

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
    x0.emplace_back(0.5*r_part, 0.5*r_part, 2.5 * r_part);
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), -0.5 * Eigen::Vector3d::UnitZ());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    contact_force_functor_t contact_force_model(x0.size(),
        k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
        mu_o, phi, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);
    binary_force_container_t
        binary_force_functors {contact_force_model};

    surface_contact_force_functor_t surface_contact_force(x0.size(),
        k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
        mu_o, phi, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);
    surface_hamaker_functor_t surface_hamaker_force(A, h0, r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    triangular_facet_t facet_model(Eigen::Vector3d::Zero(), facet, surface_contact_force, surface_hamaker_force);
    unary_force_container_t unary_force_functors {facet_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
        step_handler_instance;

    granular_system_t system(x0,
        v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
        step_handler_instance, binary_force_functors, unary_force_functors);

    for (size_t n = 0; n < n_steps; n ++) {
        system.do_step(dt);
    }

    const double ke_target = 1.0922e-21;
    const double ke = compute_total_kinetic_energy(system.get_v(), system.get_omega(), mass, inertia);

    const double ke_tolerance = 1.0; // Percent

    if (abs(ke - ke_target) / ke_target > ke_tolerance / 100.0)
        return EXIT_FAILURE;

    return 0;
}
