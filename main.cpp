#include <iostream>
#include <vector>

#include <Eigen/Eigen>

#include <libgran/granular_system/granular_system.h>
#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>

#include "writer.h"

int main() {
    const double dt = 1e-13;
    const double t_tot = 1.0e-6;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 100;
    const size_t dump_period = n_steps / n_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * cube(r_part) * rho;
    const double inertia = 2.0 / 5.0 * mass * square(r_part);

    // Parameters for the contact model
    const double k = 1000.0;
    const double gamma_n = 5.0e-10;
//    const double gamma_n = 0.0;
    const double mu = 0.1;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2e2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the Van der Waals model
    const double A = 1.0e-20;
    const double h0 = 5.0e-9;

    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0.emplace_back(-1.1 * r_part, 0.0, 0.0);
    x0.emplace_back(1.1 * r_part, 0.0, 0.0);

//    v0.emplace_back(0.0001, 0.0, 0.0);
//    v0.emplace_back(-0.0001, 0.0, 0.0);

    v0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());

    theta0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());

    omega0.emplace_back(0.0, 10000000.0, 0.0);
    omega0.emplace_back(0.0, 0.0, 0.0);

    contact_force_functor<Eigen::Vector3d, double> contact_force(x0.size(), k, gamma_n, k, gamma_t, mu, phi, k, gamma_r,
                                                                 mu, phi, k, gamma_o, mu_o, phi, r_part, mass, inertia,
                                                                 dt, Eigen::Vector3d::Zero(), 0.0);

    hamaker_functor<Eigen::Vector3d, double> hamaker(A, h0, r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half, rotational_step_handler,
        contact_force_functor<Eigen::Vector3d, double>,
        hamaker_functor<Eigen::Vector3d, double>> gran_system(x0, v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                                                                    contact_force, hamaker);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << gran_system.get_v()[0].norm() << std::endl;
            write_particles("run", gran_system.get_x(), gran_system.get_theta(), r_part);
        }

        gran_system.do_step(dt);
    }

    return 0;
}
