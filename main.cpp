#include <iostream>
#include <vector>

#include <Eigen/Eigen>

#include <libgran/granular_system/granular_system.h>
#include <libgran/contact_force/contact_force.h>

#include "writer.h"

//struct F1 {
//    std::pair<Eigen::Vector3d, Eigen::Vector3d> operator() (size_t i [[maybe_unused]],
//            size_t j [[maybe_unused]],
//            std::vector<Eigen::Vector3d> const & x [[maybe_unused]],
//            std::vector<Eigen::Vector3d> const & v [[maybe_unused]],
//            std::vector<Eigen::Vector3d> const & theta [[maybe_unused]],
//            std::vector<Eigen::Vector3d> const & omega [[maybe_unused]],
//            double t [[maybe_unused]]) const {
//        return std::make_pair(Eigen::Vector3d{0.0, 10.0, 0.0}, Eigen::Vector3d::Zero());
//    }
//};
//
//struct F2 {
//    std::pair<Eigen::Vector3d, Eigen::Vector3d> operator() (size_t i [[maybe_unused]],
//            size_t j [[maybe_unused]],
//            std::vector<Eigen::Vector3d> const & x [[maybe_unused]],
//            std::vector<Eigen::Vector3d> const & v [[maybe_unused]],
//            std::vector<Eigen::Vector3d> const & theta [[maybe_unused]],
//            std::vector<Eigen::Vector3d> const & omega [[maybe_unused]],
//            double t [[maybe_unused]]) const {
//        return std::make_pair(Eigen::Vector3d{20.0, 0.0, 0.0}, Eigen::Vector3d::Zero());
//    }
//};

int main() {
    const double dt = 0.0001;
    const double t_tot = 5.0;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 100;
    const size_t dump_period = n_steps / n_dumps;

    const double k = 1000.0;
    const double gamma_n = 0.1;
    const double r_part = 0.1;
    const double mass = 1.0;
    const double inertia = 0.01;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0.emplace_back(-2.0 * r_part, 0.0, 0.0);
    x0.emplace_back(2.0 * r_part, 0.0, 0.0);

    v0.emplace_back(0.1, 0.0, 0.0);
    v0.emplace_back(-0.1, 0.0, 0.0);

    theta0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());

    omega0.emplace_back(0.0, 1.0, 0.0);
    omega0.emplace_back(0.0, 0.0, 0.0);

    contact_force_functor<Eigen::Vector3d, double> contact_force(x0.size(),
                                                                 k,
                                                                 gamma_n,
                                                                 k,
                                                                 gamma_t,
                                                                 mu,
                                                                 phi,
                                                                 k,
                                                                 gamma_r,
                                                                 mu,
                                                                 phi,
                                                                 k,
                                                                 gamma_o,
                                                                 mu_o,
                                                                 phi,
                                                                 r_part,
                                                                 mass,
                                                                 inertia,
                                                                 dt,
                                                                 Eigen::Vector3d::Zero(),
                                                                 0.0);

    granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half, rotational_step_handler,
        contact_force_functor<Eigen::Vector3d, double>> gran_system(x0, v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                                                                    contact_force);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << gran_system.get_omega()[0].norm() << std::endl;
            write_particles("run", gran_system.get_x(), gran_system.get_theta(), r_part);
        }

        gran_system.do_step(dt);
    }

//    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
//
//    x0.emplace_back(Eigen::Vector3d::Zero());
//    v0.emplace_back(Eigen::Vector3d::Zero());
//    theta0.emplace_back(Eigen::Vector3d::Zero());
//    omega0.emplace_back(Eigen::Vector3d::Zero());

//    std::vector<double> x0, v0, theta0, omega0;
//
//    x0.emplace_back(0.0);
//    v0.emplace_back(0.0);
//    theta0.emplace_back(0.0);
//    omega0.emplace_back(0.0);

//    granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half, rotational_step_handler, F1, F2>
//            system(x0, v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0);
//
//    auto [a, alpha] = system.compute_accelerations(0, 0, system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), 0.0);
//    std::cout << a << std::endl;

    return 0;
}
