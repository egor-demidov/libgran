//
// Created by egor on 1/26/24.
//

#include "compute_energy.h"

double compute_kinetic_energy(std::vector<Eigen::Vector3d> const & v,
                              std::vector<Eigen::Vector3d> const & omega,
                              double mass, double inertia) {

    double energy = 0.0;

    for (size_t i = 0; i < v.size(); i ++) {
        energy += mass * v[i].dot(v[i]) / 2.0;
        energy += inertia * omega[i].dot(omega[i]) / 2.0;
    }

    return energy;
}

double compute_linear_momentum(std::vector<Eigen::Vector3d> const & v, double mass) {
    Eigen::Vector3d momentum = Eigen::Vector3d::Zero();

    std::for_each(v.begin(), v.end(), [&momentum] (auto const & v) {
        momentum += v;
    });

    return mass * momentum.norm();
}
