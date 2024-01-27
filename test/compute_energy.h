//
// Created by egor on 1/26/24.
//

#ifndef LIBGRAN_COMPUTE_ENERGY_H
#define LIBGRAN_COMPUTE_ENERGY_H

#include <vector>

#include <Eigen/Eigen>

double compute_kinetic_energy(std::vector<Eigen::Vector3d> const & v,
                              std::vector<Eigen::Vector3d> const & omega,
                              double mass, double inertia);

double compute_linear_momentum(std::vector<Eigen::Vector3d> const & v, double mass);

#endif //LIBGRAN_COMPUTE_ENERGY_H
