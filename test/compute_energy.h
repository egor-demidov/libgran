//
// Created by egor on 1/26/24.
//

#ifndef LIBGRAN_COMPUTE_ENERGY_H
#define LIBGRAN_COMPUTE_ENERGY_H

#include <vector>

#include <Eigen/Eigen>

double compute_total_kinetic_energy(std::vector<Eigen::Vector3d> const & v,
                                    std::vector<Eigen::Vector3d> const & omega,
                                    double mass, double inertia);

double compute_linear_momentum(std::vector<Eigen::Vector3d> const & vs, double mass);

double compute_rotational_kinetic_energy(std::vector<Eigen::Vector3d> const & omegas, double inertia);

double compute_translational_kinetic_energy(std::vector<Eigen::Vector3d> const & vs, double m);

double compute_angular_momentum(std::vector<Eigen::Vector3d> const & xs, std::vector<Eigen::Vector3d> const & vs, double m,
                                std::vector<Eigen::Vector3d> const & omegas, double inertia);


#endif //LIBGRAN_COMPUTE_ENERGY_H
