//
// Created by egor on 1/26/24.
//

#ifndef LIBGRAN_MASS_DISTRIBUTION_H
#define LIBGRAN_MASS_DISTRIBUTION_H

#include <vector>

#include <Eigen/Eigen>

Eigen::Vector3d center_of_mass(std::vector<Eigen::Vector3d> const & x);
double radius_of_gyration(std::vector<Eigen::Vector3d> const & x);

#endif //LIBGRAN_MASS_DISTRIBUTION_H
