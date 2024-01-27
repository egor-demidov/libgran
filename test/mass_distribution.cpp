//
// Created by egor on 1/26/24.
//

#include "mass_distribution.h"

Eigen::Vector3d center_of_mass(std::vector<Eigen::Vector3d> const & x) {
    Eigen::Vector3d center = Eigen::Vector3d::Zero();

    std::for_each(x.begin(), x.end(), [&center] (auto const & x) {
        center += x;
    });

    return center / double(x.size());
}

double radius_of_gyration(std::vector<Eigen::Vector3d> const & x) {
    Eigen::Vector3d center = center_of_mass(x);

    double rg = 0.0;
    std::for_each(x.begin(), x.end(), [&rg, center] (auto const & x) {
        auto distance = x - center;
        rg += distance.dot(distance);
    });

    return sqrt(rg / double(x.size()));
}
