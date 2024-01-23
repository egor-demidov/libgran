//
// Created by egor on 1/23/24.
//

#ifndef LIBGRAN_WRITER_H
#define LIBGRAN_WRITER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <Eigen/Eigen>

void write_particles(const std::string & dir, std::vector<Eigen::Vector3d> const & x, std::vector<Eigen::Vector3d> const & theta, double r_part);

#endif //LIBGRAN_WRITER_H
