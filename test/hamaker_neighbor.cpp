//
// Created by egor on 3/4/24.
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
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system_neighbor_list.h>
