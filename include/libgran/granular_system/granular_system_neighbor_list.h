//
// Created by egor on 3/4/24.
//

#ifndef LIBGRAN_GRANULAR_SYSTEM_NEIGHBOR_LIST_H
#define LIBGRAN_GRANULAR_SYSTEM_NEIGHBOR_LIST_H


#include <libtimestep/rotational_integrator/rotational_integrator.h>
#include <libtimestep/rotational_step_handler/rotational_step_handler.h>

#ifndef LIBGRAN_USE_OMP
#define binary_system_implementation rotational_binary_system
#include <libtimestep/rotational_system/rotational_binary_system.h>
#error "Neighbor lists are not supported with C++ 17 algorithms yet"
#else
#define binary_system_implementation rotational_binary_system_neighbors_omp
#include <libtimestep/rotational_system/rotational_binary_system_neighbors_omp.h>
#pragma message("libgran: using OpenMP parallelization")
#endif //LIBGRAN_USE_OMP

#include "../utility/operator_overloads.h"
#include "binary_force_container.h"
#include "unary_force_container.h"

template <
        typename field_value_t,
        typename real_t,
        template <
        typename _field_container_t,
        typename _field_value_t,
        typename _real_t,
        typename _functor_t,
        template <
        typename __field_container_t,
        typename __field_value_t>
        typename _step_handler_t>
        typename integrator_t,
        template <
        typename _field_container_t,
        typename _field_value_t>
        typename step_handler_t,
        typename binary_force_functor_container_t,
        typename unary_force_functor_container_t>
class granular_system_neighbor_list : public binary_system_implementation<field_value_t, real_t, integrator_t, step_handler_t,
        granular_system_neighbor_list<field_value_t, real_t, integrator_t, step_handler_t,
                binary_force_functor_container_t, unary_force_functor_container_t>,
        (std::tuple_size<decltype(unary_force_functor_container_t::unary_force_functors)>::value > 0)> {
public:
    typedef std::vector<field_value_t> field_container_t;

    granular_system_neighbor_list(granular_system_neighbor_list const &) = delete;

    granular_system_neighbor_list(long n_part, real_t r_verlet, field_container_t x0, field_container_t v0,
                    field_container_t theta0, field_container_t omega0,
                    real_t t0, field_value_t field_zero, real_t real_zero, step_handler_t<field_container_t, field_value_t> & step_handler,
                    binary_force_functor_container_t binary_force_functors, unary_force_functor_container_t unary_force_functors) :
            binary_system_implementation<field_value_t, real_t, integrator_t, step_handler_t,
                    granular_system_neighbor_list<field_value_t, real_t, integrator_t, step_handler_t, binary_force_functor_container_t, unary_force_functor_container_t>,
                    (std::tuple_size<decltype(unary_force_functor_container_t::unary_force_functors)>::value > 0)>
                    (n_part, r_verlet, std::move(x0), std::move(v0), std::move(theta0),
                     std::move(omega0), t0, field_zero, real_zero, *this, step_handler),
            binary_force_functors(std::move(binary_force_functors)),
            unary_force_functors(std::move(unary_force_functors)) {}

    std::pair<field_value_t, field_value_t> compute_accelerations(size_t i, size_t j, field_container_t const & x,
                                                                  field_container_t const & v,
                                                                  field_container_t const & theta,
                                                                  field_container_t const & omega,
                                                                  real_t t) {

        return binary_force_functors(i, j, x, v, theta, omega, t);
    }

    std::pair<field_value_t, field_value_t> compute_accelerations(size_t i, field_container_t const & x,
                                                                  field_container_t const & v,
                                                                  field_container_t const & theta,
                                                                  field_container_t const & omega,
                                                                  real_t t) {

        return unary_force_functors(i, x, v, theta, omega, t);
    }

private:
    binary_force_functor_container_t binary_force_functors;
    unary_force_functor_container_t unary_force_functors;
};


#endif //LIBGRAN_GRANULAR_SYSTEM_NEIGHBOR_LIST_H
