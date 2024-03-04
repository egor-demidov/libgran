//
// Created by egor on 3/4/24.
//

#ifndef LIBGRAN_BINARY_FORCE_CONTAINER_H
#define LIBGRAN_BINARY_FORCE_CONTAINER_H

#include <vector>

template <
        typename field_value_t,
        typename real_t,
        typename... binary_force_functors_t>
struct binary_force_functor_container {
    typedef std::vector<field_value_t> field_container_t;

    binary_force_functor_container(binary_force_functors_t & ... force_functors)
            : binary_force_functors{force_functors...} {}

    std::pair<field_value_t, field_value_t> operator () (size_t i, size_t j, field_container_t const & x,
                                                         field_container_t const & v,
                                                         field_container_t const & theta,
                                                         field_container_t const & omega,
                                                         real_t t) {
        static_assert(std::tuple_size<decltype(binary_force_functors)>::value > 0, "at least one binary force functor must be provided "
                                                                                   "as a template parameter");

        std::pair<field_value_t, field_value_t> accelerations = std::apply([i, j, &x, &v, &theta, &omega, t] (auto & ... e) -> std::pair<field_value_t, field_value_t> {
            return (e(i, j, x, v, theta, omega, t) + ...);
        }, binary_force_functors);
        return accelerations;
    }

    std::tuple<binary_force_functors_t & ...> binary_force_functors;
};


#endif //LIBGRAN_BINARY_FORCE_CONTAINER_H
