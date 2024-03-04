//
// Created by egor on 3/4/24.
//

#ifndef LIBGRAN_UNARY_FORCE_CONTAINER_H
#define LIBGRAN_UNARY_FORCE_CONTAINER_H

#include <vector>

template <
        typename field_value_t,
        typename real_t,
        typename... unary_force_functors_t>
struct unary_force_functor_container {
    typedef std::vector<field_value_t> field_container_t;

    unary_force_functor_container(unary_force_functors_t & ... force_functors)
            : unary_force_functors{force_functors...} {}

    std::pair<field_value_t, field_value_t> operator () (size_t i, field_container_t const & x,
                                                         field_container_t const & v,
                                                         field_container_t const & theta,
                                                         field_container_t const & omega,
                                                         real_t t) {
        static_assert(std::tuple_size<decltype(unary_force_functors)>::value > 0, "at least one unary force functor must be provided "
                                                                                  "as a template parameter");

        std::pair<field_value_t, field_value_t> accelerations = std::apply([i, &x, &v, &theta, &omega, t] (auto & ... e) -> std::pair<field_value_t, field_value_t> {
            return (e(i, x, v, theta, omega, t) + ...);
        }, unary_force_functors);
        return accelerations;
    }

    std::tuple<unary_force_functors_t & ...> unary_force_functors;
};

#endif //LIBGRAN_UNARY_FORCE_CONTAINER_H
