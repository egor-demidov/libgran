//
// Created by egor on 1/21/24.
//

#ifndef LIBGRAN_GRANULAR_SYSTEM_H
#define LIBGRAN_GRANULAR_SYSTEM_H

#include <libtimestep/rotational_integrator/rotational_integrator.h>
#include <libtimestep/rotational_step_handler/rotational_step_handler.h>
#include <libtimestep/rotational_system/rotational_system.h>

template<typename T1, typename T2>
std::pair<T1, T2> operator + (std::pair<T1, T2> const & lhs, std::pair<T1, T2> const & rhs) {
    return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

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
    typename... force_functors_t>
class granular_system : public rotational_binary_system<field_value_t, real_t, integrator_t, step_handler_t,
        granular_system<field_value_t, real_t, integrator_t, step_handler_t, force_functors_t...>> {
public:
    typedef std::vector<field_value_t> field_container_t;

    granular_system(granular_system const &) = delete;

    granular_system(field_container_t x0, field_container_t v0,
                    field_container_t theta0, field_container_t omega0,
                    real_t t0, field_value_t field_zero, real_t real_zero, step_handler_t<field_container_t, field_value_t> & step_handler,
                    force_functors_t & ... force_functors) :
            rotational_binary_system<field_value_t, real_t, integrator_t, step_handler_t,
                    granular_system<field_value_t, real_t, integrator_t, step_handler_t, force_functors_t...>>
                    (std::move(x0), std::move(v0), std::move(theta0),
                     std::move(omega0), t0, field_zero, real_zero, *this, step_handler), force_functors{force_functors...} {}

    std::pair<field_value_t, field_value_t> compute_accelerations(size_t i, size_t j, field_container_t const & x,
                                                                  field_container_t const & v,
                                                                  field_container_t const & theta,
                                                                  field_container_t const & omega,
                                                                  real_t t) {
        static_assert(std::tuple_size<decltype(force_functors)>::value > 0, "at least one force functor must be provided "
                                                                            "as a template parameter");

        std::pair<field_value_t, field_value_t> accelerations = std::apply([i, j, &x, &v, &theta, &omega, t] (auto & ... e) -> std::pair<field_value_t, field_value_t> {
            return (e(i, j, x, v, theta, omega, t) + ...);
        }, force_functors);
        return accelerations;
    }

private:
    std::tuple<force_functors_t & ...> force_functors;
};

#endif //LIBGRAN_GRANULAR_SYSTEM_H
