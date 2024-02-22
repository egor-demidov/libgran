//
// Created by egor on 2/22/24.
//

#ifndef SURFACE_CONTACT_FORCE_H
#define SURFACE_CONTACT_FORCE_H

template <typename field_value_t, typename real_t>
struct surface_contact_force_functor {
    surface_contact_force_functor(real_t k, real_t gamma_n, real_t r_part, real_t mass, field_value_t zero_field) :
        k(k), gamma_n(gamma_n), r_part(r_part), mass(mass), zero_field(zero_field) {}

    std::pair<field_value_t, field_value_t> operator () (size_t i,
                                                         field_value_t const & x_facet,
                                                         std::vector<field_value_t> const & x,
                                                         std::vector<field_value_t> const & v,
                                                         std::vector<field_value_t> const & theta [[maybe_unused]],
                                                         std::vector<field_value_t> const & omega [[maybe_unused]],
                                                         real_t t [[maybe_unused]]) const {

        real_t overlap = r_part - (x[i] - x_facet).norm();

        if (overlap < 0)
            return std::make_pair(zero_field, zero_field);

        field_value_t n = (x[i] - x_facet).normalized();

        field_value_t f = (k * overlap + v[i].dot(n) * gamma_n) * n;
        return std::make_pair(f / mass, zero_field);
    }

private:
    const real_t k, gamma_n, r_part, mass;
    const field_value_t zero_field;
};

#endif //SURFACE_CONTACT_FORCE_H
