//
// Created by egor on 1/24/24.
//

#ifndef LIBGRAN_SIMPLIFIED_HAMAKER_FORCE_H
#define LIBGRAN_SIMPLIFIED_HAMAKER_FORCE_H

template <typename field_value_t, typename real_t>
struct simplified_hamaker_functor {
    simplified_hamaker_functor(real_t A,
                               real_t h0,
                               real_t r_part,
                               real_t mass,
                               field_value_t field_zero,
                               real_t real_zero) :
       A(A), h0(h0), r_part(r_part), mass(mass), real_zero(real_zero), field_zero(field_zero) {}

    std::pair<field_value_t, field_value_t> operator () (size_t i, size_t j,
                                                         std::vector<field_value_t> const & x, std::vector<field_value_t> const & v [[maybe_unused]],
                                                         std::vector<field_value_t> const & theta [[maybe_unused]], std::vector<field_value_t> const & omega [[maybe_unused]], real_t t [[maybe_unused]]) const {

        field_value_t r = (x[j] - x[i]); // Distance vector
        real_t h = r.norm() - 2.0 * r_part; // Surface separation
        if (h < h0)
            h = h0;

        field_value_t f = 2.0 * A * r_part / square(h) * r.normalized();

        return std::make_pair(f/mass, field_zero);
    }

private:
    const real_t A, h0, r_part, mass, real_zero;
    const field_value_t field_zero;
};

#endif //LIBGRAN_SIMPLIFIED_HAMAKER_FORCE_H
