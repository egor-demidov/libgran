//
// Created by egor on 2/22/24.
//

#ifndef SURFACE_HAMAKER_FORCE_H
#define SURFACE_HAMAKER_FORCE_H

template <typename field_value_t, typename real_t>
struct surface_hamaker_functor {
    surface_hamaker_functor(real_t A,           // Hamaker constant
                    real_t h0,                  // Saturation separation
                    real_t r_part,              // Particle radius
                    real_t mass,                // Particle mass
                    field_value_t field_zero,   // Zero value of field_value_t
                    real_t real_zero) :         // Zero value of real_t
        A(A), h0(h0), r_part(r_part), mass(mass), real_zero(real_zero), field_zero(field_zero) {}

    std::pair<field_value_t, field_value_t> operator () (size_t i, field_value_t const & x_facet,
            std::vector<field_value_t> const & x, std::vector<field_value_t> const & v [[maybe_unused]],
            std::vector<field_value_t> const & theta [[maybe_unused]], std::vector<field_value_t> const & omega [[maybe_unused]], real_t t [[maybe_unused]]) const {

        field_value_t r = (x_facet - x[i]); // Distance vector
        real_t h = r.norm() - 2.0 * r_part; // Surface separation
        if (h < h0)
            h = h0;
        field_value_t f = -A / 6.0 * ((4.0*r_part+2.0*h)/(4.0*r_part+h)/h
                - 2.0/(2.0*r_part+h)
                - 4.0*square(r_part)/cube(2.0*r_part+h)
                -2.0*square(r_part)*(4.0*r_part+2.0*h)/square(4.0*r_part+h)/square(h))
                * r.normalized();

        return std::make_pair(f/mass, field_zero);
    }

private:
    const real_t A, h0, r_part, mass, real_zero;
    const field_value_t field_zero;
};

#endif //SURFACE_HAMAKER_FORCE_H
