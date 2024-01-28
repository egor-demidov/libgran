//
// Created by egor on 1/28/24.
//

#ifndef LIBGRAN_NO_UNARY_FORCE_H
#define LIBGRAN_NO_UNARY_FORCE_H

template <typename field_value_t, typename real_t>
struct no_unary_force_functor {
    explicit no_unary_force_functor(field_value_t zero_field) : zero_field(std::move(zero_field)) {}

    std::pair<field_value_t, field_value_t> operator () (size_t i [[maybe_unused]],
                                                         std::vector<field_value_t> const & x [[maybe_unused]],
                                                         std::vector<field_value_t> const & v [[maybe_unused]],
                                                         std::vector<field_value_t> const & theta [[maybe_unused]],
                                                         std::vector<field_value_t> const & omega [[maybe_unused]],
                                                         real_t t [[maybe_unused]]) const {

        return std::make_pair(zero_field, zero_field);
    }

private:
    const field_value_t zero_field;
};

#endif //LIBGRAN_NO_UNARY_FORCE_H
