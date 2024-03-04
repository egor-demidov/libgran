//
// Created by egor on 3/4/24.
//

#ifndef LIBGRAN_OPERATOR_OVERLOADS_H
#define LIBGRAN_OPERATOR_OVERLOADS_H

#include <utility>

template<typename T1, typename T2>
std::pair<T1, T2> operator + (std::pair<T1, T2> const & lhs, std::pair<T1, T2> const & rhs) {
return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

template<typename T1, typename T2>
std::pair<T1, T2> & operator += (std::pair<T1, T2> & lhs, std::pair<T1, T2> const & rhs) {
    lhs.first += rhs.first;
    lhs.second += rhs.second;
    return lhs;
}

#endif //LIBGRAN_OPERATOR_OVERLOADS_H
