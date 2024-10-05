#pragma once

#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "mpreal/mpreal.h"

namespace odes {

namespace ublas = boost::numeric::ublas;

#ifdef SNUMLIB_USE_MPFR
using real_t = mpfr::mpreal;
#else
using real_t = double;
#endif

using integer_t = size_t;
using vector_t  = ublas::vector<real_t>;
using matrix_t  = ublas::matrix<real_t>;
template <typename T>
using array_t          = std::vector<T>;
using ode_t            = std::function<vector_t(real_t, vector_t)>;
using string_t         = std::string;
using function_t       = std::function<real_t(real_t)>;
using multi_function_t = std::function<vector_t(vector_t)>;

inline vector_t operator*(matrix_t m, vector_t v)
{
    return boost::numeric::ublas::prod(m, v);
}

inline matrix_t operator*(matrix_t m1, matrix_t m2)
{
    return boost::numeric::ublas::prod(m1, m2);
}

template <typename T>
using uptr = std::unique_ptr<T>;

inline void set_precision(integer_t digits)
{
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
}

};
