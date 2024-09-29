#pragma once

#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "mpreal/mpreal.h"

namespace odes {

namespace ublas = boost::numeric::ublas;

using real_t    = mpfr::mpreal;
using integer_t = size_t;
using vector_t  = ublas::vector<real_t>;
template <typename T>
using array_t    = std::vector<T>;
using ode_t      = std::function<vector_t(real_t, vector_t)>;
using string_t   = std::string;
using function_t = std::function<real_t(real_t)>;

template <typename T>
using uptr = std::unique_ptr<T>;

};
