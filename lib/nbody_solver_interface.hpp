#pragma once

#include "types.hpp"

namespace odes {

constexpr integer_t dimention = 3;

struct body_t {
    body_t()
        : r(dimention)
        , v(dimention)
        , a(dimention)
        , j(dimention)
    {
    }

    vector_t r;
    vector_t v;
    vector_t a;
    vector_t j;
};

using nbody_t                = array_t<body_t>;
using interaction_function_t = std::function<vector_t(vector_t r, vector_t v, real_t m)>;

struct nbody_params_t {
    real_t t0;
    real_t t1;
    real_t dt;
    integer_t count;
    array_t<real_t> m;
    nbody_t state;
    interaction_function_t a;
    interaction_function_t j;
};

class inbody_solver {
public:
    virtual void step() noexcept                        = 0;
    virtual const nbody_t& current() const noexcept     = 0;
    virtual const real_t& current_time() const noexcept = 0;
    virtual const real_t& begin_time() const noexcept   = 0;
    virtual const real_t& end_time() const noexcept     = 0;
    virtual ~inbody_solver()                            = default;
};

}
