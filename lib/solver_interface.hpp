#pragma once

#include "types.hpp"

namespace odes {

struct ode_params_t {
    real_t t0;
    real_t t1;
    real_t dt;
    vector_t x0;
    ode_t ode;
};

class isolver {
public:
    virtual void step() noexcept                        = 0;
    virtual const vector_t& current() const noexcept    = 0;
    virtual const real_t& current_time() const noexcept = 0;
    virtual const real_t& begin_time() const noexcept   = 0;
    virtual const real_t& end_time() const noexcept     = 0;
    virtual ~isolver()                                  = default;
};

} // namespace odes
