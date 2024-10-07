#pragma once

#include "types.hpp"

namespace odes {

class iroot_finder {
public:
    virtual vector_t solve(multi_function_t func, vector_t initial) = 0;
    virtual real_t get_residual() = 0;
    virtual bool precision_was_reached() = 0;
    virtual ~iroot_finder() = default;
};

}
