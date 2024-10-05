#pragma once

#include "types.hpp"

namespace odes {
class ijacoby_matrix_evaluator {
public:
    virtual matrix_t evaluate(multi_function_t func, vector_t point) = 0;
    virtual ~ijacoby_matrix_evaluator() = default;
};
}
