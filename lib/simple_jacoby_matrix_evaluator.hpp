#pragma once

#include "jacoby_matrix_evaluator_interface.hpp"
#include "types.hpp"

namespace odes {
struct simple_jacoby_mattrix_evaluator_params_t {
    real_t step;
};

class simple_jacoby_mattrix_evaluator : public ijacoby_matrix_evaluator {
public:
    simple_jacoby_mattrix_evaluator(simple_jacoby_mattrix_evaluator_params_t params);
    matrix_t evaluate(multi_function_t func, vector_t point);

private:
    simple_jacoby_mattrix_evaluator_params_t params_;
};
}
