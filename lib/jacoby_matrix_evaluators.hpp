#pragma once

#include "jacoby_matrix_evaluator_interface.hpp"
#include "types.hpp"

namespace odes {
struct jacoby_mattrix_evaluator_params_t {
    real_t step;
};

class second_order_jacoby_mattrix_evaluator : public ijacoby_matrix_evaluator {
public:
    second_order_jacoby_mattrix_evaluator(jacoby_mattrix_evaluator_params_t params);
    matrix_t evaluate(multi_function_t func, vector_t point) override;

private:
    jacoby_mattrix_evaluator_params_t params_;
};

class fourth_order_jacoby_mattrix_evaluator : public ijacoby_matrix_evaluator {
public:
    fourth_order_jacoby_mattrix_evaluator(jacoby_mattrix_evaluator_params_t params);
    matrix_t evaluate(multi_function_t func, vector_t point) override;

private:
    jacoby_mattrix_evaluator_params_t params_;
};

}
