#include "jacoby_matrix_evaluators.hpp"
#include "types.hpp"

namespace odes {
second_order_jacoby_mattrix_evaluator::second_order_jacoby_mattrix_evaluator(jacoby_mattrix_evaluator_params_t params)
    : params_(params)
{
}

matrix_t second_order_jacoby_mattrix_evaluator::evaluate(multi_function_t f, vector_t x)
{
    vector_t y = f(x);
    matrix_t result(x.size(), y.size());

    for (integer_t i = 0; i < result.size1(); ++i) {
        vector_t h = x * 0.0;

        h[i] = params_.step;

        vector_t derived = (f(x + h) - f(x - h)) / (2.0 * params_.step);

        for (integer_t j = 0; j < result.size2(); ++j) {
            result(i, j) = derived[j];
        }
    }

    return result;
}

fourth_order_jacoby_mattrix_evaluator::fourth_order_jacoby_mattrix_evaluator(jacoby_mattrix_evaluator_params_t params)
    : params_(params)
{
}

matrix_t fourth_order_jacoby_mattrix_evaluator::evaluate(multi_function_t f, vector_t x)
{
    vector_t y = f(x);
    matrix_t result(x.size(), y.size());

    for (integer_t i = 0; i < result.size1(); ++i) {
        vector_t h = x * 0.0;

        h[i] = params_.step;

        vector_t derived = (-f(x + h * 2.0) + 8.0 * f(x + h) - 8.0 * f(x - h) + f(x - 2.0 * h)) / (12.0 * params_.step);

        for (integer_t j = 0; j < result.size2(); ++j) {
            result(i, j) = derived[j];
        }
    }

    return result;
}

}
