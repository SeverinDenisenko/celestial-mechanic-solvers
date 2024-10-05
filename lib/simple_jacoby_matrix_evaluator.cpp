#include "simple_jacoby_matrix_evaluator.hpp"
#include "types.hpp"

namespace odes {
simple_jacoby_mattrix_evaluator::simple_jacoby_mattrix_evaluator(simple_jacoby_mattrix_evaluator_params_t params)
    : params_(params)
{
}

matrix_t simple_jacoby_mattrix_evaluator::evaluate(multi_function_t func, vector_t point)
{
    vector_t x = func(point);
    matrix_t result(point.size(), x.size());

    vector_t a = point;
    vector_t b = point;

    for (integer_t i = 0; i < result.size1(); ++i) {
        a[i] += params_.step;
        b[i] -= params_.step;

        vector_t derived = (func(a) - func(b)) / (real_t(2) * params_.step);

        for (integer_t j = 0; j < result.size2(); ++j) {
            result(i, j) = derived[j];
        }

        a[i] -= params_.step;
        b[i] += params_.step;
    }

    return result;
}

}
