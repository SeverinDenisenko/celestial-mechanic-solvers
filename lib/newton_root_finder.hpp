#pragma once

#include "jacoby_matrix_evaluator_interface.hpp"
#include "matrix_solver_interface.hpp"
#include "root_finder_interface.hpp"
#include "types.hpp"

namespace odes {

struct newton_root_finder_params_t {
    integer_t max_interations;
    real_t precision;
    uptr<ijacoby_matrix_evaluator> jacoby_matrix_evaluator;
    uptr<imatrix_solver> matrix_solver;
};

class newton_root_finder : public iroot_finder {
public:
    newton_root_finder(newton_root_finder_params_t params);
    vector_t solve(multi_function_t func, vector_t initial) override;
    real_t get_residual() override;
    bool precision_was_reached() override;

private:
    newton_root_finder_params_t params_;
    real_t residual_;
};

}
