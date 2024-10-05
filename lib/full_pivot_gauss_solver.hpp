#pragma once

#include "matrix_solver_interface.hpp"

namespace odes {

class full_pivot_gauss_solver : public imatrix_solver {
public:
    vector_t solve(matrix_t m, vector_t b) override;
};

}
