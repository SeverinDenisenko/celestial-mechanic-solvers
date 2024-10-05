#pragma once

#include "types.hpp"

namespace odes {

class imatrix_solver {
public:
    virtual vector_t solve(matrix_t m, vector_t b) = 0;
    virtual ~imatrix_solver() = default;
};

}
