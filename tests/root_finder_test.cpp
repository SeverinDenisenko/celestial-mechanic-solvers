#include <cstdlib>
#include <memory>

#include "full_pivot_gauss_solver.hpp"
#include "jacoby_matrix_evaluator_interface.hpp"
#include "newton_root_finder.hpp"
#include "simple_jacoby_matrix_evaluator.hpp"
#include "types.hpp"

using namespace odes;

int main()
{
    vector_t a(2);
    vector_t b(2);
    vector_t c;
    a[0]    = 1;
    a[1]    = 1;
    b[0]    = 0;
    b[1]    = 0;
    multi_function_t func = [](vector_t a) -> vector_t{
        return a;
    };

    odes::simple_jacoby_mattrix_evaluator_params_t jacoby_mattrix_evaluator_params { .step = 1e-5 };
    uptr<odes::ijacoby_matrix_evaluator> jacoby_matrix_evaluator
        = std::make_unique<odes::simple_jacoby_mattrix_evaluator>(jacoby_mattrix_evaluator_params);
    odes::newton_root_finder_params_t params { .precision       = 1e-5,
                                               .max_interations = 10,
                                               .matrix_solver   = std::make_unique<odes::full_pivot_gauss_solver>(),
                                               .jacoby_matrix_evaluator = std::move(jacoby_matrix_evaluator) };
    odes::newton_root_finder root_finder(std::move(params));

    c = root_finder.solve(func, a);

    if (norm(c - b) < 1e-10) {
        std::cout << "Good." << c << std::endl;
        exit(0);
    } else {
        std::cout << "Bad." << c << std::endl;
        exit(1);
    }
}
