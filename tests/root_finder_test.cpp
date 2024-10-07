#include <memory>

#include "full_pivot_gauss_solver.hpp"
#include "jacoby_matrix_evaluator_interface.hpp"
#include "jacoby_matrix_evaluators.hpp"
#include "newton_root_finder.hpp"
#include "types.hpp"

using namespace odes;

void test(multi_function_t function, vector_t initial, vector_t expected)
{
    odes::jacoby_mattrix_evaluator_params_t jacoby_mattrix_evaluator_params { .step = 1e-5 };
    uptr<odes::ijacoby_matrix_evaluator> jacoby_matrix_evaluator
        = std::make_unique<odes::fourth_order_jacoby_mattrix_evaluator>(jacoby_mattrix_evaluator_params);
    odes::newton_root_finder_params_t params { .max_interations         = 100,
                                               .precision               = 1e-10,
                                               .jacoby_matrix_evaluator = std::move(jacoby_matrix_evaluator),
                                               .matrix_solver = std::make_unique<odes::full_pivot_gauss_solver>() };
    odes::newton_root_finder root_finder(std::move(params));

    vector_t solution = root_finder.solve(function, initial);

    if (norm(solution - expected) < 1e-10) {
        std::cout << "Good." << solution << std::endl;
    } else {
        std::cout << "Bad." << solution << std::endl;
    }
}

int main()
{
    test(
        [](vector_t a) -> vector_t { return a; }, make_vector<real_t>({ 1.0, 1.0 }), make_vector<real_t>({ 0.0, 0.0 }));
    test(
        [](vector_t a) -> vector_t {
            vector_t res(1);
            res[0] = a[0] * a[0] - 4.0;
            return res;
        },
        make_vector<real_t>({ 3.0 }),
        make_vector<real_t>({ 2.0 }));
    test(
        [](vector_t a) -> vector_t {
            vector_t res(2);
            res[0] = a[0] * a[0] + a[1] * a[1] - 4.0;
            res[1] = a[0] * a[1] - 1.0;
            return res;
        },
        make_vector<real_t>({ 2.0, 0.5 }),
        make_vector<real_t>({ 1.93185165257813657350, 0.51763809020504152470 }));

    return 0;
}
