#include "full_pivot_gauss_solver.hpp"
#include "jacoby_matrix_evaluator_interface.hpp"
#include "newton_root_finder.hpp"
#include "types.hpp"

#include <iomanip>
#include <memory>

class jacoby_mattrix_evaluator : public odes::ijacoby_matrix_evaluator {
public:
    odes::matrix_t evaluate(odes::multi_function_t, odes::vector_t u) override
    {
        using namespace odes;
        real_t x = u[0];
        real_t y = u[1];

        matrix_t m(2, 2);

        m(0, 0) = 3 * pow(x, 2) - 3 * pow(y, 2);
        m(0, 1) = -6 * x * y;
        m(1, 0) = 6 * x * y;
        m(1, 1) = -3 * pow(y, 2) + 3 * pow(x, 2);

        return m;
    }
};

int main()
{
    odes::set_precision(50);
    odes::newton_root_finder solver(
        odes::newton_root_finder_params_t { .max_interations         = 1000,
                                            .precision               = 1e-30,
                                            .jacoby_matrix_evaluator = std::make_unique<jacoby_mattrix_evaluator>(),
                                            .matrix_solver = std::make_unique<odes::full_pivot_gauss_solver>() });

    using namespace odes;

    multi_function_t equation = [](vector_t u) -> vector_t {
        real_t x = u[0];
        real_t y = u[1];

        real_t z = pow(x, 3) - 3 * x * pow(y, 2) - 1;
        real_t m = -pow(y, 3) + 3 * y * pow(x, 2);

        vector_t w = make_vector<real_t>({ z, m });

        return w;
    };

    std::cerr << equation(make_vector<real_t>({ 1, 0 })) << std::endl;
    std::cerr << equation(make_vector<real_t>({ -0.5, sqrt(real_t(3)) / 2.0 })) << std::endl;
    std::cerr << equation(make_vector<real_t>({ -0.5, -sqrt(real_t(3)) / 2.0 })) << std::endl;

    real_t x1 = -5, x2 = 5;
    real_t y1 = -5, y2 = 5;
    integer_t n = 4000;
    real_t dx   = (x2 - x1) / n;
    real_t dy   = (y2 - y1) / n;

    array_t<array_t<vector_t>> results(n);
    for (auto& row : results) {
        row = array_t<vector_t>(n);
    }

    for (integer_t i = 0; i < n; ++i) {
        for (integer_t j = 0; j < n; ++j) {
            real_t x = x1 + dx * i;
            real_t y = y1 + dy * j;

            vector_t u    = make_vector<real_t>({ x, y });
            results[i][j] = solver.solve(equation, u);
            if (!solver.precision_was_reached()) {
                std::cerr << "(" << x << ", " << y << ") "
                          << "Precision loss" << std::endl;
            }
        }
    }

    array_t<array_t<real_t>> angles(n);
    for (auto& row : angles) {
        row = array_t<real_t>(n);
    }

    for (integer_t i = 0; i < n; ++i) {
        for (integer_t j = 0; j < n; ++j) {
            angles[i][j] = atan2(results[i][j][0], results[i][j][1]);

            std::cout << std::fixed << std::setw(7) << std::setprecision(5) << angles[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
