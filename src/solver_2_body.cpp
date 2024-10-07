#include "adams_extrapolation_solver.hpp"
#include "adams_interpolation_solver.hpp"
#include "adams_predictor_corrector_solver.hpp"

#include "boole_intergrator.hpp"
#include "full_pivot_gauss_solver.hpp"
#include "jacoby_matrix_evaluators.hpp"
#include "newton_root_finder.hpp"
#include "root_finder_interface.hpp"
#include "solver_interface.hpp"
#include "types.hpp"

#include <memory>

odes::real_t calc_orbital_period(odes::real_t r, odes::real_t v)
{
    odes::real_t a = odes::real_t(1) / (odes::real_t(2) / r - v * v);
    odes::real_t t = pow(a, odes::real_t(3) / odes::real_t(2)) * odes::real_t(2) * odes::real_t(M_PI);
    return t;
};

int main()
{
    size_t digits = 50;
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));

    odes::ode_t ode = [](odes::real_t t, odes::vector_t x) -> odes::vector_t {
        odes::vector_t y(4);

        // pow is overrided for mpfr::mpreal
        odes::real_t z = pow(x[0] * x[0] + x[1] * x[1], 3.0 / 2.0);

        y[2] = -x[0] / z;
        y[3] = -x[1] / z;
        y[0] = x[2];
        y[1] = x[3];

        return y;
    };

    odes::integer_t order = 3;
    odes::boole_integrator_params_t integrator_params { .order = 1'000'000 };

    odes::adams_interpolation_coefficients_params_t interpolation_coefficients_params {
        .order = order, .integrator = std::make_unique<odes::boole_integrator>(integrator_params)
    };
    odes::adams_interpolation_coefficients interpolation_coefficients(std::move(interpolation_coefficients_params));

    odes::adams_extrapolation_coefficients_params_t extrapolation_coefficients_params {
        .order = order, .integrator = std::make_unique<odes::boole_integrator>(integrator_params)
    };
    odes::adams_extrapolation_coefficients extrapolation_coefficients(std::move(extrapolation_coefficients_params));

    odes::jacoby_mattrix_evaluator_params_t jacoby_mattrix_evaluator_params { .step = 1e100 };
    odes::uptr<odes::ijacoby_matrix_evaluator> matrix_evaluator
        = std::make_unique<odes::fourth_order_jacoby_mattrix_evaluator>(jacoby_mattrix_evaluator_params);

    odes::newton_root_finder_params_t newton_root_finder_params { .max_interations = 100,
                                                                  .precision       = 1e-45,
                                                                  .jacoby_matrix_evaluator
                                                                  = std::move(matrix_evaluator),
                                                                  .matrix_solver
                                                                  = std::make_unique<odes::full_pivot_gauss_solver>() };
    odes::uptr<odes::iroot_finder> root_finder
        = std::make_unique<odes::newton_root_finder>(std::move(newton_root_finder_params));

    odes::adams_interpolation_solver_params_t interpolation_solver_params {
        .order = order, .coefficients = std::move(interpolation_coefficients), .root_finder = std::move(root_finder)
    };

    odes::adams_extrapolation_solver_params_t extrapolation_solver_params { .order = order,
                                                                            .coefficients
                                                                            = std::move(extrapolation_coefficients) };

    odes::vector_t x0(4);
    x0[0] = 1.0;
    x0[1] = 0.0;
    x0[2] = 0.0;
    x0[3] = 0.5;

    odes::real_t r = sqrt(x0[0] * x0[0] + x0[1] * x0[1]);
    odes::real_t v = sqrt(x0[2] * x0[2] + x0[3] * x0[3]);
    odes::real_t t = calc_orbital_period(r, v);

    odes::ode_params_t ode_params { .t0 = 0.0, .t1 = t, .dt = 0.0001, .x0 = x0, .ode = ode };

    odes::adams_predictor_corrector_solver_params_t solver_params {
        .interpolation = std::move(interpolation_solver_params), .extrapolation = std::move(extrapolation_solver_params)
    };

    std::unique_ptr<odes::isolver> solver
        = std::make_unique<odes::adams_predictor_corrector_solver>(ode_params, std::move(solver_params));

    odes::integer_t n             = static_cast<odes::integer_t>((ode_params.t1 - ode_params.t0) / ode_params.dt);
    odes::integer_t print_every_n = std::max(odes::integer_t(1), n / 1000);

    odes::solve(std::move(solver), odes::array_t<odes::string_t> { "t", "x", "y", "vx", "vy" }, print_every_n);

    return 0;
}
