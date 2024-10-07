#include "adams_extrapolation_solver.hpp"
#include "adams_interpolation_solver.hpp"
#include "adams_iterative_predictor_corrector_solver.hpp"

#include "boole_intergrator.hpp"
#include "rk4_solver.hpp"
#include "solver_interface.hpp"
#include "types.hpp"

#include <memory>

odes::real_t calc_orbital_period(odes::real_t r, odes::real_t v)
{
    odes::real_t a = odes::real_t(1) / (odes::real_t(2) / r - v * v);
    odes::real_t t = pow(a, odes::real_t(3) / odes::real_t(2)) * odes::real_t(2) * odes::real_t(M_PI);
    return t;
};

inline void
solve(std::unique_ptr<odes::isolver> solver, odes::array_t<odes::string_t> names, odes::integer_t print_every_n = 1)
{
    using namespace odes;

    integer_t precision  = 25;
    integer_t addiaional = 5;

    std::cout << names << std::endl;
    std::cout << std::setw(precision + addiaional) << std::fixed << std::setprecision(precision)
              << solver->current_time() << solver->current() << std::endl;

    integer_t counter = 0;
    while (solver->current_time() < solver->end_time()) {
        solver->step();

        if (++counter % print_every_n == 0) {
            std::cout << std::setw(precision + addiaional) << std::fixed << std::setprecision(precision)
                      << solver->current_time() << solver->current() << std::endl;
        }
    }

    std::cout << std::setw(precision + addiaional) << std::fixed << std::setprecision(precision)
              << solver->current_time() << solver->current() << std::endl;
}

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

    odes::vector_t x0(4);
    x0[0] = 1.0;
    x0[1] = 0.0;
    x0[2] = 0.0;
    x0[3] = 0.5;

    odes::real_t r  = sqrt(x0[0] * x0[0] + x0[1] * x0[1]);
    odes::real_t v  = sqrt(x0[2] * x0[2] + x0[3] * x0[3]);
    odes::real_t t  = calc_orbital_period(r, v);
    odes::real_t dt = 0.001;

    odes::real_t dt_norm = t / static_cast<odes::integer_t>(t / dt);

    odes::ode_params_t ode_params { .t0 = 0.0, .t1 = t, .dt = dt_norm, .x0 = x0, .ode = ode };

    odes::integer_t order = 4;
    odes::boole_integrator_params_t integrator_params { .order = 1'000'000 };

    odes::adams_interpolation_coefficients interpolation_coefficients(odes::adams_interpolation_coefficients_params_t {
        .order = order, .integrator = std::make_unique<odes::boole_integrator>(integrator_params) });

    odes::adams_extrapolation_coefficients extrapolation_coefficients(odes::adams_extrapolation_coefficients_params_t {
        .order = order, .integrator = std::make_unique<odes::boole_integrator>(integrator_params) });

    odes::adams_iterative_predictor_corrector_solver_params_t solver_params {
        .order                      = order,
        .iterations                 = 3,
        .interpolation_coefficients = std::move(interpolation_coefficients),
        .extrapolation_coefficients = std::move(extrapolation_coefficients),
        .initial_solver             = std::make_unique<odes::rk4_solver>(ode_params, odes::rk4_solver_params_t {})
    };

    std::unique_ptr<odes::isolver> solver
        = std::make_unique<odes::adams_iterative_predictor_corrector_solver>(ode_params, std::move(solver_params));

    odes::integer_t n             = static_cast<odes::integer_t>((ode_params.t1 - ode_params.t0) / ode_params.dt);
    odes::integer_t print_every_n = std::max(odes::integer_t(1), n / 1000);

    solve(std::move(solver), odes::array_t<odes::string_t> { "t", "x", "y", "vx", "vy" }, print_every_n);

    return 0;
}
