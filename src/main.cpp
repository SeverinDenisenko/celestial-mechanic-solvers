#include "adams.hpp"
#include "rk.hpp"
#include "solver_interface.hpp"

#include <cmath>

int main()
{
    size_t digits = 100;
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));

    odes::ode_t ode = [](odes::real_t t, odes::vector_t x) -> odes::vector_t {
        odes::vector_t y(4);

        y[0] = x[2];
        y[1] = x[3];

        odes::real_t r = sqrt(x[0] * x[0] + x[1] * x[1]);
        y[2]           = -x[0] / pow(r, 3);
        y[3]           = -x[1] / pow(r, 3);

        return y;
    };

    odes::adams_solver_params_t adams_solver_params { .order = 1 };
    odes::rk4_solver_params_t rk_solver_params {};

    odes::vector_t x0(4);
    x0[0] = 1.0;
    x0[1] = 0.0;
    x0[2] = 0.0;
    x0[3] = 0.5;

    odes::ode_params_t ode_params { .t0 = 0.0, .t1 = 1.0, .dt = 0.1, .x0 = x0, .ode = ode };

    std::unique_ptr<odes::isolver> solver = std::make_unique<odes::adams_solver>(ode_params, adams_solver_params);

    odes::solve(std::move(solver), odes::array_t<odes::string_t> { "t", "x", "y", "vx", "vy" });

    return 0;
}
