#include "adams.hpp"
#include "rk.hpp"
#include "solver_interface.hpp"

#include <cmath>

odes::real_t calc_orbital_period(odes::real_t r, odes::real_t v)
{
    odes::real_t a = odes::real_t(1) / (odes::real_t(2) / r - v * v);
    odes::real_t t = pow(a, odes::real_t(3) / odes::real_t(2)) * odes::real_t(2) * mpfr::const_pi();
    return t;
};

int main()
{
    size_t digits = 20;
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

    odes::adams_solver_params_t adams_solver_params { .order = 4 };
    odes::rk4_solver_params_t rk_solver_params {};

    odes::vector_t x0(4);
    x0[0] = 1.0;
    x0[1] = 0.0;
    x0[2] = 0.0;
    x0[3] = 0.5;

    odes::real_t r = sqrt(x0[0] * x0[0] + x0[1] * x0[1]);
    odes::real_t v = sqrt(x0[2] * x0[2] + x0[3] * x0[3]);
    odes::real_t t = calc_orbital_period(r, v);

    odes::ode_params_t ode_params { .t0 = 0.0, .t1 = t, .dt = 0.001, .x0 = x0, .ode = ode };

    std::unique_ptr<odes::isolver> solver = std::make_unique<odes::adams_solver>(ode_params, adams_solver_params);

    odes::integer_t n             = static_cast<odes::integer_t>((ode_params.t1 - ode_params.t0) / ode_params.dt);
    odes::integer_t print_every_n = std::max(odes::integer_t(1), n / 1000);

    odes::solve(std::move(solver), odes::array_t<odes::string_t> { "t", "x", "y", "vx", "vy" }, print_every_n);

    return 0;
}
