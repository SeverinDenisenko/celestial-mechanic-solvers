#include "adams.hpp"
#include "rk.hpp"
#include "solver_interface.hpp"

int main()
{
    odes::ode_t ode = [](odes::real_t t, odes::vector_t x) -> odes::vector_t {
        odes::vector_t y(4);

        y[0] = x[0];
        y[1] = x[1];
        y[2] = x[2];
        y[3] = x[3];

        return y;
    };

    odes::adams_solver_params_t adams_solver_params { .order = 10 };
    odes::rk4_solver_params_t rk_solver_params {};

    odes::vector_t x0(4);
    x0[0] = 1.0;
    x0[1] = 0.0;
    x0[2] = 0.0;
    x0[3] = 0.5;

    odes::ode_params_t ode_params { .t0 = 0.0, .t1 = 10.0, .dt = 0.1, .x0 = x0, .ode = ode };

    std::unique_ptr<odes::isolver> solver = std::make_unique<odes::rk4_solver>(ode_params, rk_solver_params);

    odes::solve(std::move(solver));

    return 0;
}
