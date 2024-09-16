#include "adams.hpp"
#include "solver_interface.hpp"

int main()
{
    odes::ode_t ode = [](odes::real_t t, odes::vector_t x) -> odes::vector_t { return {}; };

    odes::adams_solver_params_t solver_params { .order = 10 };
    odes::ode_params_t ode_params { .t0 = 0.0, .t1 = 10.0, .dt = 0.1, .ode = ode };

    std::unique_ptr<odes::isolver> solver = std::make_unique<odes::adams_solver>(ode_params, solver_params);

    odes::solve(std::move(solver));

    return 0;
}
