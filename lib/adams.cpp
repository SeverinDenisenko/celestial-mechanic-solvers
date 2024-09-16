#include "adams.hpp"

namespace odes {

adams_solver::adams_solver(ode_params_t ode_params, adams_solver_params_t adams_solver_params)
    : ode_params_(ode_params)
    , adams_solver_params_(adams_solver_params)
{
}

void adams_solver::step() noexcept
{
    t_ += ode_params_.dt;
}

const vector_t& adams_solver::current() const noexcept
{
    return x_;
}

const real_t& adams_solver::current_time() const noexcept
{
    return t_;
}

const real_t& adams_solver::begin_time() const noexcept
{
    return ode_params_.t0;
}

const real_t& adams_solver::end_time() const noexcept
{
    return ode_params_.t1;
}

}
