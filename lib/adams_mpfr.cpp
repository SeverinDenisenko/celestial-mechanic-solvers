#include "adams_mpfr.hpp"

namespace odes {

adams_solver_mpfr::adams_solver_mpfr(ode_params_t ode_params, adams_solver_mpfr_params_t adams_solver_mpfr_params)
    : ode_params_(ode_params)
    , adams_solver_mpfr_params_(adams_solver_mpfr_params)
{
}

void adams_solver_mpfr::step() noexcept
{
    // TODO
}

const vector_t& adams_solver_mpfr::current() const noexcept
{
    return current_;
}

const real_t& adams_solver_mpfr::current_time() const noexcept
{
    return time_;
}

const real_t& adams_solver_mpfr::begin_time() const noexcept
{
    return ode_params_.t0;
}

const real_t& adams_solver_mpfr::end_time() const noexcept
{
    return ode_params_.t1;
}

}
