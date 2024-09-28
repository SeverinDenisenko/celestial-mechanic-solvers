#include "rk4_solver.hpp"

namespace odes {

rk4_solver::rk4_solver(ode_params_t ode_params, rk4_solver_params_t rk4_solver_params)
    : ode_params_(ode_params)
    , rk4_solver_params_(rk4_solver_params)
    , x_(ode_params_.x0)
{
}

vector_t rk4_solver::k1(real_t t, vector_t x) const noexcept
{
    return ode_params_.ode(t, x) * ode_params_.dt;
}

vector_t rk4_solver::k2(real_t t, vector_t x, vector_t k2) const noexcept
{
    return ode_params_.ode(t + ode_params_.dt / 2.0, x + k2 / 2.0) * ode_params_.dt;
}

vector_t rk4_solver::k3(real_t t, vector_t x, vector_t k3) const noexcept
{
    return ode_params_.ode(t + ode_params_.dt / 2.0, x + k3 / 2.0) * ode_params_.dt;
}

vector_t rk4_solver::k4(real_t t, vector_t x, vector_t k4) const noexcept
{
    return ode_params_.ode(t + ode_params_.dt, x + k4) * ode_params_.dt;
}

void rk4_solver::step() noexcept
{
    vector_t _k1 = k1(t_, x_);
    vector_t _k2 = k2(t_, x_, _k1);
    vector_t _k3 = k3(t_, x_, _k2);
    vector_t _k4 = k4(t_, x_, _k3);

    x_ += 1.0 / 6.0 * (_k1 + 2 * _k2 + 2 * _k3 + _k4);
    t_ += ode_params_.dt;
}

const vector_t& rk4_solver::current() const noexcept
{
    return x_;
}

const real_t& rk4_solver::current_time() const noexcept
{
    return t_;
}

const real_t& rk4_solver::begin_time() const noexcept
{
    return ode_params_.t0;
}

const real_t& rk4_solver::end_time() const noexcept
{
    return ode_params_.t1;
}

}
