#include "adams.hpp"

#include <cmath>

#include "rk.hpp"

namespace odes {

adams_solver::adams_solver(ode_params_t ode_params, adams_solver_params_t adams_solver_params)
    : ode_params_(ode_params)
    , adams_solver_params_(adams_solver_params)
    , x_(ode_params.x0)
{
    // Compute initial values
    rk4_solver_params_t rk4_solver_params {};
    rk4_solver solver(ode_params_, rk4_solver_params);
    for (size_t i = 0; i < adams_solver_params_.order; ++i) {
        solver.step();
        initial_.push_back(solver.current());
    }
    x_ = solver.current();
    t_ = solver.current_time();

    // Compute polinomial
    std::function<real_t(size_t, real_t)> integrand = [this](size_t j, real_t z) -> real_t {
        real_t res = 1;

        for (size_t i = 0; i < adams_solver_params_.order; ++i) {
            if (i == j) {
                continue;
            }
            res *= z + i;
        }

        return res;
    };

    for (size_t j = 0; j < adams_solver_params_.order; ++j) {
        real_t a = pow(-1.0, j) / (factorial(j) * factorial(adams_solver_params_.order - 1 - j));
        a *= integrate(integrand, j);
        a_.push_back(a);
    }
}

real_t adams_solver::integrate(std::function<real_t(size_t, real_t)> integrand, size_t j)
{
    size_t order = 100'000;
    real_t res   = 0;

    for (size_t i = 0; i < order; ++i) {
        res += integrand(j, (real_t)(i + 1) / order) / order;
    }

    return res;
}

size_t adams_solver::factorial(size_t x)
{
    size_t res = 1;

    for (size_t i = 2; i <= x; ++i) {
        res *= i;
    }

    return res;
}

void adams_solver::step() noexcept
{
    for (size_t j = 0; j < adams_solver_params_.order; ++j) {
        x_ += a_[adams_solver_params_.order - 1 - j] * initial_[j] * ode_params_.dt;
    }
    t_ += ode_params_.dt;

    // Shift the vector of intials
    for (size_t i = 0; i < adams_solver_params_.order - 1; ++i) {
        initial_[i] = std::move(initial_[i + 1]);
    }
    initial_[adams_solver_params_.order - 1] = ode_params_.ode(t_, x_);
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
