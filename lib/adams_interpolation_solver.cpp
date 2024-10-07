#include "adams_interpolation_solver.hpp"

#include "rk4_solver.hpp"
#include "types.hpp"
#include "utils.hpp"
#include <algorithm>

namespace odes {
adams_interpolation_coefficients::adams_interpolation_coefficients(adams_interpolation_coefficients_params_t params)
    : params_(std::move(params))
{
    compute_polinomial();
}

const real_t& adams_interpolation_coefficients::operator[](integer_t num) const
{
    return a_[num];
}

void adams_interpolation_coefficients::compute_polinomial()
{
    for (integer_t j = 0; j < params_.order + 1; ++j) {
        function_t integrand = [this, j](real_t z) -> real_t {
            real_t res = 1;

            for (integer_t i = 0; i < params_.order + 1; ++i) {
                if (i == j) {
                    continue;
                }
                res *= z + real_t(i) - real_t(1);
            }

            return res;
        };

        real_t a;
        a = params_.integrator->integrate(integrand, real_t(0), real_t(1));
        a /= real_t(factorial(j) * factorial(params_.order - j));
        a *= j % 2 == 0 ? real_t(1.0) : real_t(-1.0);
        a_.push_back(a);
    }
}

adams_interpolation_solver::adams_interpolation_solver(
    ode_params_t ode_params, adams_interpolation_solver_params_t params)
    : ode_params_(ode_params)
    , params_(std::move(params))
    , x_(ode_params.x0)
{
    compute_initial_values();

    solved_func_ = [this](vector_t y) -> vector_t {
        vector_t res;

        res = params_.coefficients[0] * ode_params_.ode(t_, y) * ode_params_.dt;
        for (integer_t j = 0; j < params_.order; ++j) {
            res += params_.coefficients[j + 1] * initial_[params_.order - 1 - j] * ode_params_.dt;
        }

        res += x_ - y;

        return res;
    };
}

void adams_interpolation_solver::compute_initial_values()
{
    for (size_t i = 0; i < params_.order; ++i) {
        params_.initial_solver->step();
        initial_.push_back(ode_params_.ode(params_.initial_solver->current_time(), params_.initial_solver->current()));
    }

    x_ = params_.initial_solver->current();
    t_ = params_.initial_solver->current_time();
}

void adams_interpolation_solver::step() noexcept
{
    x_ = params_.root_finder->solve(solved_func_, x_);
    t_ += ode_params_.dt;

    if (!params_.root_finder->precision_was_reached()) {
        std::cerr << "Residual is bugger than expected: " << params_.root_finder->get_residual() << std::endl;
    }

    std::shift_left(initial_.begin(), initial_.end(), 1);
    initial_.back() = ode_params_.ode(t_, x_);
}

const vector_t& adams_interpolation_solver::current() const noexcept
{
    return x_;
}

const real_t& adams_interpolation_solver::current_time() const noexcept
{
    return t_;
}

const real_t& adams_interpolation_solver::begin_time() const noexcept
{
    return ode_params_.t0;
}

const real_t& adams_interpolation_solver::end_time() const noexcept
{
    return ode_params_.t1;
}

}
