#include "adams_iterative_predictor_corrector_solver.hpp"

#include "rk4_solver.hpp"
#include "types.hpp"

namespace odes {

adams_iterative_predictor_corrector_solver::adams_iterative_predictor_corrector_solver(
    ode_params_t ode_params, adams_iterative_predictor_corrector_solver_params_t params)
    : ode_params_(ode_params)
    , params_(std::move(params))
    , x_(ode_params.x0)
{
    compute_initial_values();
}

void adams_iterative_predictor_corrector_solver::compute_initial_values()
{
    rk4_solver_params_t rk4_solver_params {};
    rk4_solver solver(ode_params_, rk4_solver_params);
    for (size_t i = 0; i < params_.order; ++i) {
        solver.step();
        initial_.push_back(ode_params_.ode(solver.current_time(), solver.current()));
    }

    x_ = solver.current();
    t_ = solver.current_time();
}

void adams_iterative_predictor_corrector_solver::step() noexcept
{
    // Predictor step
    vector_t x = x_;
    for (size_t j = 0; j < params_.order; ++j) {
        x += params_.extrapolation_coefficients[params_.order - 1 - j] * initial_[j] * ode_params_.dt;
    }

    // Corrector steps
    for (integer_t i = 0; i < params_.iterations; ++i) {
        x = x_ + params_.interpolation_coefficients[0] * ode_params_.ode(t_, x) * ode_params_.dt;
        for (integer_t j = 0; j < params_.order; ++j) {
            x += params_.interpolation_coefficients[j + 1] * initial_[params_.order - 1 - j] * ode_params_.dt;
        }
    }

    x_ = x;
    t_ += ode_params_.dt;

    std::shift_left(initial_.begin(), initial_.end(), 1);
    initial_.back() = ode_params_.ode(t_, x_);
}

const vector_t& adams_iterative_predictor_corrector_solver::current() const noexcept
{
    return x_;
}

const real_t& adams_iterative_predictor_corrector_solver::current_time() const noexcept
{
    return t_;
}

const real_t& adams_iterative_predictor_corrector_solver::begin_time() const noexcept
{
    return ode_params_.t0;
}

const real_t& adams_iterative_predictor_corrector_solver::end_time() const noexcept
{
    return ode_params_.t1;
}

}
