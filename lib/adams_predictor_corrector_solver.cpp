#include "adams_predictor_corrector_solver.hpp"

#include "rk4_solver.hpp"
#include "types.hpp"

namespace odes {

adams_predictor_corrector_solver::adams_predictor_corrector_solver(
    ode_params_t ode_params, adams_predictor_corrector_solver_params_t params)
    : ode_params_(ode_params)
    , params_(std::move(params))
    , x_(ode_params.x0)
{
    compute_initial_values();

    solved_func_ = [this](vector_t y) -> vector_t {
        vector_t res;

        res = params_.interpolation_coefficients[0] * ode_params_.ode(t_, y) * ode_params_.dt;
        for (integer_t j = 0; j < params_.order; ++j) {
            res += params_.interpolation_coefficients[j + 1] * initial_[params_.order - 1 - j] * ode_params_.dt;
        }

        res += x_ - y;

        return res;
    };
}

void adams_predictor_corrector_solver::compute_initial_values()
{
    for (size_t i = 0; i < params_.order; ++i) {
        params_.initial_solver->step();
        initial_.push_back(ode_params_.ode(params_.initial_solver->current_time(), params_.initial_solver->current()));
    }

    x_ = params_.initial_solver->current();
    t_ = params_.initial_solver->current_time();
}

void adams_predictor_corrector_solver::step() noexcept
{
    // Predictor step
    vector_t x = x_;
    for (size_t j = 0; j < params_.order; ++j) {
        x += params_.extrapolation_coefficients[params_.order - 1 - j] * initial_[j] * ode_params_.dt;
    }

    // Corrector step
    x_ = params_.root_finder->solve(solved_func_, x);

    t_ += ode_params_.dt;

    if (!params_.root_finder->precision_was_reached()) {
        std::cerr << "Residual is bugger than expected: " << params_.root_finder->get_residual() << std::endl;
    }

    std::shift_left(initial_.begin(), initial_.end(), 1);
    initial_.back() = ode_params_.ode(t_, x_);
}

const vector_t& adams_predictor_corrector_solver::current() const noexcept
{
    return x_;
}

const real_t& adams_predictor_corrector_solver::current_time() const noexcept
{
    return t_;
}

const real_t& adams_predictor_corrector_solver::begin_time() const noexcept
{
    return ode_params_.t0;
}

const real_t& adams_predictor_corrector_solver::end_time() const noexcept
{
    return ode_params_.t1;
}

}
