#include "adams_predictor_corrector_solver.hpp"

#include "rk4_solver.hpp"

namespace odes {

adams_predictor_corrector_solver::adams_predictor_corrector_solver(
    ode_params_t ode_params,
    adams_interpolation_solver_params_t interpolation_params,
    adams_extrapolation_solver_params_t extrapolation_params)
    : ode_params_(ode_params)
    , interpolation_params_(std::move(interpolation_params))
    , extrapolation_params_(std::move(extrapolation_params))
    , x_(ode_params.x0)
{
    compute_initial_values();

    solved_func_ = [this](vector_t y) -> vector_t {
        vector_t res;

        res = interpolation_params_.coefficients[0] * ode_params_.ode(t_, y) * ode_params_.dt;
        for (integer_t j = 0; j < interpolation_params_.order; ++j) {
            res += interpolation_params_.coefficients[j + 1] * initial_[interpolation_params_.order - 1 - j]
                * ode_params_.dt;
        }

        res += x_ - y;

        return res;
    };
}

void adams_predictor_corrector_solver::compute_initial_values()
{
    rk4_solver_params_t rk4_solver_params {};
    rk4_solver solver(ode_params_, rk4_solver_params);
    for (size_t i = 0; i < extrapolation_params_.order; ++i) {
        solver.step();
        initial_.push_back(ode_params_.ode(solver.current_time(), solver.current()));
    }

    x_ = solver.current();
    t_ = solver.current_time();
}

void adams_predictor_corrector_solver::step() noexcept
{
    // Predictor step
    for (size_t j = 0; j < extrapolation_params_.order; ++j) {
        x_ += extrapolation_params_.coefficients[extrapolation_params_.order - 1 - j] * initial_[j] * ode_params_.dt;
    }

    // Corrector step
    x_ = interpolation_params_.root_finder->solve(solved_func_, x_);

    t_ += ode_params_.dt;

    if (!interpolation_params_.root_finder->precision_was_reached()) {
        std::cerr << "Residual is bugger than expected: " << interpolation_params_.root_finder->get_residual()
                  << std::endl;
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
