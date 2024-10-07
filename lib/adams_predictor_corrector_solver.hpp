#pragma once

#include "adams_extrapolation_solver.hpp"
#include "adams_interpolation_solver.hpp"
#include "solver_interface.hpp"

namespace odes {

struct adams_predictor_corrector_solver_params_t {
    adams_interpolation_solver_params_t interpolation;
    adams_extrapolation_solver_params_t extrapolation;
};

class adams_predictor_corrector_solver : public isolver {
public:
    explicit adams_predictor_corrector_solver(
        ode_params_t ode_params, adams_predictor_corrector_solver_params_t params);

    void step() noexcept override;
    const vector_t& current() const noexcept override final;
    const real_t& current_time() const noexcept override final;
    const real_t& begin_time() const noexcept override final;
    const real_t& end_time() const noexcept override final;

private:
    void compute_initial_values();

    ode_params_t ode_params_;
    adams_predictor_corrector_solver_params_t params_;
    array_t<vector_t> initial_;
    vector_t x_;
    real_t t_;
    multi_function_t solved_func_;
};

}
