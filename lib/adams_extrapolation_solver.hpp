#pragma once

#include "integrator_interface.hpp"
#include "solver_interface.hpp"
#include "types.hpp"

namespace odes {

struct adams_extrapolation_coefficients_params_t {
    integer_t order;
    uptr<iintegrator> integrator;
};

class adams_extrapolation_coefficients {
public:
    adams_extrapolation_coefficients(adams_extrapolation_coefficients_params_t params);
    const real_t& operator[](integer_t num) const;

private:
    void compute_polinomial();
    integer_t factorial(integer_t x);

    adams_extrapolation_coefficients_params_t params_;
    array_t<real_t> a_;
};

struct adams_extrapolation_solver_params_t {
    integer_t order;
    adams_extrapolation_coefficients coefficients;
};

class adams_extrapolation_solver : public isolver {
public:
    adams_extrapolation_solver(
        ode_params_t ode_params, adams_extrapolation_solver_params_t params);

    void step() noexcept override;
    const vector_t& current() const noexcept override final;
    const real_t& current_time() const noexcept override final;
    const real_t& begin_time() const noexcept override final;
    const real_t& end_time() const noexcept override final;

private:
    void compute_initial_values();

    ode_params_t ode_params_;
    adams_extrapolation_solver_params_t params_;
    array_t<vector_t> initial_;
    vector_t x_;
    real_t t_;
};

} // namespace odes
