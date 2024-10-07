#pragma once

#include "integrator_interface.hpp"
#include "root_finder_interface.hpp"
#include "solver_interface.hpp"
#include "types.hpp"

namespace odes {

struct adams_interpolation_coefficients_params_t {
    integer_t order;
    uptr<iintegrator> integrator;
};

class adams_interpolation_coefficients {
public:
    explicit adams_interpolation_coefficients(adams_interpolation_coefficients_params_t params);
    const real_t& operator[](integer_t num) const;

private:
    void compute_polinomial();

    adams_interpolation_coefficients_params_t params_;
    array_t<real_t> a_;
};

struct adams_interpolation_solver_params_t {
    integer_t order;
    adams_interpolation_coefficients coefficients;
    uptr<iroot_finder> root_finder;
};

class adams_interpolation_solver : public isolver {
public:
    adams_interpolation_solver(ode_params_t ode_params, adams_interpolation_solver_params_t params);

    void step() noexcept override;
    const vector_t& current() const noexcept override final;
    const real_t& current_time() const noexcept override final;
    const real_t& begin_time() const noexcept override final;
    const real_t& end_time() const noexcept override final;

private:
    void compute_initial_values();

    ode_params_t ode_params_;
    adams_interpolation_solver_params_t params_;
    array_t<vector_t> initial_;
    vector_t x_;
    real_t t_;
    multi_function_t solved_func_;
};

} // namespace odes
