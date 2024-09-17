#pragma once

#include "solver_interface.hpp"

namespace odes {

struct adams_solver_params_t {
    integer_t order;
};

class adams_solver : public isolver {
public:
    adams_solver(ode_params_t ode_params, adams_solver_params_t adams_solver_params);

    void step() noexcept override;
    const vector_t& current() const noexcept override final;
    const real_t& current_time() const noexcept override final;
    const real_t& begin_time() const noexcept override final;
    const real_t& end_time() const noexcept override final;

private:
    real_t integrate(std::function<real_t(size_t, real_t)> integrand, size_t j);
    size_t factorial(size_t x);

    ode_params_t ode_params_;
    adams_solver_params_t adams_solver_params_;
    array_t<vector_t> initial_;
    vector_t x_;
    real_t t_;
    array_t<real_t> a_;
};

} // namespace odes
