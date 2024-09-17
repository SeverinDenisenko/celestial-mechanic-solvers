#pragma once

#include "solver_interface.hpp"

namespace odes {

struct rk4_solver_params_t {
};

class rk4_solver : public isolver {
public:
    rk4_solver(ode_params_t ode_params, rk4_solver_params_t rk4_solver_params);

    void step() noexcept override;
    const vector_t& current() const noexcept override final;
    const real_t& current_time() const noexcept override final;
    const real_t& begin_time() const noexcept override final;
    const real_t& end_time() const noexcept override final;

private:
    vector_t k1(real_t t, vector_t x) const noexcept;
    vector_t k2(real_t t, vector_t x, vector_t k2) const noexcept;
    vector_t k3(real_t t, vector_t x, vector_t k3) const noexcept;
    vector_t k4(real_t t, vector_t x, vector_t k4) const noexcept;

    ode_params_t ode_params_;
    rk4_solver_params_t rk4_solver_params_;
    vector_t x_;
    real_t t_;
};

} // namespace odes
