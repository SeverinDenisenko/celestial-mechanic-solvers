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
    ode_params_t ode_params_;
    adams_solver_params_t adams_solver_params_;
    vector_t x_;
    real_t t_;
};

} // namespace odes
