#pragma once

#include "solver_interface.hpp"

namespace odes {

struct adams_solver_mpfr_params_t {
    integer_t order;
};

class adams_solver_mpfr : public isolver {
public:
    adams_solver_mpfr(ode_params_t ode_params, adams_solver_mpfr_params_t adams_solver_mpfr_params);

    void step() noexcept override;
    const vector_t& current() const noexcept override final;
    const real_t& current_time() const noexcept override final;
    const real_t& begin_time() const noexcept override final;
    const real_t& end_time() const noexcept override final;

private:
    ode_params_t ode_params_;
    adams_solver_mpfr_params_t adams_solver_mpfr_params_;
    vector_t current_;
    real_t time_;
};

}
