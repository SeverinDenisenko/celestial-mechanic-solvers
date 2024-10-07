#pragma once

#include "integrator_interface.hpp"
#include "types.hpp"

namespace odes {

struct boole_integrator_params_t {
    integer_t order;
};

class boole_integrator : public iintegrator {
public:
    explicit boole_integrator(boole_integrator_params_t params);
    real_t integrate(function_t& func, real_t a, real_t b) noexcept override;

private:
    boole_integrator_params_t params_;
};

}
