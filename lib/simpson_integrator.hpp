#pragma once

#include "integrator_interface.hpp"
#include "types.hpp"

namespace odes {

struct simpson_integrator_params_t {
    integer_t order;
};

class simpson_integrator : public iintegrator {
public:
    simpson_integrator(simpson_integrator_params_t params);
    real_t integrate(function_t& func, real_t a, real_t b) noexcept override;

private:
    simpson_integrator_params_t params_;
};

}
