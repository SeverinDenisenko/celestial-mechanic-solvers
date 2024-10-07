#pragma once

#include "integrator_interface.hpp"
#include "types.hpp"

namespace odes {

struct quad_integrator_params_t {
    integer_t order;
};

class quad_integrator : public iintegrator {
public:
    explicit quad_integrator(quad_integrator_params_t quad_integrator_params);

    real_t integrate(function_t& func, real_t a, real_t b) noexcept override;

private:
    quad_integrator_params_t quad_integrator_params_;
};

}
