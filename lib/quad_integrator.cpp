#include "quad_integrator.hpp"

namespace odes {

quad_integrator::quad_integrator(quad_integrator_params_t quad_integrator_params)
	: quad_integrator_params_(quad_integrator_params)
{
}

real_t quad_integrator::integrate(function_t& func, real_t a, real_t b) noexcept
{
    real_t res = 0;
    real_t h = (b - a) / real_t(quad_integrator_params_.order);

    for (integer_t i = 0; i < quad_integrator_params_.order; ++i) {
        real_t x = a + (b - a) * real_t(i) / real_t(quad_integrator_params_.order);
        res += func(x) * h;
    }

    return res;
}

}
