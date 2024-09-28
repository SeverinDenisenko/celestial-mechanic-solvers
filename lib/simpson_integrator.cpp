#include "simpson_integrator.hpp"
#include "types.hpp"

namespace odes {

simpson_integrator::simpson_integrator(simpson_integrator_params_t params)
    : params_(params)
{
}

real_t simpson_integrator::integrate(function_t& func, real_t a, real_t b) noexcept
{
    real_t res;
    real_t s1 = 0, s2 = 0;
    real_t h = (b - a) / real_t(params_.order);

    for (integer_t i = 3; i < params_.order - 1; i += 2) {
        s2 = s2 + func(h * i);
    }

    for (integer_t i = 2; i < params_.order; i += 2) {
        s1 = s1 + func(h * i);
    }

    real_t sp  = func(a) + func(b) + 4 * s1 + 2 * s2;
    res = (h / 3) * sp;

    return res;
}

}
