#include "boole_intergrator.hpp"

#include <stdexcept>

#include "types.hpp"

namespace odes {

boole_integrator::boole_integrator(boole_integrator_params_t params)
    : params_(params)
{
    if (params_.order % 4 != 0) {
        throw std::runtime_error("order of boole_integrator must be divisible by 4");
    }
}

real_t boole_integrator::integrate(function_t& func, real_t a, real_t b) noexcept
{
    real_t res;
    real_t sp;
    real_t s1 = 0, s2 = 0, s3 = 0;
    real_t h = (b - a) / real_t(params_.order);

    for (integer_t i = 1; i <= params_.order - 1; i += 2) {
        s1 += func(h * i);
    }

    for (integer_t i = 2; i <= params_.order - 2; i += 4) {
        s2 += func(h * i);
    }

    for (integer_t i = 4; i <= params_.order - 4; i += 4) {
        s3 += func(h * i);
    }

    sp  = (func(a) + func(b)) * 7 + 32 * s1 + 12 * s2 + 14 * s3;
    res = sp * (2 * h) / 45;

    return res;
}

}
