#pragma once

#include "types.hpp"

namespace odes {

class iintegrator {
public:
	virtual real_t integrate(function_t& func, real_t a, real_t b) noexcept = 0;
};

}
