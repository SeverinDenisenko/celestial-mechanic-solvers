#pragma once

#include "types.hpp"

namespace odes {

inline std::ostream& operator<<(std::ostream& stream, const vector_t& vector)
{
    for (size_t i = 0; i < vector.size(); ++i) {
        stream << " " << vector[i];
    }

    return stream;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& stream, const array_t<T>& vector)
{
    for (size_t i = 0; i < vector.size(); ++i) {
        stream << " " << vector[i];
    }

    return stream;
}

struct ode_params_t {
    real_t t0;
    real_t t1;
    real_t dt;
    vector_t x0;
    ode_t ode;
};

class isolver {
public:
    virtual void step() noexcept                        = 0;
    virtual const vector_t& current() const noexcept    = 0;
    virtual const real_t& current_time() const noexcept = 0;
    virtual const real_t& begin_time() const noexcept   = 0;
    virtual const real_t& end_time() const noexcept     = 0;
    virtual ~isolver()                                  = default;
};

inline void solve(std::unique_ptr<isolver> solver, array_t<string_t> names, odes::integer_t print_every_n = 1)
{
    std::cout << names << std::endl;
    std::cout << std::setw(12) << std::fixed << std::setprecision(9) << solver->current_time() << solver->current()
              << std::endl;
    odes::integer_t counter = 0;
    while (solver->current_time() < solver->end_time()) {
        solver->step();
        if (++counter % print_every_n == 0) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(14) << solver->current_time()
                      << solver->current() << std::endl;
        }
    }
    std::cout << std::setw(12) << std::fixed << std::setprecision(14) << solver->current_time()
              << solver->current() << std::endl;
}

} // namespace odes
