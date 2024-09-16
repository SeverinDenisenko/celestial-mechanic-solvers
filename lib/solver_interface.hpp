#pragma once

#include <functional>
#include <iostream>
#include <iomanip>
#include <memory>

#include <boost/numeric/ublas/vector.hpp>

namespace odes {

namespace ublas = boost::numeric::ublas;

using real_t    = double;
using integer_t = int;
using vector_t  = ublas::vector<real_t>;
using ode_t     = std::function<vector_t(real_t, vector_t)>;

inline std::ostream& operator<<(std::ostream& stream, const vector_t& vector)
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

inline void solve(std::unique_ptr<isolver> solver)
{
    while (solver->current_time() < solver->end_time()) {
        solver->step();
        std::cout << std::setw(12) << std::fixed << std::setprecision(9) << solver->current_time() << solver->current() << std::endl;
    }
}

} // namespace odes
