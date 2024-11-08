#include "hermite_solver.hpp"
#include "nbody_solver_interface.hpp"
#include "types.hpp"

#include <complex>
#include <iomanip>
#include <memory>

odes::real_t calc_orbital_period(odes::real_t r, odes::real_t v)
{
    odes::real_t a = odes::real_t(1) / (odes::real_t(2) / r - v * v);
    odes::real_t t = pow(a, odes::real_t(3) / odes::real_t(2)) * odes::real_t(2) * odes::real_t(M_PI);
    return t;
};

inline void solve(std::unique_ptr<odes::inbody_solver> solver, odes::array_t<odes::string_t> names)
{
    using namespace odes;

    integer_t precision  = 25;
    integer_t addiaional = 5;

    std::cout << names << std::endl;
    std::cout << std::setw(precision + addiaional) << std::fixed << std::setprecision(precision)
              << solver->current_time() << solver->current()[1].r - solver->current()[0].r
              << solver->current()[1].v - solver->current()[0].v << std::endl;

    while (solver->current_time() < solver->end_time()) {
        std::cout << std::setw(precision + addiaional) << std::fixed << std::setprecision(precision)
                  << solver->current_time() << solver->current()[1].r - solver->current()[0].r
                  << solver->current()[1].v - solver->current()[0].v << std::endl;

        solver->step();
    }

    std::cout << std::setw(precision + addiaional) << std::fixed << std::setprecision(precision)
              << solver->current_time() << solver->current()[1].r - solver->current()[0].r
              << solver->current()[1].v - solver->current()[0].v << std::endl;
}

int main()
{
    using odes::real_t;
    using odes::vector_t;

    size_t digits = 50;
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));

    odes::array_t<real_t> m = { 0.5, 0.5 };

    odes::nbody_t bodies(2);
    bodies[0].r = odes::make_vector<real_t>({ 0.0, 0.0, 0.0 });
    bodies[0].v = odes::make_vector<real_t>({ 0.0, 0.0, 0.0 });

    bodies[1].r = odes::make_vector<real_t>({ 1.0, 0.0, 0.0 });
    bodies[1].v = odes::make_vector<real_t>({ 0.0, 0.5, 0.0 });

    real_t r  = 1.0;
    real_t v  = 0.5;
    real_t t  = calc_orbital_period(r, v);
    real_t dt = 0.0001;

    odes::interaction_function_t a = [](vector_t r, vector_t v, real_t m) -> vector_t {
        real_t r3 = pow(odes::norm(r), 3);
        return -r * m / r3;
    };

    odes::interaction_function_t j = [](vector_t r, vector_t v, real_t m) -> vector_t {
        real_t r3 = pow(odes::norm(r), 3);
        real_t r5 = pow(odes::norm(r), 5);
        return -m * (v / r3 - 3 * odes::dot(r, v) * r / r5);
    };

    odes::nbody_params_t params { .t0 = 0.0, .t1 = t, .dt = dt, .state = bodies, .m = m, .a = a, .j = j };

    std::unique_ptr<odes::inbody_solver> solver = std::make_unique<odes::hermite_solver>(params);

    solve(std::move(solver), odes::array_t<odes::string_t> { "t", "x", "y", "z", "vx", "vy", "vz" });

    return 0;
}
