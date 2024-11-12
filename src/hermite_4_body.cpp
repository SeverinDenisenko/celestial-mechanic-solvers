#include "hermite_solver.hpp"
#include "mpreal/mpreal.h"
#include "nbody_solver_interface.hpp"
#include "types.hpp"

#include <cmath>
#include <iomanip>
#include <memory>

void print(std::unique_ptr<odes::inbody_solver>& solver)
{
    using namespace odes;

    integer_t precision  = 25;
    integer_t addiaional = 5;

    std::cout << std::setw(precision + addiaional) << std::fixed << std::setprecision(precision)
              << solver->current_time() << solver->current()[0].r << solver->current()[0].v << solver->current()[1].r
              << solver->current()[1].v << solver->current()[2].r << solver->current()[2].v << std::endl;
}

inline void
solve(std::unique_ptr<odes::inbody_solver> solver, odes::array_t<odes::string_t> names, odes::integer_t print_n)
{
    using namespace odes;

    std::cout << names << std::endl;
    print(solver);

    odes::integer_t i = 0;
    while (solver->current_time() < solver->end_time()) {
        solver->step();
        ++i;

        if (i % print_n == 0) {
            print(solver);
        }
    }

    print(solver);
}

int main()
{
    using odes::integer_t;
    using odes::real_t;
    using odes::vector_t;

    size_t digits = 50;
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));

    odes::array_t<real_t> m = { 0.333, 0.333, 0.333 };

    odes::nbody_t bodies(3);
    bodies[0].r = odes::make_vector<real_t>({ 0.0, 0.8, 0.0 });
    bodies[0].v = odes::make_vector<real_t>({ -0.8, 0, 0.0 });
    bodies[1].r = odes::make_vector<real_t>({ 0.8 * cos(2.0 * M_PI / 3), 0.8 * sin(2.0 * M_PI / 3), 0.0 });
    bodies[1].v = odes::make_vector<real_t>({ 0.8 * sin(2.0 * M_PI / 3), 0.8 * cos(2.0 * M_PI / 3), 0.0 });
    bodies[2].r = odes::make_vector<real_t>({ 0.8 * cos(4.0 * M_PI / 3), 0.8 * sin(4.0 * M_PI / 3), 0.0 });
    bodies[2].v = odes::make_vector<real_t>({ 0.8 * sin(4.0 * M_PI / 3), 0.8 * cos(4.0 * M_PI / 3), 0.0 });

    real_t t          = 10 * M_PI;
    real_t dt         = 0.001;
    real_t dt_norm    = t / static_cast<integer_t>(t / dt);
    integer_t print_n = static_cast<integer_t>(t / dt_norm) / 1000;

    odes::interaction_function_t a = [](vector_t r, vector_t v, real_t m) -> vector_t {
        real_t r3 = pow(odes::norm(r), 3);
        return -r * m / r3;
    };

    odes::interaction_function_t j = [](vector_t r, vector_t v, real_t m) -> vector_t {
        real_t r3 = pow(odes::norm(r), 3);
        real_t r5 = pow(odes::norm(r), 5);
        return -m * (v / r3 - 3 * odes::dot(v, r) * r / r5);
    };

    odes::nbody_params_t params {
        .t0 = 0.0, .t1 = t, .dt = dt_norm, .count = 3, .m = m, .state = bodies, .a = a, .j = j
    };

    std::unique_ptr<odes::inbody_solver> solver = std::make_unique<odes::hermite_solver>(params);

    solve(
        std::move(solver),
        odes::array_t<odes::string_t> { "t",
                                        "x0",
                                        "y0",
                                        "z0",
                                        "vx0",
                                        "vy0",
                                        "vz0",
                                        "x1",
                                        "y1",
                                        "z1",
                                        "vx1",
                                        "vy1",
                                        "vz1",
                                        "x2",
                                        "y2",
                                        "z2",
                                        "vx2",
                                        "vy2",
                                        "vz2" },
        print_n);

    return 0;
}
