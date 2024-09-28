#include "adams_extrapolation_solver.hpp"
#include "integrator_interface.hpp"
#include "rk4_solver.hpp"
#include "simpson_integrator.hpp"
#include "solver_interface.hpp"
#include "types.hpp"

#include <cmath>

using namespace odes;

class vector_to_system_mapper {
public:
    vector_to_system_mapper(vector_t& vec)
        : vec_(vec)
        , c_(vec_.size() / 4)
    {
    }

    void set_x(integer_t n, real_t x)
    {
        vec_[n] = x;
    }

    void set_y(integer_t n, real_t x)
    {
        vec_[n + c_] = x;
    }

    void set_vx(integer_t n, real_t x)
    {
        vec_[n + 2 * c_] = x;
    }

    void set_vy(integer_t n, real_t x)
    {
        vec_[n + 3 * c_] = x;
    }

    real_t get_x(integer_t n)
    {
        return vec_[n];
    }

    real_t get_y(integer_t n)
    {
        return vec_[n + c_];
    }

    real_t get_vx(integer_t n)
    {
        return vec_[n + 2 * c_];
    }

    real_t get_vy(integer_t n)
    {
        return vec_[n + 3 * c_];
    }

    real_t get_r(integer_t a, integer_t b)
    {
        return pow(pow(get_x(a) - get_x(b), 2) + pow(get_y(a) - get_y(b), 2), 0.5);
    }

    real_t get_v(integer_t n)
    {
        return pow(pow(get_vx(n), 2) + pow(get_vy(n), 2), 0.5);
    }

private:
    vector_t& vec_;
    integer_t c_;
};

vector_t create_x0()
{
    integer_t n = 3;
    vector_t x0(n * 4);
    vector_to_system_mapper mapper(x0);

    mapper.set_x(0, -0.07880227334416882);
    mapper.set_y(0, 0.5570371897354746);
    mapper.set_vx(0, 0.15998292728488323);
    mapper.set_vy(0, 1.1593418791674066);

    mapper.set_x(1, 0.5940359608209828);
    mapper.set_y(1, 0.383319210563721);
    mapper.set_vx(1, -0.5557289806160467);
    mapper.set_vy(1, -0.9029539156799118);

    mapper.set_x(2, -0.5152336874768139);
    mapper.set_y(2, -0.9403564002991956);
    mapper.set_vx(2, 0.39574605333116347);
    mapper.set_vy(2, -0.2563879634874948);

    return x0;
}

int main()
{
    size_t digits = 15;
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));

    integer_t n = 3;

    odes::ode_t ode = [n](odes::real_t t, odes::vector_t x) -> odes::vector_t {
        odes::vector_t y(x.size());
        vector_to_system_mapper ym(y);
        vector_to_system_mapper xm(x);

        for (integer_t i = 0; i < n; i++) {
            ym.set_x(i, xm.get_vx(i));
            ym.set_y(i, xm.get_vy(i));
        }

        for (integer_t i = 0; i < n; i++) {
            real_t vx = 0;
            real_t vy = 0;
            for (integer_t j = 0; j < n; ++j) {
                if (i == j) {
                    continue;
                }
                vx -= (xm.get_x(i) - xm.get_x(j)) / pow(xm.get_r(i, j), 3);
                vy -= (xm.get_y(i) - xm.get_y(j)) / pow(xm.get_r(i, j), 3);
            }
            ym.set_vx(i, vx);
            ym.set_vy(i, vy);
        }

        return y;
    };

    odes::integer_t order = 3;
    odes::simpson_integrator_params_t integrator_params { .order = 10'000'000 };
    odes::uptr<odes::iintegrator> integrator = std::make_unique<odes::simpson_integrator>(integrator_params);

    odes::adams_extrapolation_coefficients_params_t coefficients_params { .order      = order,
                                                                          .integrator = std::move(integrator) };
    odes::adams_extrapolation_coefficients coefficients(std::move(coefficients_params));

    odes::adams_extrapolation_solver_params_t solver_params { .order = order, .coefficients = std::move(coefficients) };

    odes::vector_t x0 = create_x0();

    odes::ode_params_t ode_params { .t0 = 0.0, .t1 = M_PI * 4, .dt = 0.00001, .x0 = x0, .ode = ode };

    std::unique_ptr<odes::isolver> solver
        = std::make_unique<odes::adams_extrapolation_solver>(ode_params, std::move(solver_params));

    odes::integer_t steps         = static_cast<odes::integer_t>((ode_params.t1 - ode_params.t0) / ode_params.dt);
    odes::integer_t print_every_n = std::max(odes::integer_t(1), steps / 1000);

    odes::solve(
        std::move(solver),
        odes::array_t<odes::string_t> {
            "t", "x1", "x2", "x3", "y1", "y2", "y3", "vx1", "vx2", "vx3", "vy1", "vy2", "vy3" },
        print_every_n);

    return 0;
}
