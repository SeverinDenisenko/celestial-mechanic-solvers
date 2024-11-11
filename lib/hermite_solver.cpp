#include "hermite_solver.hpp"
#include "nbody_solver_interface.hpp"
#include "types.hpp"

namespace odes {

hermite_solver::hermite_solver(nbody_params_t params)
    : params_(std::move(params))
    , bodies_(std::move(params_.state))
    , bodies_next_(bodies_)
    , t_(params_.t0)
{
    calculate_a_and_j(bodies_);
}

void hermite_solver::calculate_a_and_j(nbody_t& target)
{
    for (integer_t i = 0; i < target.size(); ++i) {
        vector_t a(target[i].a.size());
        vector_t j(target[i].j.size());

        for (integer_t k = 0; k < target.size(); ++k) {
            if (i == k) {
                continue;
            }

            vector_t r = target[i].r - target[k].r;
            vector_t v = target[i].v - target[k].v;

            a += params_.a(r, v, params_.m[k]);
            j += params_.j(r, v, params_.m[k]);
        }

        target[i].a = a;
        target[i].j = j;
    }
}

void hermite_solver::step() noexcept
{
    real_t dt = params_.dt;

    real_t dt2 = dt * dt;
    real_t dt3 = dt2 * dt;
    real_t dt4 = dt3 * dt;
    real_t dt5 = dt4 * dt;

    for (integer_t i = 0; i < bodies_.size(); ++i) {
        bodies_next_[i].r = bodies_[i].r + bodies_[i].v * dt + real_t(1) / real_t(2) * bodies_[i].a * dt2
            + real_t(1) / real_t(6) * bodies_[i].j * dt3;
        bodies_next_[i].v = bodies_[i].v + bodies_[i].a * dt + real_t(1) / real_t(2) * bodies_[i].j * dt2;
    }

    calculate_a_and_j(bodies_next_);

    for (integer_t i = 0; i < bodies_.size(); ++i) {
        vector_t a2 = (-6 * (bodies_[i].a - bodies_next_[i].a) - dt * (4 * bodies_[i].j + 2 * bodies_next_[i].j)) / dt2;
        vector_t a3 = (12 * (bodies_[i].a - bodies_next_[i].a) + 6 * dt * (bodies_[i].j + bodies_next_[i].j)) / dt3;
        bodies_[i].r = bodies_next_[i].r + dt4 * a2 / 24 + dt5 * a3 / 120;
        bodies_[i].v = bodies_next_[i].v + dt3 * a2 / 6 + dt4 * a3 / 24;
    }

    calculate_a_and_j(bodies_);

    t_ += dt;
}

const nbody_t& hermite_solver::current() const noexcept
{
    return bodies_;
}

const real_t& hermite_solver::current_time() const noexcept
{
    return t_;
}

const real_t& hermite_solver::begin_time() const noexcept
{
    return params_.t0;
}

const real_t& hermite_solver::end_time() const noexcept
{
    return params_.t1;
}

}
