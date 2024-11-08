#pragma once

#include "nbody_solver_interface.hpp"
#include "types.hpp"

namespace odes {

class hermite_solver : public inbody_solver {
public:
    hermite_solver(nbody_params_t params);

    void step() noexcept override;
    const nbody_t& current() const noexcept override;
    const real_t& current_time() const noexcept override;
    const real_t& begin_time() const noexcept override;
    const real_t& end_time() const noexcept override;

private:
    void calculate_a_and_j(nbody_t& object);

    nbody_params_t params_;
    nbody_t bodies_;
    nbody_t bodies_next_;
    real_t t_;
};

}
