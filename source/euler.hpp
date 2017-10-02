#ifndef JFLOW_EULER_HPP
#define JFLOW_EULER_HPP

#include "jflow.hpp"
#include <array>

namespace jflow {

// TODO: Factor out hard-coded gas model (currently pg, g=1.4)
struct euler {

    // Type definitions
    using state    = vector4;
    using flux     = vector4;
    using jacobian = matrix44;

    // Field indices
    struct field {
        static const std::size_t density      = 0;
        static const std::size_t momentum_x   = 1;
        static const std::size_t momentum_y   = 2;
        static const std::size_t total_energy = 3;
    };

    // Function for debug/test
    static auto compute_flux(const state& q, const vector2& n) -> flux;
    static auto compute_freestream_flux(const state& q, const vector2& n) -> flux;
    static auto compute_wall_flux(const state& q, const vector2& n) -> flux;
    static auto compute_jump_flux(const state& ql, const state& qr, const vector2& n) -> flux;
};

struct perfect_gas {
    static const double gamma;
    static auto compute_pressure(double e, double rho) -> double {
        return (gamma - 1) * rho * e;
    }
    static auto compute_sound_speed(double e, double rho) -> double {
        return std::sqrt(gamma * (gamma - 1) * e);
    }
};

}  // namespace jflow

#endif
