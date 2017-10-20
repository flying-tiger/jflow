#ifndef JFLOW_EULER_HPP
#define JFLOW_EULER_HPP

#include "jflow/common.hpp"
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

    // Flux functions
    static auto compute_flux(const state& q, const vector2& n) -> flux;
    static auto compute_freestream_flux(const state& q, const vector2& n) -> flux;
    static auto compute_wall_flux(const state& q, const vector2& n) -> flux;
    static auto compute_jump_flux(const state& ql, const state& qr, const vector2& n) -> flux;

    // Global configuration
    static auto set_freestream(double p, double T, double vx, double vy) -> void;

  private:
    // TODO: Get rid of global static state
    static state freestream_;
};

struct perfect_gas {
    static auto compute_energy(double T) -> double {
        auto& R = perfect_gas::specific_gas_constant_;
        auto& g = perfect_gas::specific_heat_ratio_;
        return R * T / (g - 1);
    }
    static auto compute_density(double p, double T) -> double {
        auto& R = perfect_gas::specific_gas_constant_;
        return p / (R * T);
    }
    static auto compute_pressure(double e, double rho) -> double {
        auto& g = perfect_gas::specific_heat_ratio_;
        return (g - 1) * rho * e;
    }
    static auto compute_sound_speed(double e, double rho) -> double {
        auto& g = perfect_gas::specific_heat_ratio_;
        return std::sqrt(g * (g - 1) * e);
    }

    // Global configuration
    static auto set_gas_props(double g, double R) -> void;

  private:
    // TODO: Get rid of global static state
    static double specific_heat_ratio_;
    static double specific_gas_constant_;
};

}  // namespace jflow

#endif
