#include <array>
#include "jflow.hpp"

namespace jflow {

// TODO: Factor out hard-coded gas model (currently pg, g=1.4)
struct euler {

    // Type definitions
    using state     = vector4;
    using flux      = vector4;
    using jacobian  = matrix44;
    using fluxes    = std::array<flux,2>;
    using jacobians = std::array<jacobian,2>;

    // Field indices
    struct field {
        static const std::size_t density      = 0;
        static const std::size_t momentum_x   = 1;
        static const std::size_t momentum_y   = 2;
        static const std::size_t total_energy = 3;
    };

    // Function for debug/test
    static fluxes compute_fluxes(const state& q);
    static flux compute_interface_flux(const state& ql, const state& qr, const vector2& area);

};

// This is a rough outline of how to factor out the gas model.
//class perfect_gas {
//
//    static const gamma = 1.4;
//
//    double compute_pressure(double rho, double e) {
//        return (gamma - 1)*rho*e;
//    }
//
//    std::array<double,2> compute_pressure_derivatives(double rho, double e) {
//        return {
//            (gamma - 1)*e,
//            (gamma - 1)*rho,
//        };
//    }
//
//}


}
