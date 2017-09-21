#include "euler.hpp"

namespace jflow {

    euler::fluxes euler::compute_fluxes(const euler::state& q) {

        // Create sensible names
        auto& rho  = q[euler::field::density];
        auto& rhou = q[euler::field::momentum_x];
        auto& rhov = q[euler::field::momentum_y];
        auto& rhoE = q[euler::field::total_energy];

        // Primitive state
        auto u = rhou/rho;
        auto v = rhov/rho;
        auto p = 0.4*(rhoE - 0.5*(rhou*u + rhov*v));
        auto rhoH = rhoE + p;

        // Compute fluxes
        return {
            flux{ rhou, rhou*u+p, rhov*u,   rhoH*u },
            flux{ rhov, rhou*v,   rhov*v+p, rhoH*v }
        };

    }

    // Utilizes the AUSM+ scheme of Liou 1996
    euler::flux euler::compute_interface_flux(
        const euler::state& ql,     // Fluid state on left side of interface
        const euler::state& qr,     // Fluid state on right side of interface
        const vector2& area         // Interface area vector
    ) {




        return flux{ 0.0, 0.0, 0.0, 0.0 };
    }
}
