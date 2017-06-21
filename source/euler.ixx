import std.core;

module euler;

namespace jflow::euler {

using State  = std::array<double,4>;
using Flux   = std::array<double,4>;
using Fluxes = std::array<Flux,2>;

Fluxes compute_fluxes(const State& q) {

    // TODO: Factor out assumption gamma = 1.4 (air)

    // Create sensible names
    auto& rho  = q[0];
    auto& rhou = q[1];
    auto& rhov = q[2];
    auto& rhoE = q[3];

    // Primitive state
    auto u = rhou/rho;
    auto v = rhov/rho;
    auto p = 0.4*(rhoE - 0.5*(rhou*u + rhov*v));
    auto rhoH = rhoE + p;

    // Compute fluxes
    return {
        Flux{ rhou, rhou*u+p, rhov*u,   rhoH*u },
        Flux{ rhov, rhou*v,   rhov*v+p, rhoH*v }
    };

}

}
