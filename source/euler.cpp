#include "euler.hpp"

// TODO: Try implementing a smooth 'abs' and 'max' so flux is differentiable

namespace jflow {

const double perfect_gas::gamma = 1.4;

double spectral_radius(const euler::state& q, const vector2& n) {

    // Shorthand
    auto& rho  = q[euler::field::density];
    auto& rhou = q[euler::field::momentum_x];
    auto& rhov = q[euler::field::momentum_y];
    auto& rhoE = q[euler::field::total_energy];

    // Primitive state
    auto u = rhou / rho;
    auto v = rhov / rho;
    auto e = rhoE / rho - 0.5 * (u * u + v * v);
    auto c = perfect_gas::compute_sound_speed(e, rho);

    // Spectral radius
    return c + std::abs(u * n[0] + u * n[1]);
}

euler::flux euler::compute_flux(const euler::state& q, const vector2& n) {

    // Shorthand
    auto& rho  = q[field::density];
    auto& rhou = q[field::momentum_x];
    auto& rhov = q[field::momentum_y];
    auto& rhoE = q[field::total_energy];

    // Primitive state
    auto u = rhou / rho;
    auto v = rhov / rho;
    auto e = rhoE / rho - 0.5 * (u * u + v * v);
    auto p = perfect_gas::compute_pressure(e, rho);

    // Compute fluxes
    auto un = u * n[0] + v * n[1];
    return { un * rho, un * rhou + p * n[0], un * rhov + p * n[1], un * (rhoE + p) };
}

euler::flux euler::compute_jump_flux(
    const euler::state& ql,  // Fluid state on left side of interface
    const euler::state& qr,  // Fluid state on right side of interface
    const vector2& n         // Unit vector in direction of computed flux
) {

    // Estimate spectral radius of flux jacobian
    auto lam = std::max(spectral_radius(ql, n), spectral_radius(qr, n));

    // Compute flux w/ dissipation
    auto fl = compute_flux(ql, n);
    auto fr = compute_flux(qr, n);
    return 0.5 * (fl + fr - lam * (ql - qr));
}

}  // namespace jflow
