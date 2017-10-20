#include "euler.hpp"

// TODO: Try implementing a smooth 'abs' and 'max' so flux is differentiable

namespace jflow {

// Static data
double perfect_gas::specific_heat_ratio_   = 1.400;
double perfect_gas::specific_gas_constant_ = 287.058;
euler::state euler::freestream_;

// Implementation helpers
auto spectral_radius(const euler::state& q, const vector2& n) -> double {

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

// API functions
auto euler::compute_flux(const euler::state& q, const vector2& n) -> flux {

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

    // Compute flux
    auto un = u * n[0] + v * n[1];
    return { un * rho, un * rhou + p * n[0], un * rhov + p * n[1], un * (rhoE + p) };
}
auto euler::compute_freestream_flux(const state& q, const vector2& n) -> flux {
    return euler::compute_flux(freestream_, n);
}
auto euler::compute_wall_flux(const state& q, const vector2& n) -> flux {

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

    // Compute flux
    return { 0.0, p * n[0], p * n[1], 0.0 };
}
auto euler::compute_jump_flux(
    const euler::state& ql,  // Fluid state on left side of interface
    const euler::state& qr,  // Fluid state on right side of interface
    const vector2& n         // Unit vector in direction of computed flux
    ) -> euler::flux {

    // Estimate spectral radius of flux jacobian
    auto lam = std::max(spectral_radius(ql, n), spectral_radius(qr, n));

    // Compute flux w/ dissipation
    auto fl = compute_flux(ql, n);
    auto fr = compute_flux(qr, n);
    return 0.5 * (fl + fr - lam * (ql - qr));
}
auto euler::set_freestream(double p, double T, double vx, double vy) -> void {
    double rho  = perfect_gas::compute_density(p, T);
    double E    = perfect_gas::compute_energy(T) + 0.5 * (vx * vx + vy * vy);
    freestream_ = { rho, rho * vx, rho * vy, rho * E };
}
auto perfect_gas::set_gas_props(double g, double R) -> void {
    specific_heat_ratio_   = g;
    specific_gas_constant_ = R;
}

}  // namespace jflow
