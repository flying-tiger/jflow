#include "jflow/finite_volume.hpp"

namespace jflow {

auto finite_volume::compute_residual(const solution& U, residual& R) -> void {

    // Shorthand
    using std::size_t;
    const auto& jump_flux = physics::compute_jump_flux;
    const auto& imin_flux = physics::compute_flux;  // Extrapolate
    const auto& imax_flux = physics::compute_flux;  // Extrapolate
    const auto& jmin_flux = physics::compute_wall_flux;
    const auto& jmax_flux = physics::compute_freestream_flux;

    // Interior flux
    for (const auto& f : grid_.interior_ifaces()) {
        auto left  = f.cell(0).id();
        auto right = f.cell(1).id();
        auto flux  = jump_flux(U[left], U[right], f.area());
        R[left] -= flux;
        R[right] += flux;
    }
    for (const auto& f : grid_.interior_jfaces()) {
        auto left  = f.cell(0).id();
        auto right = f.cell(1).id();
        auto flux  = jump_flux(U[left], U[right], f.area());
        R[left] -= flux;
        R[right] += flux;
    }

    // Boundary flux
    for (const auto& f : grid_.min_ifaces()) {
        auto id = f.cell(1).id();
        R[id] += imin_flux(U[id], f.area());
    }
    for (const auto& f : grid_.max_ifaces()) {
        auto id = f.cell(0).id();
        R[id] -= imax_flux(U[id], f.area());
    }
    for (const auto& f : grid_.min_jfaces()) {
        auto id = f.cell(1).id();
        R[id] += jmin_flux(U[id], f.area());
    }
    for (const auto& f : grid_.max_jfaces()) {
        auto id = f.cell(0).id();
        R[id] -= jmax_flux(U[id], f.area());
    }
}

}  // namespace jflow
