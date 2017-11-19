#include "jflow/finite_volume.hpp"

namespace jflow {

auto finite_volume::compute_rhs(double t, const state& U) const -> residual {

    // Shorthand
    using std::size_t;
    const auto& jump_flux = physics::compute_jump_flux;
    const auto& imin_flux = physics::compute_flux;  // Extrapolate
    const auto& imax_flux = physics::compute_flux;  // Extrapolate
    const auto& jmin_flux = physics::compute_wall_flux;
    const auto& jmax_flux = physics::compute_freestream_flux;

    // Allocate the residual vector
    // TODO: Profile this! Should we pass in the buffer? What sort of a hit are we taking here?
    residual R = make_residual_vector();

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

    return R * inverse_volumes_;
}
auto finite_volume::make_residual_vector() const -> residual {
    return residual(grid_.size_cell());
}
auto finite_volume::make_state_vector(physics::state init) const -> state {
    return state(grid_.size_cell(), init);
}
auto finite_volume::update_inverse_volumes() -> void {
    inverse_volumes_.resize(grid_.size_cell());
    inverse_volumes_.shrinkToFit();
    auto it = inverse_volumes_.begin();
    for (auto cell : grid_.cells()) {
        *it = 1.0 / cell.volume();
        ++it;
    }
}
}  // namespace jflow
