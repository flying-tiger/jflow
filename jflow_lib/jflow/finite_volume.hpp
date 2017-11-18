#ifndef JFLOW_FINITE_VOLUME_HPP
#define JFLOW_FINITE_VOLUME_HPP

#include "jflow/common.hpp"
#include "jflow/euler.hpp"
#include "jflow/structured_grid.hpp"

namespace jflow {

class finite_volume {

  public:
    using physics  = euler;
    using state    = blaze::DynamicVector<physics::state>;
    using residual = blaze::DynamicVector<physics::flux>;

    // TODO: This is error prone... can silently create dangling reference if grid is a
    // temporary... Could explicitly delete an r-value constructors. Should FV own it's
    // grid? If so, should this be a move-only argument to avoid duplicating the grid?
    // What about a shared_ptr? That would be the the most robust way, but then we have
    // shared mutable state, which I've so far managed to avoid.
    finite_volume(const structured_grid& grid)
        : grid_(grid) {}

    auto make_residual() const -> residual;
    auto make_state(physics::state init = physics::state()) const -> state;
    auto compute_rhs(double t, const state& U) const -> residual;

  private:
    // Should finite_volume be a stateless class like Euler? You provide the grid and the
    // state and we provide the RHS/LHS; no configuration, etc. Will have to come up with
    // a good architecture for assigning boundary conditions and other problem-dependent
    // data.
    const structured_grid& grid_;
};

}  // namespace jflow

#endif
