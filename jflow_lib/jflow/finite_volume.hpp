#ifndef JFLOW_FINITE_VOLUME_HPP
#define JFLOW_FINITE_VOLUME_HPP

#include "jflow/common.hpp"
#include "jflow/euler.hpp"
#include "jflow/structured_grid.hpp"

namespace jflow {

class finite_volume {

  public:
    using physics  = euler;
    using state    = dynamic_vector<physics::state>;
    using residual = dynamic_vector<physics::flux>;

    // Constructor
    finite_volume(structured_grid grid)
        : grid_(std::move(grid)) {
        update_inverse_volumes();
    }

    // Accessors
    auto grid() const -> const structured_grid& {
        return grid_;
    }

    // API Functions
    auto compute_rhs(double t, const state& U) const -> residual;
    auto make_residual_vector() const -> residual;
    auto make_state_vector(physics::state init = physics::state()) const -> state;

  private:
    // Initializers
    auto update_inverse_volumes() -> void;

    // Members
    structured_grid grid_;
    dynamic_vector<double> inverse_volumes_;
};

}  // namespace jflow

#endif
