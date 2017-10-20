#ifndef JFLOW_FINITE_VOLUME_HPP
#define JFLOW_FINITE_VOLUME_HPP

#include "euler.hpp"
#include "jflow.hpp"
#include "structured_grid.hpp"

namespace jflow {

class finite_volume {

  public:
    using physics  = euler;
    using solution = std::vector<euler::state>;
    using residual = std::vector<euler::flux>;

    finite_volume(const structured_grid& grid)
        : grid_(grid) {}

    auto compute_residual(const solution& U, residual& R) -> void;

  private:
    // Should finite_volume be a stateless class like Euler? You provide the grid and the
    // state and we provide the RHS/LHS; no configuration, etc. Will have to come up with
    // a good architecture for assigning boundary conditions and other problem-dependent
    // data.
    const structured_grid& grid_;
};

}  // namespace jflow

#endif