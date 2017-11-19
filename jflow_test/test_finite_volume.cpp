#include "catch.hpp"
#include "finite_volume.hpp"
#include "grid_util.hpp"

#include <iostream>

//----------------------------------------------------------------------
// This test cases verifies the flux integration over a 2-by-3 cartesian
// mesh. Standard boundary conditions apply: imin,imax = extrapolate,
// jmin = wall, jmax = freestream. The residual is checked for accuracy
// when the flow is parallel and perpendicular to the wall.
//----------------------------------------------------------------------
TEST_CASE("Verify flux integration") {

    using namespace jflow;

    // Create finite volume instance
    auto xrange = vector2{ -1.0, 1.0 };
    auto yrange = vector2{ -1.0, 1.0 };
    auto grid   = make_cartesian_grid(xrange, yrange, size2{ 3, 4 });
    auto fv     = finite_volume(grid);

    SECTION("Check flow parallel to the wall had zero residual") {

        // Define fluid state
        double p = 1000.0, T = 300.0, u = 500.0, v = 0.0;
        auto test_state = euler::make_state(p, T, u, v);

        // Set the freestream and initial solution vector
        euler::set_freestream(test_state);
        auto U = fv.make_state(test_state);

        // Compute residual and check results
        auto time = 0.0;
        auto res  = fv.compute_rhs(time, U);
        for (auto& cell_res : res) {
            for (auto& field_res : cell_res) {
                REQUIRE(field_res == Approx(0.0));
            }
        }
    }

    SECTION("Check flow perpendicular to the wall ") {

        // Initialize solution state vector
        double p = 1000.0, T = 300.0, u = 0.0, v = 500.0;
        auto interior   = euler::make_state(p, T, u, v);
        auto freestream = euler::make_state(p, T, u, 2 * v);

        // Set the freestream and initial solution vector
        euler::set_freestream(freestream);
        auto U = fv.make_state(interior);

        // Calculate reference fluxes
        auto rho = perfect_gas::compute_density(p, T);
        auto E   = perfect_gas::compute_energy(T) + 0.5 * (u * u + v * v);
        auto H   = E + p / rho;
        auto hi  = euler::flux{ -rho * v, 0.0, -3 * rho * v * v, -rho * v * (H + 3 * v * v) };
        auto mid = euler::flux{ 0.0, 0.0, 0.0, 0.0 };
        auto lo  = euler::flux{ -rho * v, 0.0, -rho * v * v, -rho * H * v };

        // Compute and check residual
        auto time = 0.0;
        auto res  = fv.compute_rhs(time, U);
        REQUIRE(norm(res[0] - lo) == Approx(0.0));
        REQUIRE(norm(res[1] - mid) == Approx(0.0));
        REQUIRE(norm(res[2] - hi) == Approx(0.0));
        REQUIRE(norm(res[3] - lo) == Approx(0.0));
        REQUIRE(norm(res[4] - mid) == Approx(0.0));
        REQUIRE(norm(res[5] - hi) == Approx(0.0));
    }
}
