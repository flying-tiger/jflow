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
    auto xlim = vector2{ 0.0, 1.0 };
    auto ylim = vector2{ 0.0, 1.0 };
    auto size = size2{ 3, 4 };
    auto fv   = finite_volume(make_cartesian_grid(xlim, ylim, size));

    // Get relevant metric terms
    auto ivol = 1.0 / fv.grid().cell(0, 0).volume();
    auto area = fv.grid().cell(0, 0).jface(0).area()[1];

    SECTION("Check flow parallel to the wall had zero residual") {

        // Define fluid state
        double p = 1000.0, T = 300.0, u = 500.0, v = 0.0;
        auto test_state = euler::make_state(p, T, u, v);

        // Set the freestream and initial solution vector
        euler::set_freestream(test_state);
        auto U = fv.make_state_vector(test_state);

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
        auto U = fv.make_state_vector(interior);

        // Calculate flux differences across each cell
        auto rho   = perfect_gas::compute_density(p, T);
        auto E     = perfect_gas::compute_energy(T) + 0.5 * (u * u + v * v);
        auto H     = E + p / rho;
        auto diff0 = euler::flux{ -rho * v, 0.0, -rho * v * v, -rho * H * v };
        auto diff1 = euler::flux{ 0.0, 0.0, 0.0, 0.0 };
        auto diff2 = euler::flux{ -rho * v, 0.0, -3 * rho * v * v, -rho * v * (H + 3 * v * v) };

        // Compute and check residual
        auto time = 0.0;
        auto res  = fv.compute_rhs(time, U);
        REQUIRE(norm(res[0] - ivol * area * diff0) == Approx(0.0));
        REQUIRE(norm(res[1] - ivol * area * diff1) == Approx(0.0));
        REQUIRE(norm(res[2] - ivol * area * diff2) == Approx(0.0));
        REQUIRE(norm(res[3] - ivol * area * diff0) == Approx(0.0));
        REQUIRE(norm(res[4] - ivol * area * diff1) == Approx(0.0));
        REQUIRE(norm(res[5] - ivol * area * diff2) == Approx(0.0));
    }
}
