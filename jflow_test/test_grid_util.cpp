#include "catch.hpp"
#include "jflow/grid_util.hpp"

TEST_CASE("test elliptic grid generator") {

    using vec2 = jflow::vector2;
    using jflow::constants::pi;

    // Grid parameters
    double a        = 2.0;
    vec2 mu_range   = { 0.0, 1.0 };
    vec2 nu_range   = { pi / 6, pi / 3 };
    auto exact_area = pi * a * a * sinh(2.0) / 24;

    // Create the grid
    auto grid = jflow::make_elliptic_grid(a, mu_range, nu_range, { 21, 17 });

    // Verify the volume (area for 2D) is what we expect
    double volume = 0.0;
    for (auto cell : grid.cells()) {
        volume += cell.volume();
    }
    REQUIRE(volume == Approx(exact_area).margin(0.001));
}
