#include "catch.hpp"
#include "jflow/grid_util.hpp"

TEST_CASE("test elliptic grid generator") {

    using vec2  = jflow::vector2;
    using size2 = jflow::size2;
    using jflow::constants::pi;

    // Grid parameters
    double a        = 2.0;
    auto mu_range   = vec2{ 0.0, 1.0 };
    auto nu_range   = vec2{ pi / 6, pi / 3 };
    auto exact_area = pi * a * a * sinh(2.0) / 24;

    // Create the grid
    auto grid = jflow::make_elliptic_grid(a, mu_range, nu_range, size2{ 21, 17 });

    // Verify the volume (area for 2D) is what we expect
    double volume = 0.0;
    for (auto cell : grid.cells()) {
        volume += cell.volume();
    }
    REQUIRE(volume == Approx(exact_area).margin(0.001));
}
TEST_CASE("test hyperbolic forebody grid genereator") {

    using jflow::constants::pi;
    using jflow::size2;
    const std::size_t n = 10;  // Number of cell per coordinate

    auto grid = jflow::make_hyperbolic_forebody_grid(
        2.0,                   // Length          [m]
        1.0,                   // Base radius     [m]
        0.2,                   // Nose radius     [m]
        pi / 4,                // Boundary angle  [rad]
        size2{ n + 1, n + 1 }  // Grid size
    );
    REQUIRE(grid.vertex(0, 0)[0] == Approx(0.0));
    REQUIRE(grid.vertex(0, 0)[1] == Approx(0.0));
    REQUIRE(grid.vertex(0, n)[0] == Approx(-7.136646549690036e-01));
    REQUIRE(grid.vertex(0, n)[1] == Approx(0.0));
    REQUIRE(grid.vertex(n, n)[0] == Approx(9.295030175464944e-01));
    REQUIRE(grid.vertex(n, n)[1] == Approx(2.738612787525831e+00));
}
