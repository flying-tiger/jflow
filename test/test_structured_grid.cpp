#include <iostream>
#include "catch.hpp"
#include "structured_grid.hpp"

TEST_CASE("test structured_grid class") {
    // Creates a simple 5-by-3 point test gird and verifies that our
    // indexing operations work as intended. Grid looks like this:
    //
    //        x=-2.0                  x=2.0
    //
    //  j=2     2-----5-----8-----11----14    y=1.0
    //          |     |     |     |     |
    //          |     |     |     |     |
    //  j=1     1-----4-----7-----10----13
    //          |     |     |     |     |
    //          |     |     |     |     |
    //  j=0     0-----3-----6-----9-----12    y=-1.0
    //
    //         i=0   i=1   i=2   i=3   i=4

    using dbl2 = jflow::double2;

    // Don't change parameters! Tests will break
    dbl2 xrange = { -2.0, 2.0 };
    dbl2 yrange = { -1.0, 1.0 };
    std::size_t nx = 5, ny = 3;
    auto grid = jflow::make_cartesian_grid(xrange, yrange, nx, ny);

    SECTION("Test vertex indexing (linear)") {
        REQUIRE((grid.vertex(0) == dbl2{ -2.0, -1.0 }));
        REQUIRE((grid.vertex(14) == dbl2{ 2.0, 1.0 }));
        REQUIRE_THROWS(grid.vertex(-1));
        REQUIRE_THROWS(grid.vertex(15));
    }

    SECTION("Test vertex indexing via (i,j) pairs") {        
        REQUIRE((grid.vertex(2,1) == dbl2{ 0.0, 0.0 }));
        REQUIRE((grid.vertex(0,2) == dbl2{ -2.0, 1.0 }));
        REQUIRE((grid.vertex(4,0) == dbl2{ 2.0, -1.0 }));
        REQUIRE_THROWS(grid.vertex(0, -1));
        REQUIRE_THROWS(grid.vertex(0, 3));
        REQUIRE_THROWS(grid.vertex(-1, 0)); 
        REQUIRE_THROWS(grid.vertex(5, 0));
    }

    SECTION("Test iface indexing (linear)") {
        // This *looks* backwards but is correct. We want the edge tangent 
        // vector (t = v(1)-v(0)) cross the the edge normal vector to yield
        // the out-of-plane vector via righ-hand-rule. To achieve this, 
        // vertex(0) must be the +j vertex.
        REQUIRE((grid.iface(0).vertex(0) == grid.vertex(0,1)));
        REQUIRE((grid.iface(0).vertex(1) == grid.vertex(0,0)));
        REQUIRE((grid.iface(9).vertex(0) == grid.vertex(4,2)));
        REQUIRE((grid.iface(9).vertex(1) == grid.vertex(4,1)));
        REQUIRE_THROWS(grid.iface(-1));
        REQUIRE_THROWS(grid.iface(10));
    }

    SECTION("Test iface indexing via (i,j) pairs") {
        REQUIRE((grid.iface(0, 0) == grid.iface(0)));
        REQUIRE((grid.iface(2, 1) == grid.iface(5)));
        REQUIRE_THROWS(grid.iface(-1, 0));
        REQUIRE_THROWS(grid.iface(5, 0));
        REQUIRE_THROWS(grid.iface(0, -1));
        REQUIRE_THROWS(grid.iface(0, 2));
    }

    SECTION("Test jface indexing (linear)") {        
        REQUIRE((grid.jface(0).vertex(0) == grid.vertex(0, 0)));
        REQUIRE((grid.jface(0).vertex(1) == grid.vertex(1, 0)));
        REQUIRE((grid.jface(11).vertex(0) == grid.vertex(3, 2)));
        REQUIRE((grid.jface(11).vertex(1) == grid.vertex(4, 2)));
        REQUIRE_THROWS(grid.jface(-1));
        REQUIRE_THROWS(grid.jface(12));
    }

    SECTION("Test jface indexing via (i,j) pairs") {
        REQUIRE((grid.jface(0, 0) == grid.jface(0)));
        REQUIRE((grid.jface(2, 1) == grid.jface(7)));
        REQUIRE_THROWS(grid.jface(-1, 0));
        REQUIRE_THROWS(grid.jface(4, 0));
        REQUIRE_THROWS(grid.jface(0, -1));
        REQUIRE_THROWS(grid.jface(0, 3));
    }

}