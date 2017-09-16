#include <iostream>
#include "catch.hpp"
#include "structured_grid.hpp"

TEST_CASE("test cross product functoin") {
    jflow::double2 x( 3.0, 4.0 );
    jflow::double2 y( 5.0, 6.0 );
    REQUIRE( jflow::cross(x,y) == -2.0 );
}

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
        REQUIRE_THROWS(grid.vertex(grid.num_vertex()));
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
        REQUIRE_THROWS(grid.iface(grid.num_iface()));
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
        REQUIRE_THROWS(grid.jface(grid.num_jface()));
    }

    SECTION("Test jface indexing via (i,j) pairs") {
        REQUIRE((grid.jface(0, 0) == grid.jface(0)));
        REQUIRE((grid.jface(2, 1) == grid.jface(7)));
        REQUIRE_THROWS(grid.jface(-1, 0));
        REQUIRE_THROWS(grid.jface(4, 0));
        REQUIRE_THROWS(grid.jface(0, -1));
        REQUIRE_THROWS(grid.jface(0, 3));
    }

    SECTION("Test iface/jface area_vector() accessor") {
        REQUIRE((grid.iface(0).area_vector() == dbl2{1.0, 0.0}));
        REQUIRE((grid.iface(9).area_vector() == dbl2{1.0, 0.0}));
        REQUIRE((grid.jface(0).area_vector() == dbl2{0.0, 1.0}));
        REQUIRE((grid.jface(11).area_vector() == dbl2{0.0, 1.0}));
    }

    SECTION("Test cell indexing (linear)") {
        REQUIRE((grid.cell(0).vertex(0) == grid.vertex(0, 0)));
        REQUIRE((grid.cell(5).vertex(2) == grid.vertex(3, 2)));
        REQUIRE_THROWS((grid.cell(-1)));
        REQUIRE_THROWS((grid.cell(grid.num_cell())));
    }

    SECTION("Test cell indexing via (i,j) pairs") {
        REQUIRE((grid.cell(1,1) == grid.cell(3)));
        REQUIRE((grid.cell(3,0) == grid.cell(6)));
        REQUIRE_THROWS((grid.cell(-1,0)));
        REQUIRE_THROWS((grid.cell(4,0)));
        REQUIRE_THROWS((grid.cell(0,-1)));
        REQUIRE_THROWS((grid.cell(0,2)));
    }

    SECTION("Test cell face indexing") {
        REQUIRE((grid.cell(2).iface(0) == grid.iface(3)));
        REQUIRE((grid.cell(2).iface(1) == grid.iface(4)));
        REQUIRE((grid.cell(5).jface(0) == grid.jface(5)));
        REQUIRE((grid.cell(5).jface(1) == grid.jface(7)));
        REQUIRE_THROWS((grid.cell(0).iface(-1)));
        REQUIRE_THROWS((grid.cell(0).iface(2)));
        REQUIRE_THROWS((grid.cell(0).jface(-1)));
        REQUIRE_THROWS((grid.cell(0).jface(2)));
    }

    SECTION("Test cell volume calculation") {
        REQUIRE((grid.cell(0).volume() == 1.0));
        REQUIRE((grid.cell(grid.num_cell() - 1).volume() == 1.0));
    }

}
