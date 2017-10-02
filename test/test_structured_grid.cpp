#include "catch.hpp"
#include "structured_grid.hpp"

TEST_CASE("test structured_grid class") {
    // Creates a simple 5-by-3 point test gird and verifies that our
    // grid methods work as intended. Grid looks like this:
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

    using vec2 = jflow::vector2;

    // Don't change parameters! Tests will break
    vec2 xrange = { -2.0, 2.0 };
    vec2 yrange = { -1.0, 1.0 };
    auto grid   = jflow::make_cartesian_grid(xrange, yrange, { 5, 3 });

    SECTION("Test vertex indexing via (i,j) coordinates") {
        REQUIRE((grid.vertex(2, 1) == vec2{ 0.0, 0.0 }));
        REQUIRE((grid.vertex(0, 2) == vec2{ -2.0, 1.0 }));
        REQUIRE((grid.vertex(4, 0) == vec2{ 2.0, -1.0 }));
        REQUIRE_THROWS(grid.vertex(0, -1));
        REQUIRE_THROWS(grid.vertex(0, 3));
        REQUIRE_THROWS(grid.vertex(-1, 0));
        REQUIRE_THROWS(grid.vertex(5, 0));
    }
    SECTION("Test iface indexing via (i,j) coordinates)") {
        // This *looks* backwards but is correct. We want the edge tangent
        // vector (t = v(1)-v(0)) cross the the edge normal vector to yield
        // the out-of-plane vector via righ-hand-rule. To achieve this,
        // vertex(0) must be the +j vertex.
        REQUIRE((grid.iface(0, 0).vertex(0) == grid.vertex(0, 1)));
        REQUIRE((grid.iface(0, 0).vertex(1) == grid.vertex(0, 0)));
        REQUIRE((grid.iface(4, 1).vertex(0) == grid.vertex(4, 2)));
        REQUIRE((grid.iface(4, 1).vertex(1) == grid.vertex(4, 1)));
        REQUIRE_THROWS(grid.iface(-1, 0));
        REQUIRE_THROWS(grid.iface(5, 0));
        REQUIRE_THROWS(grid.iface(0, -1));
        REQUIRE_THROWS(grid.iface(0, 2));
    }
    SECTION("Test jface indexing via (i,j) coordinates)") {
        REQUIRE((grid.jface(0, 0).vertex(0) == grid.vertex(0, 0)));
        REQUIRE((grid.jface(0, 0).vertex(1) == grid.vertex(1, 0)));
        REQUIRE((grid.jface(3, 2).vertex(0) == grid.vertex(3, 2)));
        REQUIRE((grid.jface(3, 2).vertex(1) == grid.vertex(4, 2)));
        REQUIRE_THROWS(grid.jface(-1, 0));
        REQUIRE_THROWS(grid.jface(4, 0));
        REQUIRE_THROWS(grid.jface(0, -1));
        REQUIRE_THROWS(grid.jface(0, 3));
    }
    SECTION("Test iface/jface area() accessor") {
        REQUIRE((grid.iface(0, 0).area() == vec2{ 1.0, 0.0 }));
        REQUIRE((grid.iface(4, 1).area() == vec2{ 1.0, 0.0 }));
        REQUIRE((grid.jface(0, 0).area() == vec2{ 0.0, 1.0 }));
        REQUIRE((grid.jface(3, 2).area() == vec2{ 0.0, 1.0 }));
    }
    SECTION("Test iface/jface cell() accessor") {
        REQUIRE((grid.iface(1, 1).cell(0) == grid.cell(0, 1)));
        REQUIRE((grid.iface(1, 1).cell(1) == grid.cell(1, 1)));
        REQUIRE_THROWS(grid.iface(0, 0).cell(0));
        REQUIRE_THROWS(grid.iface(4, 0).cell(1));
        REQUIRE_THROWS(grid.iface(1, 1).cell(2));

        REQUIRE((grid.jface(1, 1).cell(0) == grid.cell(1, 0)));
        REQUIRE((grid.jface(1, 1).cell(1) == grid.cell(1, 1)));
        REQUIRE_THROWS(grid.jface(0, 0).cell(0));
        REQUIRE_THROWS(grid.jface(2, 2).cell(1));
        REQUIRE_THROWS(grid.jface(1, 1).cell(2));
    }
    SECTION("Test cell indexing via (i,j) coordinates)") {
        REQUIRE((grid.cell(0, 0).vertex(0) == grid.vertex(0, 0)));
        REQUIRE((grid.cell(2, 1).vertex(2) == grid.vertex(3, 2)));
        REQUIRE_THROWS((grid.cell(-1, 0)));
        REQUIRE_THROWS((grid.cell(4, 0)));
        REQUIRE_THROWS((grid.cell(0, -1)));
        REQUIRE_THROWS((grid.cell(0, 2)));
    }
    SECTION("Test cell face indexing") {
        REQUIRE((grid.cell(1, 0).iface(0) == grid.iface(1, 0)));
        REQUIRE((grid.cell(1, 0).iface(1) == grid.iface(2, 0)));
        REQUIRE((grid.cell(2, 1).jface(0) == grid.jface(2, 1)));
        REQUIRE((grid.cell(2, 1).jface(1) == grid.jface(2, 2)));
        REQUIRE_THROWS((grid.cell(0, 0).iface(-1)));
        REQUIRE_THROWS((grid.cell(0, 0).iface(2)));
        REQUIRE_THROWS((grid.cell(0, 0).jface(-1)));
        REQUIRE_THROWS((grid.cell(0, 0).jface(2)));
    }
    SECTION("Test cell volume calculation") {
        REQUIRE((grid.cell(0, 0).volume() == 1.0));
        REQUIRE((grid.cell(3, 1).volume() == 1.0));
    }
    SECTION("Test range/iterator objects") {

        auto vbegin = grid.vertices().begin();
        auto vend   = grid.vertices().end();
        REQUIRE(grid.vertex(0, 0) == *vbegin);
        ++vbegin;
        REQUIRE(grid.vertex(0, 1) == *vbegin);
        --vend;
        REQUIRE(grid.vertex(4, 2) == *vend);
        --vend;
        REQUIRE(grid.vertex(4, 1) == *vend);

        auto cbegin = grid.cells().begin();
        auto cend   = grid.cells().end();
        REQUIRE(grid.cell(0, 0) == *cbegin);
        ++cbegin;
        REQUIRE(grid.cell(0, 1) == *cbegin);
        --cend;
        REQUIRE(grid.cell(3, 1) == *cend);
        --cend;
        REQUIRE(grid.cell(3, 0) == *cend);

        auto ibegin = grid.ifaces().begin();
        auto iend   = grid.ifaces().end();
        REQUIRE(grid.iface(0, 0) == *ibegin);
        ++ibegin;
        REQUIRE(grid.iface(0, 1) == *ibegin);
        --iend;
        REQUIRE(grid.iface(4, 1) == *iend);
        --iend;
        REQUIRE(grid.iface(4, 0) == *iend);

        auto jbegin = grid.jfaces().begin();
        auto jend   = grid.jfaces().end();
        REQUIRE(grid.jface(0, 0) == *jbegin);
        ++jbegin;
        REQUIRE(grid.jface(0, 1) == *jbegin);
        --jend;
        REQUIRE(grid.jface(3, 2) == *jend);
        --jend;
        REQUIRE(grid.jface(3, 1) == *jend);
    }
    SECTION("Test basic range-based for looping") {

        std::size_t sum;

        sum = 0;
        for (const auto& v : grid.vertices())
            ++sum;
        REQUIRE(sum == grid.num_vertex());

        sum = 0;
        for (const auto& c : grid.cells())
            ++sum;
        REQUIRE(sum == grid.num_cell());

        sum = 0;
        for (const auto& f : grid.ifaces())
            ++sum;
        REQUIRE(sum == grid.num_iface());

        sum = 0;
        for (const auto& f : grid.jfaces())
            ++sum;
        REQUIRE(sum == grid.num_jface());
    }
    SECTION("Test looping over min/max boundaries") {
        std::size_t sum, ans;

        sum = ans = 0;
        ans += grid.iface(0, 0).id();
        ans += grid.iface(0, 1).id();
        for (const auto& f : grid.min_ifaces()) {
            sum += f.id();
        }
        REQUIRE(sum == ans);

        sum = ans = 0;
        ans += grid.iface(4, 0).id();
        ans += grid.iface(4, 1).id();
        for (const auto& f : grid.max_ifaces()) {
            sum += f.id();
        }
        REQUIRE(sum == ans);

        sum = 0;
        for (const auto& f : grid.interior_ifaces()) {
            ++sum;
        }
        REQUIRE(sum == 6);

        sum = ans = 0;
        ans += grid.jface(0, 0).id();
        ans += grid.jface(1, 0).id();
        ans += grid.jface(2, 0).id();
        ans += grid.jface(3, 0).id();
        for (const auto& f : grid.min_jfaces()) {
            sum += f.id();
        }
        REQUIRE(sum == ans);

        sum = ans = 0;
        ans += grid.jface(0, 2).id();
        ans += grid.jface(1, 2).id();
        ans += grid.jface(2, 2).id();
        ans += grid.jface(3, 2).id();
        for (const auto& f : grid.max_jfaces()) {
            sum += f.id();
        }
        REQUIRE(sum == ans);

        sum = 0;
        for (const auto& f : grid.interior_jfaces()) {
            ++sum;
        }
        REQUIRE(sum == 4);
    }
}
