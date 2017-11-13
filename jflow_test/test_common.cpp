#include "catch.hpp"
#include "jflow/common.hpp"
#include <cmath>

#ifdef NDEBUG
#define REQUIRE_THROWS_IF_DEBUG REQUIRE_NOTHROW
#else
#define REQUIRE_THROWS_IF_DEBUG REQUIRE_THROWS
#endif

TEST_CASE("test check_precondition function") {
    jflow::check_precondition(true, "this should not throw");
    REQUIRE_THROWS_IF_DEBUG(jflow::check_precondition(false, "this should throw"));
}

TEST_CASE("test 2D cross product function") {

    // Works for statically sized vectors
    jflow::vector2 x{ 3.0, 4.0 };
    jflow::vector2 y{ 5.0, 6.0 };
    REQUIRE(cross2d(x, y) == -2.0);

    // Works for dynamically sized vectors
    blaze::DynamicVector<double> u{ 1.0, 2.0 };
    blaze::DynamicVector<double> v{ 2.0, 1.0 };
    REQUIRE(cross2d(u, v) == -3.0);

    // Works for mixed types
    REQUIRE(cross2d(x, v) == -5.0);
    REQUIRE(cross2d(u, y) == -4.0);

    // Works for expressions
    REQUIRE(cross2d(x + 2 * y, u - v) == 29.0);
}

TEST_CASE("test norm function") {

    // Works for statically sized vectors
    jflow::vector2 x{ 1.0, 2.0 };
    jflow::vector3 y{ 1.0, 2.0, 3.0 };
    REQUIRE(norm(x) == std::sqrt(5.0));
    REQUIRE(norm(y) == std::sqrt(14.0));

    // Works for dynamically sized vectors
    blaze::DynamicVector<double> z{ 2.0, 3.0, 4.0 };
    REQUIRE(norm(z) == std::sqrt(29.0));

    // Work for vector expressions
    REQUIRE(norm(2 * y + z) == std::sqrt(165));
}

TEST_CASE("test blaze structured bindings support") {

    // Verify get<i> function
    blaze::StaticVector<double, 2> x{ 1.0, 2.0 };
    REQUIRE(std::get<0>(x) == 1.0);
    REQUIRE(std::get<1>(x) == 2.0);

    blaze::StaticVector<float, 3> y{ 11.0f, 12.0f, 13.0f };
    REQUIRE(std::get<2>(y) == 13.0f);

    // Verify tuple_size function
    REQUIRE(std::tuple_size<jflow::vector2>::value == 2);
    REQUIRE(std::tuple_size<jflow::vector3>::value == 3);

    // Verify structured bindings by value
    auto [x0, x1] = x;
    x1 += 1.0;
    REQUIRE(x0 == 1.0);
    REQUIRE(x1 == 3.0);
    REQUIRE(x[0] == 1.0);
    REQUIRE(x[1] == 2.0);

    // Verify structured bindings by mutable reference
    auto& [y0, y1, y2] = y;
    y0 += 1.0f;
    y1 += 2.0f;
    y2 += 3.0f;
    REQUIRE(y0 == 12.0f);
    REQUIRE(y1 == 14.0f);
    REQUIRE(y2 == 16.0f);
    REQUIRE(y[0] == 12.0f);
    REQUIRE(y[1] == 14.0f);
    REQUIRE(y[2] == 16.0f);

    // Verify structured binding via const reference
    /*
    TODO: Test this on another compiler. This is causing an error in VC++.
    Might be a bug on my side, but compiler should still give a better warning.
    blaze::StaticVector<int, 2> z{ 3, 5 };
    const auto& [z0, z1] = z;
    REQUIRE_THROWS(z0++);
    REQUIRE(z0 == 3);
    REQUIRE(z1 == 5);
    */
}
