#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <Eigen/Core>
#include <blaze/Math.h>

TEST_CASE("Verfiy Blaze is up and running") {
    using blaze::DynamicVector;
    using blaze::StaticVector;
    StaticVector<int, 3UL> a{ 4, -2, 5 };
    DynamicVector<int> b{ 2, 5, -3 };
    auto c = a + b;
    REQUIRE(c[0] == 6);
    REQUIRE(c[1] == 3);
    REQUIRE(c[2] == 2);
}

TEST_CASE("Verify Catch is up and running") {
    REQUIRE(1 == 1);
}

TEST_CASE("Verify Eigen is up and running") {
    using Eigen::Vector3d;
    Vector3d x(1.0, -3.0, 0.0);
    Vector3d y(4.0, 2.0, 0.0);
    Vector3d z = x + y;
    REQUIRE(z[0] == 5.0);
    REQUIRE(z[1] == -1.0);
    REQUIRE(z[2] == 0.0);
}
