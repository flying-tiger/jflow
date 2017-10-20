#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <Eigen/Core>

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
