#define CATCH_CONFIG_MAIN
#include "blaze/Math.h"
#include "catch.hpp"

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
