import jflow.euler;
#define _CRT_NO_VA_START_VALIDATION
#include "catch.hpp"

TEST_CASE("jflow.euler module can be imported") {
    REQUIRE( jflow::euler::ndim == 2 );
    REQUIRE( jflow::euler::ndof == 4 );
}

