#include "catch.hpp"
#include "jflow.hpp"

#ifdef NDEBUG
  #define REQUIRE_THROWS_IF_DEBUG  REQUIRE_NOTHROW
#else
  #define REQUIRE_THROWS_IF_DEBUG  REQUIRE_THROWS
#endif

TEST_CASE("test check_precondition function") {
    jflow::check_precondition(true, "this should not throw");
    REQUIRE_THROWS_IF_DEBUG(jflow::check_precondition(false, "this should throw"));
}

TEST_CASE("test cross product function (vector2 overload)") {
    jflow::vector2 x(3.0, 4.0);
    jflow::vector2 y(5.0, 6.0);
    REQUIRE(jflow::cross(x, y) == -2.0);
}

// TEST_CASE("test cross product function (vector3 overload)") {
//    jflow::vector3 x( 3.0, 4.0, 0.0 );
//    jflow::vector3 y( 5.0, 6.0, 0.0 );
//    REQUIRE(jflow::cross(x,y)[2] == -2.0);
//}
