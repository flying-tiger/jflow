#include "catch.hpp"
#include "euler.hpp"

TEST_CASE("Verify flux vector calculation") {

    using namespace jflow;

    double gamma = 1.4;

    // Primitive state
    double rho =  1.0;
    double u   =  5.0;
    double v   = -2.0;
    double p   =  1000.0;
    double E   = p/rho/(gamma-1) + 0.5*(u*u + v*v);

    // Conservative state & fluxes
    euler::state q(rho,   rho*u,     rho*v,     rho*E      );
    euler::flux fx(rho*u, rho*u*u+p, rho*u*v,   u*(rho*E+p));
    euler::flux fy(rho*v, rho*v*u,   rho*v*v+p, v*(rho*E+p));

    // Verify function output
    auto fluxes = euler::compute_fluxes(q);
    REQUIRE( (fluxes[0] - fx).norm() < 1e-12 );
    REQUIRE( (fluxes[1] - fy).norm() < 1e-12 );

}
