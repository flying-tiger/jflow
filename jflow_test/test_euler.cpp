#include "catch.hpp"
#include "jflow/euler.hpp"

TEST_CASE("Verify Euler flux vector calculation") {

    using namespace jflow;

    const double gamma = 1.4;
    const double Rgas  = 287.0;
    const double tol   = 1e-12;

    // Primitive state
    const double rho = 1.0;
    const double u   = 5.0;
    const double v   = -2.0;
    const double p   = 1000.0;
    const double E   = p / rho / (gamma - 1) + 0.5 * (u * u + v * v);

    // Conservative state & fluxes
    const euler::state q(rho, rho * u, rho * v, rho * E);
    const euler::flux fx(rho * u, rho * u * u + p, rho * u * v, u * (rho * E + p));
    const euler::flux fy(rho * v, rho * v * u, rho * v * v + p, v * (rho * E + p));

    // Configure perfect gas model
    perfect_gas::set_gas_props(gamma, Rgas);

    SECTION("Verify basic flux calculation") {
        auto fx_calc = euler::compute_flux(q, { 1, 0 });
        auto fy_calc = euler::compute_flux(q, { 0, 1 });
        REQUIRE((fx_calc - fx).norm() < 1e-12);
        REQUIRE((fy_calc - fy).norm() < 1e-12);
    }

    SECTION("Verify jump flux calculation") {
        // F_jump(q,q) should return F(q)
        auto fx_calc = euler::compute_jump_flux(q, q, { 1, 0 });
        auto fy_calc = euler::compute_jump_flux(q, q, { 0, 1 });
        REQUIRE((fx_calc - fx).norm() < 1e-12);
        REQUIRE((fy_calc - fy).norm() < 1e-12);
    }
}
