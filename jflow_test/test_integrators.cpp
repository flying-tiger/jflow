#include "catch.hpp"
#include "jflow/common.hpp"
#include "jflow/integrators.hpp"
#include <cmath>
#include <utility>

#include <iostream>

//------------------------------------------------------------------------------
// Dynamic system for testing
//------------------------------------------------------------------------------
// This class implements the classical orbital mechanic equation:
//
//      ddot{r} = - mu * r / norm(r)^3
//
// Which has the closed form solution:
//
//      r(theta) = p / (1 + e*cos(theta))
//
// See https://en.wikipedia.org/wiki/Orbital_mechanics for details.

struct orbital_body {

    // Parameters
    double mu;  // Gravitational Constant: G*(m1 + m2)

    // Definitions
    using state = jflow::vector4;
    using rhs   = jflow::vector4;

    // Methods
    auto compute_rhs(double t, const state& u) const -> rhs {

        // Shorthand
        auto& px = u[0];
        auto& py = u[1];
        auto& vx = u[2];
        auto& vy = u[3];

        // Central acceleration
        auto radius = std::sqrt(px * px + py * py);
        auto accel  = mu / (radius * radius);
        auto ax     = -accel * px / radius;
        auto ay     = -accel * py / radius;

        // Compute rhs
        return rhs{ vx, vy, ax, ay };
    }
};

//------------------------------------------------------------------------------
// Helper functions
//------------------------------------------------------------------------------
template <typename func>
double convergence_rate(std::vector<std::size_t> steps_list, func calc_error) {
    double sx  = 0.0;  // These are running sums so we can compute the slope of
    double sy  = 0.0;  // the trend line without storing each (x,y) pair. See:
    double sxy = 0.0;  // https://en.wikipedia.org/wiki/Simple_linear_regression
    double sxx = 0.0;
    for (auto steps : steps_list) {
        auto y = std::log10(calc_error(steps));
        auto x = std::log10(1.0 / steps);
        sx += x;
        sy += y;
        sxy += x * y;
        sxx += x * x;
    }
    size_t n = steps_list.size();
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

//------------------------------------------------------------------------------
// Testing
//------------------------------------------------------------------------------
TEST_CASE("test integrator accuracy") {
    using namespace jflow;
    auto pi = constants::pi;

    // Define the orbit
    double rp = 1.0;  // Radius at periapsis
    double ra = 3.0;  // Radius at apoapsis
    double mu = 1.0;  // Gravitational constant

    // Compute miscellaneous orbital parameters
    auto a  = 0.5 * (rp + ra);                     // Semi major axis
    auto T  = 2 * pi * std::sqrt(a * a * a / mu);  // Orbital period
    auto vp = std::sqrt(mu * (2 / rp - 1 / a));    // Velocity at periapsis
    auto va = std::sqrt(mu * (2 / ra - 1 / a));    // Velocity at apoapsis
    auto tp = 0.0;                                 // Time at periapsis
    auto ta = 0.5 * T;                             // Time at apoapsis

    // Instantiate object and define initial/final state
    auto body  = orbital_body{ mu };
    auto begin = orbital_body::state{ -rp, 0.0, 0.0, vp };
    auto end   = orbital_body::state{ ra, 0.0, 0.0, -va };
    auto tspan = vector2{ tp, ta };

    SECTION("Propagate to apoapsis using RK4; verify 4th-order accuracy") {
        auto rate = convergence_rate({ 100, 200, 400 }, [=](auto steps) {
            auto result = integrate<rk4_integrator>(body, begin, tspan, steps);
            return norm(result.state - end);
        });
        REQUIRE(rate == Approx(4.10).margin(0.05));
    }
    SECTION("Propagate to apoapsis using Shu-Osher; verify 2nd-order accuracy") {
        auto rate = convergence_rate({ 100, 200, 400 }, [=](auto steps) {
            auto result = integrate<shu_osher_integrator>(body, begin, tspan, steps);
            return norm(result.state - end);
        });
        REQUIRE(rate == Approx(2.03).margin(0.05));
    }
    SECTION("Propagate to apoapsis using Euler forward; verify 1st-order accuracy") {
        auto rate = convergence_rate({ 200, 400, 800 }, [=](auto steps) {
            auto result = integrate<euler_integrator>(body, begin, tspan, steps);
            return norm(result.state - end);
        });
        REQUIRE(rate == Approx(0.93).margin(0.05));
    }
}
