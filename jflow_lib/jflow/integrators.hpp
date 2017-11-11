#ifndef JFLOW_INTEGRATORS_HPP
#define JFLOW_INTEGRATORS_HPP
#include "jflow/common.hpp"

namespace jflow {

//-----------------------------------------------------------------------------
// "system" class
//-----------------------------------------------------------------------------
// A "system" is a class that represents a system of first order differential 
// eqautions of the form:
//    
//      dx/dt = R(t,x)
//
// More specifically, a "system" is any class that:
//
//  1. defines a subclass called "state", typically a vector-like type,  which
//     completely specifies the current state of the system.
//
//  2. defines a subclass called "residual", also a vector-like type, which 
//     specifies the instanteous rate of change in the state variables.
//
//  3. defines a method compute_rhs(t,x) which computes the system residual 
//     for a given time and system state.
//
// Both the state vector and the residual must support basic arithmetic 
// operations with other instances of it's class: add, subtract, multiply 
// by scalar, and divide by scalar. 
//
// Furthermore, the state and residual classes must support addition and
// subtraction with each other.
//
// TODO: Express these requirements in code!
//-----------------------------------------------------------------------------

// Integrators
struct euler_integrator {
    template <class system>
    static void update(const system& sys, double dt, double t, typename system::state& x) {
        x += dt * sys.compute_rhs(t, x);
    }
};
struct shu_osher_integrator {
    template <class system>
    static void update(const system& sys, double dt, double t, typename system::state& x) {
        auto rhs1 = sys.compute_rhs(t, x);
        auto rhs2 = sys.compute_rhs(t + dt, x + dt * rhs1);
        x += dt * (rhs1 + rhs2) / 2.0;
    }
};
struct rk4_integrator {
    template <class system>
    static void update(const system& sys, double dt, double t, typename system::state& x) {
        auto rhs1 = sys.compute_rhs(t, x);
        auto rhs2 = sys.compute_rhs(t + 0.5 * dt, x + 0.5 * dt * rhs1);
        auto rhs3 = sys.compute_rhs(t + 0.5 * dt, x + 0.5 * dt * rhs2);
        auto rhs4 = sys.compute_rhs(t + 1.0 * dt, x + 1.0 * dt * rhs3);
        x += dt * (rhs1 + 2.0 * rhs2 + 2.0 * rhs3 + rhs4) / 6.0;
    }
};

// Integration Functions
template <typename integrator, typename system>
auto integrate(const system& sys, typename system::state x, vector2 tspan, std::size_t nsteps) {
    auto t  = tspan[0];
    auto dt = (tspan[1] - tspan[0]) / nsteps;
    while (nsteps--) {
        integrator::update(sys, dt, t, x);
        t += dt;
    }
    struct result {
        double time;
        typename system::state state;
    };
    return result{ t, std::move(x) };
};
}
#endif
