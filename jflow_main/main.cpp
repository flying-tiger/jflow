#include "jflow/common.hpp"
#include "jflow/euler.hpp"
#include "jflow/finite_volume.hpp"
#include "jflow/grid_util.hpp"
#include "jflow/integrators.hpp"
#include "yaml-cpp/yaml.h"
#include <cstdlib>
#include <iostream>

//-----------------------------------------------------------------------
// YAML Serializer for jflow::size2
//-----------------------------------------------------------------------
namespace YAML {
template <>
struct convert<jflow::size2> {
    // TODO: Generalize for any blaze vector
    using Vector = jflow::size2;

    static Node encode(const Vector& v) {
        Node node;
        for (auto element : v) {
            node.push_back(element);
        }
        return node;
    }

    static bool decode(const Node& node, Vector& v) {
        if (!node.IsSequence() || node.size() != v.size()) {
            return false;
        }
        for (std::size_t i = 0; i < node.size(); ++i) {
            v[i] = node[i].as<Vector::ElementType>();
        }
        return true;
    }
};
}

//----------------------------------------------------------------------
// Helper functions
//----------------------------------------------------------------------
YAML::Node parse_arguments(int argc, char* argv[]) {
    using namespace std;
    if (argc != 2) {
        cout << "useage: jflow <input_file>\n\n";
        cout << "Run the jflow solver using specfied input file\n\n";
        cout << "Arguments:\n";
        cout << "  input_file  Path to input file to be run\n\n";
        exit(1);
    }
    return YAML::LoadFile(argv[1]);
}

//----------------------------------------------------------------------
// Main program
//----------------------------------------------------------------------
// TODO: Report residual
// TODO: Cleanup input handling
// TODO: Add banners, etc.
// TODO:
int main(int argc, char* argv[]) {
    using namespace jflow;
    auto d2r  = constants::pi / 180;
    auto args = parse_arguments(argc, argv);

    auto freestream = euler::make_state(
        args["freestream"]["pressure"].as<double>(),
        args["freestream"]["temperature"].as<double>(),
        args["freestream"]["u_velocity"].as<double>(),
        args["freestream"]["v_velocity"].as<double>());

    auto system = finite_volume(make_hyperbolic_forebody_grid(
        args["grid"]["body_length"].as<double>(),           //
        args["grid"]["base_radius"].as<double>(),           //
        args["grid"]["nose_radius"].as<double>(),           //
        args["grid"]["boundary_angle"].as<double>() * d2r,  //
        args["grid"]["size"].as<size2>()));

    euler::set_freestream(freestream);
    auto U = system.make_state_vector(freestream);

    auto dt     = args["solver"]["timestep"].as<double>();
    auto nsteps = args["solver"]["iterations"].as<std::size_t>();
    auto tspan  = vector2{ 0.0, nsteps * dt };
    integrate<euler_integrator>(system, U, tspan, nsteps);

    return 0;
}
