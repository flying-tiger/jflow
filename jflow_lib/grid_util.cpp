#include "jflow/grid_util.hpp"
#include <cmath>

namespace jflow {

auto make_cartesian_grid(vector2 xrange, vector2 yrange, size2 size) -> structured_grid {
    check_precondition(size[0] >= 2, "nx is too small.");
    check_precondition(size[1] >= 2, "ny is too small.");
    std::vector<vector2> vertices;
    vertices.reserve(size[0] * size[1]);
    auto dx = (xrange[1] - xrange[0]) / (size[0] - 1);
    auto dy = (yrange[1] - yrange[0]) / (size[1] - 1);
    for (auto i = 0u; i < size[0]; ++i) {
        for (auto j = 0u; j < size[1]; ++j) {
            auto x = xrange[0] + i * dx;
            auto y = yrange[0] + j * dy;
            vertices.push_back(vector2{ x, y });
        }
    }
    return structured_grid(size, std::move(vertices));
}
auto make_elliptic_grid(double eccentricity, vector2 mu_range, vector2 nu_range, size2 size)
    -> structured_grid {
    using namespace std;
    const auto& a = eccentricity;  // Shorthand
    check_precondition(a >= 0, "eccentricity must be positive.");
    check_precondition(size[0] >= 2, "nx is too small.");
    check_precondition(size[1] >= 2, "ny is too small.");
    vector<vector2> vertices;
    vertices.reserve(size[0] * size[1]);
    auto dmu = (mu_range[1] - mu_range[0]) / (size[0] - 1);
    auto dnu = (nu_range[1] - nu_range[0]) / (size[1] - 1);
    for (auto i = 0u; i < size[0]; ++i) {
        for (auto j = 0u; j < size[1]; ++j) {
            auto mu = mu_range[0] + i * dmu;
            auto nu = nu_range[0] + j * dnu;
            auto x  = a * cosh(mu) * cos(nu);
            auto y  = a * sinh(mu) * sin(nu);
            vertices.push_back(vector2{ x, y });
        }
    }
    return structured_grid(size, std::move(vertices));
}
auto make_hyperbolic_forebody_grid(
    double length, double base_radius, double nose_radius, double boundary_angle, size2 size)
    -> structured_grid {
    using namespace std;
    check_precondition(length > 0.0, "Body length must be >0.");
    check_precondition(base_radius > 0.0, "Base radius must be >0.");
    check_precondition(nose_radius > 0.0, "Nose radius must be >0.");
    check_precondition(boundary_angle > 0.0, "Boundary angle must be >0.");
    check_precondition(boundary_angle < constants::pi / 2, "Boundary angle must be < pi/2.");

    // Shorthand
    auto& l     = length;
    auto& r     = base_radius;
    auto& rho   = nose_radius;
    auto& theta = boundary_angle;

    // Compute mu_max for the grid
    auto beta = r * r / (l * rho);
    check_precondition(beta >= 2.0, "Invalid parameters: R^2/(L*rho) must be 2.0 or greater.");
    auto mu_max = acosh(beta - 1);

    // Compute linear eccentricity
    auto a = l / (cosh(mu_max) - 1);
    auto b = r / sinh(mu_max);
    auto c = sqrt(a * a + b * b);

    // Compute nu range for the grid
    auto nu_min = atan(b / a);
    auto nu_max = atan(tan(theta) * tanh(mu_max));

    // Construct the grid
    auto grid = make_elliptic_grid(c, { 0.0, mu_max }, { nu_min, nu_max }, size);
    grid.translate(-1.0 * grid.vertex(0,0)); // Make nosetip (0,0)
    return grid;
}

}  // namespace jflow
