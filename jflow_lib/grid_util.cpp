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

}  // namespace jflow
