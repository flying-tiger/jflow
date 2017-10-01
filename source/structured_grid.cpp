#include "structured_grid.hpp"

// TODO: Replace indexed loops over iface, jface, cell with iterator traversal

namespace jflow {

auto structured_grid::init_face_areas() -> void {
    iface_areas_.clear();
    iface_areas_.reserve(num_iface());
    for (auto n = 0u; n < num_iface(); ++n) {
        auto f  = iface(compute_iface_coordinates(n));
        auto v0 = f.vertex(0);
        auto v1 = f.vertex(1);
        iface_areas_.push_back({ -(v1[1] - v0[1]), v1[0] - v0[0] });
    }
    jface_areas_.clear();
    jface_areas_.reserve(num_jface());
    for (auto n = 0u; n < num_jface(); ++n) {
        auto f  = jface(compute_jface_coordinates(n));
        auto v0 = f.vertex(0);
        auto v1 = f.vertex(1);
        jface_areas_.push_back({ -(v1[1] - v0[1]), v1[0] - v0[0] });
    }
}

auto structured_grid::init_cell_volumes() -> void {
    cell_volumes_.clear();
    cell_volumes_.reserve(num_cell());
    for (auto n = 0u; n < num_cell(); ++n) {
        auto c  = cell(compute_cell_coordinates(n));
        auto v0 = c.vertex(0);
        auto v1 = c.vertex(1);
        auto v2 = c.vertex(2);
        auto v3 = c.vertex(3);

        // This calculation of the volume (really the area in 2D) is exact
        // when all four vertices are co-planar. Don't use this to compute
        // the area of a general 3D quadrilateral!
        cell_volumes_.push_back(0.5 * (cross(v1 - v0, v3 - v0) + cross(v3 - v2, v1 - v2)));
    }
}

auto make_cartesian_grid(vector2 xrange,  vector2 yrange, size2 size) -> structured_grid {
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

}  // namespace jflow
