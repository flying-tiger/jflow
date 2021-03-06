#include "jflow/structured_grid.hpp"
#include <fstream>
#include <iomanip>

namespace jflow {

// TODO: Optimize read/write routines; try to make them cache-friendly

auto structured_grid::read(const std::string& filename) -> structured_grid {
    return structured_grid::read(std::ifstream(filename));
}
auto structured_grid::read(std::istream& in) -> structured_grid {
    // Read from plot3d format, which is packed column-major.

    std::size_t nblock;
    size2 size;
    in >> nblock;
    in >> size[0] >> size[1];
    std::vector<vector2> vertices(size[0] * size[1], vector2{ 0.0, 0.0 });

    // X-coordinate array
    for (auto j = 0u; j < size[1]; ++j) {
        for (auto i = 0u; i < size[0]; ++i) {
            in >> vertices[i * size[1] + j][0];
        }
    }

    // Y-coordindate array
    for (auto j = 0u; j < size[1]; ++j) {
        for (auto i = 0u; i < size[0]; ++i) {
            in >> vertices[i * size[1] + j][1];
        }
    }

    return structured_grid(size, std::move(vertices));
}
auto structured_grid::translate(vector2 offset) -> void {
    for (auto& v : this->vertices_) {
        v += offset;
    }
    // No update_face_area: areas invariant under translation
}
auto structured_grid::update_face_areas() -> void {
    iface_areas_.clear();
    iface_areas_.reserve(size_iface());
    for (const auto& f : ifaces()) {
        auto v0 = f.vertex(0);
        auto v1 = f.vertex(1);
        iface_areas_.push_back(vector2{ -(v1[1] - v0[1]), v1[0] - v0[0] });
    }
    jface_areas_.clear();
    jface_areas_.reserve(size_jface());
    for (const auto& f : jfaces()) {
        auto v0 = f.vertex(0);
        auto v1 = f.vertex(1);
        jface_areas_.push_back(vector2{ -(v1[1] - v0[1]), v1[0] - v0[0] });
    }
}
auto structured_grid::write(const std::string& filename) const -> void {
    std::ofstream ofs(filename);
    if (!ofs) {
        throw runtime_error("Failure opening '" + filename + "': " + strerror(errno));
    }
    this->write(std::ofstream(filename));
}
auto structured_grid::write(std::ostream& out) const -> void {
    // Serialize to plot3d format, which is packed column-major.

    const std::size_t values_per_line = 4;

    out << std::setw(15) << 1 << "\n";
    out << std::setw(15) << size_vertex(0);
    out << std::setw(15) << size_vertex(1);
    out << "\n";
    out << std::scientific << std::setprecision(15);

    // X-coordinate array
    auto counter = values_per_line;
    for (auto j = 0u; j < size_vertex(1); ++j) {
        for (auto i = 0u; i < size_vertex(0); ++i) {
            out << std::setw(24) << vertex(i, j)[0];
            if (--counter == 0) {
                counter = values_per_line;
                out << "\n";
            }
        }
    }
    if (counter != values_per_line) {
        out << "\n";
    }

    // Y-coordinate array
    counter = values_per_line;
    for (auto j = 0u; j < size_vertex(1); ++j) {
        for (auto i = 0u; i < size_vertex(0); ++i) {
            out << std::setw(24) << vertex(i, j)[1];
            if (--counter == 0) {
                counter = values_per_line;
                out << "\n";
            }
        }
    }
    if (counter != values_per_line) {
        out << "\n";
    }
}

}  // namespace jflow
