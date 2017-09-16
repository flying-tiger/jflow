#include "structured_grid.hpp"

namespace jflow {

    double cross(const double2& x, const double2& y) {
        return x[0]*y[1] - x[1]*y[0];
    }

    void structured_grid::init_area_vectors() {
        iface_area_vectors.clear();
        iface_area_vectors.reserve(num_iface());
        for (auto n = 0u; n < num_iface(); ++n) {
            auto& v0 = this->iface(n).vertex(0);
            auto& v1 = this->iface(n).vertex(1);
            iface_area_vectors.push_back({
                -(v1[1] - v0[1]),
                v1[0] - v0[0]
            });
        }
        jface_area_vectors.clear();
        jface_area_vectors.reserve(num_jface());
        for (auto n = 0u; n < num_jface(); ++n) {
            auto& v0 = this->jface(n).vertex(0);
            auto& v1 = this->jface(n).vertex(1);
            jface_area_vectors.push_back({
                -(v1[1] - v0[1]),
                v1[0] - v0[0]
            });
        }
    }

    void structured_grid::init_cell_volumes() {
        cell_volumes.clear();
        cell_volumes.reserve(num_cell());
        for (auto n = 0u; n < num_cell(); ++n) {

            // Shorthand
            auto& v0 = this->cell(n).vertex(0);
            auto& v1 = this->cell(n).vertex(1);
            auto& v2 = this->cell(n).vertex(2);
            auto& v3 = this->cell(n).vertex(3);

            // This calculation of the volume (really the area in 2D) is exact
            // when all four vertices are co-planar. Don't use this to compute
            // the area of a general 3D quadrilateral!
            cell_volumes.push_back(0.5 * (
                cross(v1 - v0, v3 - v0) +
                cross(v3 - v2, v1 - v2)
            ));
        }
    }

    structured_grid make_cartesian_grid(
        double2 xrange,  // {xmin, xmax}
        double2 yrange,  // {ymin, ymax}
        std::size_t nx,  // Number of points in x
        std::size_t ny   // Number of points in y
    ) {
        check_precondition(nx >= 2, "nx is too small.");
        check_precondition(ny >= 2, "ny is too small.");
        std::vector<double2> vertices;
        vertices.reserve(nx*ny);
        auto dx = (xrange[1] - xrange[0]) / (nx - 1);
        auto dy = (yrange[1] - yrange[0]) / (ny - 1);
        for (auto i = 0u; i < nx; ++i) {
            for (auto j = 0u; j < ny; ++j) {
                auto x = xrange[0] + i*dx;
                auto y = yrange[0] + j*dy;
                vertices.push_back(double2{ x,y });
            }
        }
        return structured_grid(nx, ny, std::move(vertices));
    }

    // Use this function instead of assert so that we throw() instead of abort().
    // Makes things easier to test and allows error handling if kept for production.
    void check_precondition(bool check, const std::string& what) {
        #ifndef NDEBUG
        if (!check) {
            throw(precondition_error(what));
        }
        #endif
    }
}
