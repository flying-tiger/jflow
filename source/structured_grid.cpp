#include "structured_grid.hpp"

namespace jflow {

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