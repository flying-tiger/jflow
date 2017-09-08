#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace jflow {

    //-------------------------------------------------------------------------
    // Declarations (this should eventually be in a common module...)
    //-------------------------------------------------------------------------
    using double2 = std::array<double, 2>;

    struct runtime_error : public std::runtime_error {
        runtime_error(const std::string& what) : std::runtime_error(what) {}
    };

    struct precondition_error : public runtime_error {
        precondition_error(const std::string& what) : runtime_error(what) {}
    };

    void check_precondition(bool check, const std::string& what);


    //-------------------------------------------------------------------------
    // Helper Classes
    //-------------------------------------------------------------------------
    class structured_grid;

    class structured_grid_iface {
        friend structured_grid;

    public:
        bool operator==(const structured_grid_iface& other) const;
        //const double2& area_vector() const;
        const double2& vertex(std::size_t i) const;

    private:
        // Only create object via structured_grid::iface()
        structured_grid_iface(const structured_grid& parent, std::size_t index)
            : parent(parent), index(index) {}

        const structured_grid& parent;  // Parent grid containing this face
        const std::size_t index;        // Location in the set of i-faces

    };

    class structured_grid_jface {
        friend structured_grid;

    public:
        bool operator==(const structured_grid_jface& other) const;
        //const double2& area_vector() const;
        const double2& vertex(std::size_t i) const;

    private:
        // Only create object via structured_grid::jface()
        structured_grid_jface(const structured_grid& parent, std::size_t index)
            : parent(parent), index(index) {}

        const structured_grid& parent;  // Parent grid containing this face
        const std::size_t index;        // Location in the set of i-faces

    };

    //-------------------------------------------------------------------------
    // Main Class
    //-------------------------------------------------------------------------
    // A class to represent a 2D structured grid. The figure belows shows how
    // nodes and faces are numbered (if => iface, jf => jface).
    //
    //                  jf=2             jf=5             jf=8
    //   j=2     2----------------5----------------8----------------11
    //           |                |                |                |
    //           |                |                |                |
    //          if=1             if=3             if=5             if=7
    //           |                |                |                |
    //           |                |                |                |
    //           |      jf=1      |      jf=4      |      jf=7      |
    //   j=1     1----------------4----------------7----------------10
    //           |                |                |                |
    //           |                |                |                |
    //          if=0             if=2            if=4              if=6
    //           |                |                |                |
    //           |                |                |                |
    //           |      jf=0      |      jf=3      |      jf=6      |
    //   j=0     0----------------3----------------6----------------9
    //
    //          i=0              i=1              i=2              i=2
    //
    class structured_grid {
        friend structured_grid_iface;
        friend structured_grid_jface;

    public:
        structured_grid(std::size_t ni, std::size_t nj, const std::vector<double2>& vertices)
            : ni(ni), nj(nj), vertices(vertices) {
            check_precondition(ni > 0, "Vertex count ni must be positive.");
            check_precondition(nj > 0, "Vertex count nj must be positive.");
            check_precondition(vertices.size() == ni*nj, "Length of vertex vector doesn't match ni, nj.");
        }

        structured_grid(std::size_t ni, std::size_t nj, std::vector<double2>&& vertices)
            : ni(ni), nj(nj), vertices(vertices) {
            check_precondition(ni > 0, "Vertex count ni must be positive.");
            check_precondition(nj > 0, "Vertex count nj must be positive.");
            check_precondition(vertices.size() == ni*nj, "Length of vertex vector doesn't match ni, nj.");
        }

        const double2& vertex(std::size_t n) const {
            check_precondition(0 <= n && n < ni*nj, "Vertex index is out of range.");
            return vertices[n];
        }

        const double2& vertex(std::size_t i, std::size_t j) const {
            check_precondition(0 <= i && i < ni, "Vertex i-index is out of range.");
            check_precondition(0 <= j && j < nj, "Vertex j-index is out of range.");
            return vertices[i*nj + j];
        }

        structured_grid_iface iface(size_t n) const {
            check_precondition(0 <= n && n < ni*(nj - 1), "Iface index is out of range.");
            return structured_grid_iface(*this, n);
        }

        structured_grid_iface iface(size_t i, size_t j) const {
            check_precondition(0 <= i && i < ni, "Iface i-index is out of range.");
            check_precondition(0 <= j && j < nj - 1, "Iface j-index is out of range.");
            return structured_grid_iface(*this, i*(nj - 1) + j);
        }

        structured_grid_jface jface(size_t n) const {
            check_precondition(0 <= n && n < (ni - 1)*nj, "Jface index is out of range.");
            return structured_grid_jface(*this, n);
        }

        structured_grid_jface jface(size_t i, size_t j) const {
            check_precondition(0 <= i && i < ni - 1, "Jface i-index is out of range.");
            check_precondition(0 <= j && j < nj, "Jface j-index is out of range.");
            return structured_grid_jface(*this, i*nj + j);
        }

    private:
        std::size_t ni;                           // Number of vertices in i-coordinate
        std::size_t nj;                           // Number of vertices in j-coordinate
        std::vector<double2> vertices;            // Vertices defining the mesh
        //vector<double2> iface_area_vectors;  // Area vector for constant-i faces
        //vector<double2> jface_area_vectors;  // Area vector for constant-j faces

    };


    //-------------------------------------------------------------------------
    // Helper Class Implementation (Inline)
    //-------------------------------------------------------------------------
    inline bool structured_grid_iface::operator==(const structured_grid_iface& other) const {
        return  &parent == &other.parent &&  index == other.index;
    }

    inline const double2& structured_grid_iface::vertex(std::size_t i) const {
        check_precondition(0 <= i && i < 2, "Vertex index is out of range.");
        return parent.vertices[index + (1 - i) + index/(parent.nj - 1)];
        // Truncating integer division is intentional
    }

    inline bool structured_grid_jface::operator==(const structured_grid_jface& other) const {
        return  &parent == &other.parent &&  index == other.index;
    }

    inline const double2& structured_grid_jface::vertex(size_t i) const {
        check_precondition(0 <= i && i < 2, "Vertex index out of range.");
        return parent.vertices[index + i*parent.nj];
    }


    //-------------------------------------------------------------------------
    // Grid Constuction Utilites
    //-------------------------------------------------------------------------
    structured_grid make_cartesian_grid(
        double2 xrange, // {xmin, xmax}
        double2 yrange, // {ymin, ymax}
        std::size_t nx, // Number of points in x
        std::size_t ny  // Number of points in y
    );

}
