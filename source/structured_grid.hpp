#ifndef JFLOW_STRUCTURED_GRID_HPP
#define JFLOW_STRUCTURED_GRID_HPP

#include "jflow.hpp"
#include <vector>

namespace jflow {

//----------------------------------------------------------------------------------
// Helper Classes
//----------------------------------------------------------------------------------
// Design Note: Right now, iface and jface are two different types. This
// will be a problem when we implement iterators that loop over all faces
// in a memory-efficient order because we will want to interleave ifaces
// and jfaces. I see three options:
//
//   1. Create a base class. This is less than awesome, since types will be
//      heap allocated and methods will be resolved via v-table. We want
//      to keep these types light because a lot of them will get created
//      as we loop over the mesh. I think this is my least favorite option.
//
//   2. Merge both types into one class. This will require adding a tag so
//      we know which "type" of face we are and then branching on that tag
//      when indexing into the vertex array. (Actually, that's not true...
//      we could do something like odd indices are i-face, evens are j-face,
//      or the i-faces are consecutive then all the jfaces. Indexing
//      operations would be more complicated, but could be made branch
//      free... I think that's my favorite solution right now).
//
//   3. Have iterators/ranges return std::variants. This would screw with
//      the eventual goal of having the matrix assembly routines be mesh
//      agnostic, as they would have to branch on the type in the variant
//      before calling any methods. However, something like this may be
//      inevitable once we get tri faces in there...
//----------------------------------------------------------------------------------
class structured_grid;

class structured_grid_iface {
    friend structured_grid;

  public:
    bool operator==(const structured_grid_iface& other) const;
    const vector2& area_vector() const;
    const vector2& vertex(std::size_t i) const;

  private:
    // Only create object via structured_grid::iface()
    structured_grid_iface(const structured_grid& parent, std::size_t id)
        : parent(parent)
        , id(id) {}

    const structured_grid& parent;  // Parent grid containing this face
    const std::size_t id;           // Location in the set of i-faces
};

class structured_grid_jface {
    friend structured_grid;

  public:
    bool operator==(const structured_grid_jface& other) const;
    const vector2& area_vector() const;
    const vector2& vertex(std::size_t i) const;

  private:
    // Only create object via structured_grid::jface()
    structured_grid_jface(const structured_grid& parent, std::size_t id)
        : parent(parent)
        , id(id) {}

    const structured_grid& parent;  // Parent grid containing this face
    const std::size_t id;           // Location in the set of i-faces
};

class structured_grid_cell {
    friend structured_grid;

  public:
    bool operator==(const structured_grid_cell& other) const;
    const vector2& vertex(std::size_t i) const;
    structured_grid_iface iface(std::size_t i) const;
    structured_grid_jface jface(std::size_t i) const;
    double volume() const;

  private:
    // Only create objects via structured_grid::cell()
    structured_grid_cell(const structured_grid& parent, std::size_t id)
        : parent(parent)
        , id(id) {}

    const structured_grid& parent;  // Parent grid containing this face
    const std::size_t id;           // Location in the set of cells
};

//----------------------------------------------------------------------------------
// Main Class
//----------------------------------------------------------------------------------
// A class to represent a 2D structured grid. The figure belows shows how
// nodes/cells/faces are numbered (c => cell, if => iface, jf => jface).
// Numbering is column-major in the standard C/C++ fashion
//
//   j=2     2------jf=2------5------jf=5------8------jf=8------11
//           |                |                |                |
//           |                |                |                |
//          if=1    c=1      if=3    c=3      if=5    c=5      if=7
//           |                |                |                |
//           |                |                |                |
//   j=1     1------jf=1------4------jf=4------7------jf=7------10
//           |                |                |                |
//           |                |                |                |
//          if=0    c=0      if=2    c=2      if=4    c=4      if=6
//           |                |                |                |
//           |                |                |                |
//   j=0     0------jf=0------3------jf=3------6------jf=6------9
//
//          i=0              i=1              i=2              i=2
//
//----------------------------------------------------------------------------------
class structured_grid {
    friend structured_grid_cell;
    friend structured_grid_iface;
    friend structured_grid_jface;

    // TODO: Add iterator protocols

  public:
    structured_grid(std::size_t ni, std::size_t nj, const std::vector<vector2>& vertices)
        : ni(ni)
        , nj(nj)
        , vertices(vertices) {
        check_precondition(ni > 0, "Vertex count ni must be positive.");
        check_precondition(nj > 0, "Vertex count nj must be positive.");
        check_precondition(
            vertices.size() == ni * nj, "Length of vertex vector doesn't match ni, nj.");
        init_area_vectors();
        init_cell_volumes();
    }

    structured_grid(std::size_t ni, std::size_t nj, std::vector<vector2>&& vertices)
        : ni(ni)
        , nj(nj)
        , vertices(vertices) {
        check_precondition(ni > 0, "Vertex count ni must be positive.");
        check_precondition(nj > 0, "Vertex count nj must be positive.");
        check_precondition(
            vertices.size() == ni * nj, "Length of vertex vector doesn't match ni, nj.");
        init_area_vectors();
        init_cell_volumes();
    }

    std::size_t num_vertex() const {
        return ni * nj;
    }

    const vector2& vertex(std::size_t i) const {
        check_precondition(0 <= i && i < num_vertex(), "Vertex index is out of range.");
        return vertices[i];
    }

    const vector2& vertex(std::size_t i, std::size_t j) const {
        check_precondition(0 <= i && i < ni, "Vertex i-index is out of range.");
        check_precondition(0 <= j && j < nj, "Vertex j-index is out of range.");
        return vertices[i * nj + j];
    }

    std::size_t num_cell() const {
        return (ni - 1) * (nj - 1);
    }

    structured_grid_cell cell(std::size_t i) const {
        check_precondition(0 <= i && i < num_cell(), "Cell index is out of range.");
        return structured_grid_cell(*this, i);
    }

    structured_grid_cell cell(std::size_t i, std::size_t j) const {
        check_precondition(0 <= i && i < ni - 1, "Cell i-index is out of range.");
        check_precondition(0 <= j && j < nj - 1, "Cell j-index is out of range.");
        return structured_grid_cell(*this, i * (nj - 1) + j);
    }

    std::size_t num_iface() const {
        return ni * (nj - 1);
    }

    structured_grid_iface iface(std::size_t i) const {
        check_precondition(0 <= i && i < num_iface(), "Iface index is out of range.");
        return structured_grid_iface(*this, i);
    }

    structured_grid_iface iface(std::size_t i, std::size_t j) const {
        check_precondition(0 <= i && i < ni, "Iface i-index is out of range.");
        check_precondition(0 <= j && j < nj - 1, "Iface j-index is out of range.");
        return structured_grid_iface(*this, i * (nj - 1) + j);
    }

    std::size_t num_jface() const {
        return (ni - 1) * nj;
    }

    structured_grid_jface jface(std::size_t i) const {
        check_precondition(0 <= i && i < num_jface(), "Jface index is out of range.");
        return structured_grid_jface(*this, i);
    }

    structured_grid_jface jface(std::size_t i, std::size_t j) const {
        check_precondition(0 <= i && i < ni - 1, "Jface i-index is out of range.");
        check_precondition(0 <= j && j < nj, "Jface j-index is out of range.");
        return structured_grid_jface(*this, i * nj + j);
    }

  private:
    void init_area_vectors();
    void init_cell_volumes();

    std::size_t ni;                           // Number of vertices in i-coordinate
    std::size_t nj;                           // Number of vertices in j-coordinate
    std::vector<double> cell_volumes;         // Volume of each grid cell
    std::vector<vector2> iface_area_vectors;  // Area vector for constant-i faces
    std::vector<vector2> jface_area_vectors;  // Area vector for constant-j faces
    std::vector<vector2> vertices;            // Vertices defining the mesh
};

//----------------------------------------------------------------------------------
// Helper Class Implementation (Inline)
//----------------------------------------------------------------------------------
inline bool structured_grid_cell::operator==(const structured_grid_cell& other) const {
    return &parent == &other.parent && id == other.id;
}

inline const vector2& structured_grid_cell::vertex(std::size_t i) const {
    check_precondition(0 <= i && i < 4, "Vertex index is out of range");
    auto lower_left = id + id / (parent.nj - 1);  // Truncation intended
    switch (i) {
        case 0:
            return parent.vertices[lower_left];
        case 1:
            return parent.vertices[lower_left + parent.nj];
        case 2:
            return parent.vertices[lower_left + parent.nj + 1];
        default:
            return parent.vertices[lower_left + 1];
    }
    // TODO: test alternative implementation
    // std::array<std::size_t,4> increment  = { 0, 0, 1, 1 };
    // std::array<std::size_t,4> use_offset = { 0, 1, 1, 0 };
    // return parent.vertices[lower_left + increment[i] + use_offset[i]*parent.nj]
}

inline structured_grid_iface structured_grid_cell::iface(std::size_t i) const {
    check_precondition(0 <= i && i < 2, "Iface index is out of range.");
    return parent.iface(id + id / (parent.nj - 1) + i);  // Truncation intended
}

inline structured_grid_jface structured_grid_cell::jface(std::size_t i) const {
    check_precondition(0 <= i && i < 2, "Jface index is out of range.");
    return parent.jface(id + i * (parent.nj - 1));
}

inline double structured_grid_cell::volume() const {
    return parent.cell_volumes[id];
}

inline bool structured_grid_iface::operator==(const structured_grid_iface& other) const {
    return &parent == &other.parent && id == other.id;
}

inline const vector2& structured_grid_iface::area_vector() const {
    return parent.iface_area_vectors[id];
}

inline const vector2& structured_grid_iface::vertex(std::size_t i) const {
    check_precondition(0 <= i && i < 2, "Vertex index is out of range.");
    return parent.vertices[id + (1 - i) + id / (parent.nj - 1)];
    // Truncating integer division is intentional
}

inline const vector2& structured_grid_jface::area_vector() const {
    return parent.jface_area_vectors[id];
}

inline bool structured_grid_jface::operator==(const structured_grid_jface& other) const {
    return &parent == &other.parent && id == other.id;
}

inline const vector2& structured_grid_jface::vertex(size_t i) const {
    check_precondition(0 <= i && i < 2, "Vertex index out of range.");
    return parent.vertices[id + i * parent.nj];
}

//----------------------------------------------------------------------------------
// Grid Constuction Utilites
//----------------------------------------------------------------------------------
structured_grid make_cartesian_grid(
    vector2 xrange,  // {xmin, xmax}
    vector2 yrange,  // {ymin, ymax}
    std::size_t nx,  // Number of points in x
    std::size_t ny   // Number of points in y
);

}  // namespace jflow

#endif