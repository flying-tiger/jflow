#ifndef JFLOW_STRUCTURED_GRID_HPP
#define JFLOW_STRUCTURED_GRID_HPP

#include "jflow.hpp"
#include <array>
#include <vector>

#include <iostream>

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
    const vector2& vertex(std::size_t n) const;

  private:
    // Only create object via structured_grid::iface()
    structured_grid_iface(const structured_grid& parent, std::size_t id)
        : _parent(parent)
        , _id(id) {}

    const structured_grid& _parent;  // Parent grid containing this face
    const std::size_t _id;           // Location in the set of i-faces
};

class structured_grid_jface {
    friend structured_grid;

  public:
    bool operator==(const structured_grid_jface& other) const;
    const vector2& area_vector() const;
    const vector2& vertex(std::size_t n) const;

  private:
    // Only create object via structured_grid::jface()
    structured_grid_jface(const structured_grid& parent, std::size_t id)
        : _parent(parent)
        , _id(id) {}

    const structured_grid& _parent;  // Parent grid containing this face
    const std::size_t _id;           // Location in the set of i-faces
};

class structured_grid_cell {
    friend structured_grid;

  public:
    bool operator==(const structured_grid_cell& other) const;
    const vector2& vertex(std::size_t i) const;
    structured_grid_iface iface(std::size_t n) const;
    structured_grid_jface jface(std::size_t n) const;
    double volume() const;

  private:
    // Only create objects via structured_grid::cell()
    structured_grid_cell(const structured_grid& parent, std::size_t id)
        : _parent(parent)
        , _id(id) {}

    const structured_grid& _parent;  // Parent grid containing this face
    const std::size_t _id;           // Location in the set of cells
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
    structured_grid(size2 size, std::vector<vector2> vertices)
        : _size_vertex{ size }
        , _size_cell{ size[0] - 1, size[1] - 1 }
        , _size_iface{ size[0], size[1] - 1 }
        , _size_jface{ size[0] - 1, size[1] }
        , _vertices(std::move(vertices)) {
        check_precondition(size[0] >= 2, "size[0] must be 2 or more.");
        check_precondition(size[1] >= 2, "size[1] must be 2 or more.");
        check_precondition(
            _vertices.size() == size[0] * size[1],
            "Length of vertex vector doesn't match size argument.");
        init_area_vectors();
        init_cell_volumes();
    }

    std::size_t num_vertex() const {
        return _size_vertex[0] * _size_vertex[1];
    }
    std::size_t size_vertex(std::size_t dim) const {
        return _size_vertex[dim];
    }
    const vector2& vertex(std::size_t i, std::size_t j) const {
        check_precondition(0 <= i && i < _size_vertex[0], "Vertex i-index is out of range.");
        check_precondition(0 <= j && j < _size_vertex[1], "Vertex j-index is out of range.");
        return _vertices[compute_vertex_id({ i, j })];
    }
    const vector2& vertex(size2 coordinates) const {
        return vertex(coordinates[0], coordinates[1]);
    }

    std::size_t num_cell() const {
        return _size_cell[0] * _size_cell[1];
    }
    std::size_t size_cell(std::size_t dim) const {
        return _size_cell[dim];
    }
    structured_grid_cell cell(std::size_t i, std::size_t j) const {
        check_precondition(0 <= i && i < _size_cell[0], "Cell i-index is out of range.");
        check_precondition(0 <= j && j < _size_cell[1], "Cell j-index is out of range.");
        return structured_grid_cell(*this, compute_cell_id({ i, j }));
    }
    structured_grid_cell cell(size2 coordinates) const {
        return cell(coordinates[0], coordinates[1]);
    }

    std::size_t num_iface() const {
        return _size_iface[0] * _size_iface[1];
    }
    std::size_t size_iface(std::size_t dim) const {
        return _size_iface[dim];
    }
    structured_grid_iface iface(std::size_t i, std::size_t j) const {
        check_precondition(0 <= i && i < _size_iface[0], "Iface i-index is out of range.");
        check_precondition(0 <= j && j < _size_iface[1], "Iface j-index is out of range.");
        return structured_grid_iface(*this, compute_iface_id({ i, j }));
    }
    structured_grid_iface iface(size2 coordinates) const {
        return iface(coordinates[0], coordinates[1]);
    }

    std::size_t num_jface() const {
        return _size_jface[0] * _size_jface[1];
    }
    std::size_t size_jface(std::size_t dim) const {
        return _size_jface[dim];
    }
    structured_grid_jface jface(std::size_t i, std::size_t j) const {
        check_precondition(0 <= i && i < _size_jface[0], "Jface i-index is out of range.");
        check_precondition(0 <= j && j < _size_jface[1], "Jface j-index is out of range.");
        return structured_grid_jface(*this, compute_jface_id({ i, j }));
    }
    structured_grid_jface jface(size2 coordinates) const {
        return jface(coordinates[0], coordinates[1]);
    }

  private:
    void init_area_vectors();
    void init_cell_volumes();

    std::size_t compute_id(size2 coordinates, size2 size) const {
        return coordinates[0] * size[1] + coordinates[1];
    }
    std::size_t compute_vertex_id(size2 coordinates) const {
        return compute_id(coordinates, _size_vertex);
    }
    std::size_t compute_cell_id(size2 coordinates) const {
        return compute_id(coordinates, _size_cell);
    }
    std::size_t compute_iface_id(size2 coordinates) const {
        return compute_id(coordinates, _size_iface);
    }
    std::size_t compute_jface_id(size2 coordinates) const {
        return compute_id(coordinates, _size_jface);
    }

    size2 compute_coordinates(std::size_t id, size2 size) const {
        auto i = id / size[1];
        auto j = id - i * size[1];
        return { i, j };
    }
    size2 compute_vertex_coordinates(std::size_t id) const {
        return compute_coordinates(id, _size_vertex);
    }
    size2 compute_cell_coordinates(std::size_t id) const {
        return compute_coordinates(id, _size_cell);
    }
    size2 compute_iface_coordinates(std::size_t id) const {
        return compute_coordinates(id, _size_iface);
    }
    size2 compute_jface_coordinates(std::size_t id) const {
        return compute_coordinates(id, _size_jface);
    }

    size2 _size_vertex;                        // Number of vertices along each coordinate
    size2 _size_cell;                          // Number of cells along each coordinate
    size2 _size_iface;                         // Number of ifaces along each coordinate
    size2 _size_jface;                         // Number of jfaces along each coordinate
    std::vector<double> _cell_volumes;         // Volume of each grid cell
    std::vector<vector2> _iface_area_vectors;  // Area vector for constant-i faces
    std::vector<vector2> _jface_area_vectors;  // Area vector for constant-j faces
    std::vector<vector2> _vertices;            // Vertices defining the mesh
};

//----------------------------------------------------------------------------------
// Helper Class Implementation (Inline)
//----------------------------------------------------------------------------------
inline bool structured_grid_cell::operator==(const structured_grid_cell& other) const {
    return &_parent == &other._parent && _id == other._id;
}
inline const vector2& structured_grid_cell::vertex(std::size_t n) const {
    check_precondition(0 <= n && n < 4, "Vertex index is out of range");
    auto [i, j] = _parent.compute_cell_coordinates(_id);
    switch (n) {
        case 0:
            return _parent.vertex(i, j);
        case 1:
            return _parent.vertex(i + 1, j);
        case 2:
            return _parent.vertex(i + 1, j + 1);
        default:
            return _parent.vertex(i, j + 1);
    }
}
inline structured_grid_iface structured_grid_cell::iface(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Iface index is out of range.");
    auto [i, j] = _parent.compute_cell_coordinates(_id);
    return _parent.iface(i + n, j);
}
inline structured_grid_jface structured_grid_cell::jface(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Jface index is out of range.");
    auto [i, j] = _parent.compute_cell_coordinates(_id);
    return _parent.jface(i, j + n);
}
inline double structured_grid_cell::volume() const {
    return _parent._cell_volumes[_id];
}

inline bool structured_grid_iface::operator==(const structured_grid_iface& other) const {
    return &_parent == &other._parent && _id == other._id;
}
inline const vector2& structured_grid_iface::area_vector() const {
    return _parent._iface_area_vectors[_id];
}
inline const vector2& structured_grid_iface::vertex(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Vertex index is out of range.");
    auto [i, j] = _parent.compute_iface_coordinates(_id);
    return _parent.vertex(i, j + 1 - n);
}

inline bool structured_grid_jface::operator==(const structured_grid_jface& other) const {
    return &_parent == &other._parent && _id == other._id;
}
inline const vector2& structured_grid_jface::area_vector() const {
    return _parent._jface_area_vectors[_id];
}
inline const vector2& structured_grid_jface::vertex(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Vertex index out of range.");
    auto [i, j] = _parent.compute_jface_coordinates(_id);
    return _parent.vertex(i + n, j);
}

//----------------------------------------------------------------------------------
// Grid Constuction Utilites
//----------------------------------------------------------------------------------
structured_grid make_cartesian_grid(
    vector2 xrange,  // {xmin, xmax}
    vector2 yrange,  // {ymin, ymax}
    size2 size       // Number of points in x and y
);

}  // namespace jflow

#endif