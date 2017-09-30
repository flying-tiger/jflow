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
    vector2 area_vector() const;
    vector2 vertex(std::size_t n) const;

  private:
    // Only create object via structured_grid::iface()
    structured_grid_iface(const structured_grid& parent, std::size_t id)
        : parent_(parent)
        , id_(id) {}

    const structured_grid& parent_;  // Parent grid containing this face
    const std::size_t id_;           // Location in the set of i-faces
};

class structured_grid_jface {
    friend structured_grid;

  public:
    bool operator==(const structured_grid_jface& other) const;
    vector2 area_vector() const;
    vector2 vertex(std::size_t n) const;

  private:
    // Only create object via structured_grid::jface()
    structured_grid_jface(const structured_grid& parent, std::size_t id)
        : parent_(parent)
        , id_(id) {}

    const structured_grid& parent_;  // Parent grid containing this face
    const std::size_t id_;           // Location in the set of i-faces
};

class structured_grid_cell {
    friend structured_grid;

  public:
    bool operator==(const structured_grid_cell& other) const;
    vector2 vertex(std::size_t i) const;
    structured_grid_iface iface(std::size_t n) const;
    structured_grid_jface jface(std::size_t n) const;
    double volume() const;

  private:
    // Only create objects via structured_grid::cell()
    structured_grid_cell(const structured_grid& parent, std::size_t id)
        : parent_(parent)
        , id_(id) {}

    const structured_grid& parent_;  // Parent grid containing this face
    const std::size_t id_;           // Location in the set of cells
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
        : size_vertex_{ size }
        , size_cell_{ size[0] - 1, size[1] - 1 }
        , size_iface_{ size[0], size[1] - 1 }
        , size_jface_{ size[0] - 1, size[1] }
        , vertices_(std::move(vertices)) {
        check_precondition(size[0] >= 2, "size[0] must be 2 or more.");
        check_precondition(size[1] >= 2, "size[1] must be 2 or more.");
        check_precondition(
            vertices_.size() == size[0] * size[1],
            "Length of vertex vector doesn't match size argument.");
        init_area_vectors();
        init_cell_volumes();
    }

    std::size_t num_vertex() const {
        return size_vertex_[0] * size_vertex_[1];
    }
    std::size_t size_vertex(std::size_t dim) const {
        return size_vertex_[dim];
    }
    vector2 vertex(size2 coordinates) const {
        return vertices_[compute_vertex_id(coordinates)];
    }
    vector2 vertex(std::size_t i, std::size_t j) const {
        return vertex({ i, j });
    }

    std::size_t num_cell() const {
        return size_cell_[0] * size_cell_[1];
    }
    std::size_t size_cell(std::size_t dim) const {
        return size_cell_[dim];
    }
    structured_grid_cell cell(size2 coordinates) const {
        return structured_grid_cell(*this, compute_cell_id(coordinates));
    }
    structured_grid_cell cell(std::size_t i, std::size_t j) const {
        return cell({ i, j });
    }

    std::size_t num_iface() const {
        return size_iface_[0] * size_iface_[1];
    }
    std::size_t size_iface(std::size_t dim) const {
        return size_iface_[dim];
    }
    structured_grid_iface iface(size2 coordinates) const {
        return structured_grid_iface(*this, compute_iface_id(coordinates));
    }
    structured_grid_iface iface(std::size_t i, std::size_t j) const {
        return iface({ i, j });
    }

    std::size_t num_jface() const {
        return size_jface_[0] * size_jface_[1];
    }
    std::size_t size_jface(std::size_t dim) const {
        return size_jface_[dim];
    }
    structured_grid_jface jface(size2 coordinates) const {
        return structured_grid_jface(*this, compute_jface_id(coordinates));
    }
    structured_grid_jface jface(std::size_t i, std::size_t j) const {
        return jface({ i, j });
    }

  private:
    void init_area_vectors();
    void init_cell_volumes();

    size2 compute_coordinates(std::size_t id, size2 size) const {
        check_precondition(0 <= id && id < size[0] * size[1], "id is out of range.");
        auto i = id / size[1];
        auto j = id - i * size[1];
        return { i, j };
    }
    size2 compute_vertex_coordinates(std::size_t id) const {
        return compute_coordinates(id, size_vertex_);
    }
    size2 compute_cell_coordinates(std::size_t id) const {
        return compute_coordinates(id, size_cell_);
    }
    size2 compute_iface_coordinates(std::size_t id) const {
        return compute_coordinates(id, size_iface_);
    }
    size2 compute_jface_coordinates(std::size_t id) const {
        return compute_coordinates(id, size_jface_);
    }

    std::size_t compute_id(size2 coords, size2 size) const {
        check_precondition(0 <= coords[0] && coords[0] < size[0], "i-index is out of range.");
        check_precondition(0 <= coords[1] && coords[1] < size[1], "j-index is out of range.");
        return coords[0] * size[1] + coords[1];
    }
    std::size_t compute_vertex_id(size2 coords) const {
        return compute_id(coords, size_vertex_);
    }
    std::size_t compute_cell_id(size2 coords) const {
        return compute_id(coords, size_cell_);
    }
    std::size_t compute_iface_id(size2 coords) const {
        return compute_id(coords, size_iface_);
    }
    std::size_t compute_jface_id(size2 coords) const {
        return compute_id(coords, size_jface_);
    }

    size2 size_vertex_;                        // Number of vertices along each coordinate
    size2 size_cell_;                          // Number of cells along each coordinate
    size2 size_iface_;                         // Number of ifaces along each coordinate
    size2 size_jface_;                         // Number of jfaces along each coordinate
    std::vector<double> cell_volumes_;         // Volume of each grid cell
    std::vector<vector2> iface_area_vectors_;  // Area vector for constant-i faces
    std::vector<vector2> jface_area_vectors_;  // Area vector for constant-j faces
    std::vector<vector2> vertices_;            // Vertices defining the mesh
};

//----------------------------------------------------------------------------------
// Helper Class Implementation (Inline)
//----------------------------------------------------------------------------------
inline bool structured_grid_cell::operator==(const structured_grid_cell& other) const {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline vector2 structured_grid_cell::vertex(std::size_t n) const {
    check_precondition(0 <= n && n < 4, "Vertex index is out of range");
    auto [i, j] = parent_.compute_cell_coordinates(id_);
    switch (n) {
        case 0:
            return parent_.vertex(i, j);
        case 1:
            return parent_.vertex(i + 1, j);
        case 2:
            return parent_.vertex(i + 1, j + 1);
        default:
            return parent_.vertex(i, j + 1);
    }
}
inline structured_grid_iface structured_grid_cell::iface(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Iface index is out of range.");
    auto [i, j] = parent_.compute_cell_coordinates(id_);
    return parent_.iface(i + n, j);
}
inline structured_grid_jface structured_grid_cell::jface(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Jface index is out of range.");
    auto [i, j] = parent_.compute_cell_coordinates(id_);
    return parent_.jface(i, j + n);
}
inline double structured_grid_cell::volume() const {
    return parent_.cell_volumes_[id_];
}

inline bool structured_grid_iface::operator==(const structured_grid_iface& other) const {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline vector2 structured_grid_iface::area_vector() const {
    return parent_.iface_area_vectors_[id_];
}
inline vector2 structured_grid_iface::vertex(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Vertex index is out of range.");
    auto [i, j] = parent_.compute_iface_coordinates(id_);
    return parent_.vertex(i, j + 1 - n);
}

inline bool structured_grid_jface::operator==(const structured_grid_jface& other) const {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline vector2 structured_grid_jface::area_vector() const {
    return parent_.jface_area_vectors_[id_];
}
inline vector2 structured_grid_jface::vertex(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Vertex index out of range.");
    auto [i, j] = parent_.compute_jface_coordinates(id_);
    return parent_.vertex(i + n, j);
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