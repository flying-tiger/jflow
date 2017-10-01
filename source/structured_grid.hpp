#ifndef JFLOW_STRUCTURED_GRID_HPP
#define JFLOW_STRUCTURED_GRID_HPP

#include "jflow.hpp"
#include <array>
#include <vector>

#include <iostream>

namespace jflow {

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

  private:
    // Helper classes representing grid elements
    class iface_t {
      public:
        iface_t(const structured_grid& parent, std::size_t id)
            : parent_(parent)
            , id_(id) {}
        bool operator==(const iface_t& other) const;
        vector2 area() const;
        vector2 vertex(std::size_t n) const;

      private:
        const structured_grid& parent_;  // Parent grid containing this face
        const std::size_t id_;           // Location in the set of i-faces
    };
    class jface_t {
      public:
        jface_t(const structured_grid& parent, std::size_t id)
            : parent_(parent)
            , id_(id) {}
        bool operator==(const jface_t& other) const;
        vector2 area() const;
        vector2 vertex(std::size_t n) const;

      private:
        const structured_grid& parent_;  // Parent grid containing this face
        const std::size_t id_;           // Location in the set of i-faces
    };
    class cell_t {
      public:
        cell_t(const structured_grid& parent, std::size_t id)
            : parent_(parent)
            , id_(id) {}
        bool operator==(const cell_t& other) const;
        vector2 vertex(std::size_t i) const;
        iface_t iface(std::size_t n) const;
        jface_t jface(std::size_t n) const;
        double volume() const;

      private:
        const structured_grid& parent_;  // Parent grid containing this face
        const std::size_t id_;           // Location in the set of cells
    };

    // Helper class for iterating over grid elements
    template <typename T>
    class range {
      public:
        template <typename T>
        class iterator {
          public:
            iterator(const structured_grid& parent, std::size_t current, std::size_t stride)
                : parent_(parent)
                , current_(current)
                , stride_(stride) {}
            iterator& operator++() {
                current_ += stride_;
                return *this;
            }
            iterator& operator--() {
                current_ -= stride_;
                return *this;
            }
            bool operator==(const iterator& other) const {
                check_precondition(&parent_ == &other.parent_, "Element are from different grids.");
                return current_ == other.current_;
            }
            bool operator!=(const iterator& other) const {
                return !(*this == other);
            }
            T operator*() const {
                return (T(parent_, current_));
            }

          private:
            const structured_grid& parent_;
            std::size_t current_;
            const std::size_t stride_;
        };
        range(
            const structured_grid& parent, std::size_t start, std::size_t end,
            std::size_t stride = 1u)
            : parent_(parent)
            , start_(start)
            , end_(end)
            , stride_(stride) {}
        iterator<T> begin() {
            return iterator<T>(parent_, start_, stride_);
        }
        iterator<T> end() {
            return iterator<T>(parent_, end_, stride_);
        }

      public:
        const structured_grid& parent_;
        const std::size_t start_;
        const std::size_t end_;
        const std::size_t stride_;
    };

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
        init_face_areas();
        init_cell_volumes();
    }

    // Vertex methods
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
    const std::vector<vector2>& vertices() const {
        return vertices_;
    }

    // Cell methods
    std::size_t num_cell() const {
        return size_cell(0) * size_cell(1);
    }
    std::size_t size_cell(std::size_t dim) const {
        return size_cell_[dim];
    }
    cell_t cell(size2 coordinates) const {
        return cell_t(*this, compute_cell_id(coordinates));
    }
    cell_t cell(std::size_t i, std::size_t j) const {
        return cell({ i, j });
    }
    range<cell_t> cells() const {
        return range<cell_t>(*this, 0, num_cell());
    }

    // Constant-i face methods
    std::size_t num_iface() const {
        return size_iface(0) * size_iface(1);
    }
    std::size_t size_iface(std::size_t dim) const {
        return size_iface_[dim];
    }
    iface_t iface(size2 coordinates) const {
        return iface_t(*this, compute_iface_id(coordinates));
    }
    iface_t iface(std::size_t i, std::size_t j) const {
        return iface({ i, j });
    }
    range<iface_t> ifaces() const {
        return range<iface_t>(*this, 0, num_iface());
    }
    // range<iface_t> min_ifaces() const {
    //    return range<iface_t>(*this, 0, size_iface(1));
    //}
    // range<iface_t> max_ifaces() const {
    //    return range<iface_t>(*this, num_iface() - size_iface(1), num_iface());
    //}
    // range<iface_t> interior_ifaces() const {
    //    return range<iface_t>(*this, size_iface(1), num_iface() - size_iface(1));
    //}

    // Constant-j face methods
    std::size_t num_jface() const {
        return size_jface(0) * size_jface(1);
    }
    std::size_t size_jface(std::size_t dim) const {
        return size_jface_[dim];
    }
    jface_t jface(size2 coordinates) const {
        return jface_t(*this, compute_jface_id(coordinates));
    }
    jface_t jface(std::size_t i, std::size_t j) const {
        return jface({ i, j });
    }
    range<jface_t> jfaces() const {
        return range<jface_t>(*this, 0, num_jface());
    }
    // range<jface_t> min_jfaces() const {
    //    return range<jface_t>(*this, 0, num_jface(), size_jface(1));
    //}
    // range<jface_t> max_jfaces() const {
    //    return range<jface_t>(*this, size_jface(1) - 1, size_jface(1) + num_jface(),
    //    size_jface(1));
    //}
    // range<jface_t> interior_jfaces() const {
    // And this is where my strided approach fails
    // Option 1: Add jump_count and jump_offset to iterator so we can have regular gaps
    //   in the sequence. This will add branching and more state to the iterator. Might
    //   make sense to create two (three?) range types: linear_range, strided_range,
    //   jumping_range, etc. That way we can maximize efficiency
    // Also, consider having iterators point back to it's parent range so it doesn't have
    //   to save as much state. Extra level of indirection? Or would it be compiled out?
    //}

  private:
    void init_face_areas();
    void init_cell_volumes();

    // id => i,j calculations
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

    // i,j => id calculations
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

    size2 size_vertex_;                 // Number of vertices along each coordinate
    size2 size_cell_;                   // Number of cells along each coordinate
    size2 size_iface_;                  // Number of ifaces along each coordinate
    size2 size_jface_;                  // Number of jfaces along each coordinate
    std::vector<double> cell_volumes_;  // Volume of each grid cell
    std::vector<vector2> iface_areas_;  // Area vector for constant-i faces
    std::vector<vector2> jface_areas_;  // Area vector for constant-j faces
    std::vector<vector2> vertices_;     // Vertices defining the mesh
};

//----------------------------------------------------------------------------------
// Helper Class Implementation (Inline)
//----------------------------------------------------------------------------------
inline bool structured_grid::cell_t::operator==(const cell_t& other) const {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline vector2 structured_grid::cell_t::vertex(std::size_t n) const {
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
inline structured_grid::iface_t structured_grid::cell_t::iface(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Iface index is out of range.");
    auto [i, j] = parent_.compute_cell_coordinates(id_);
    return parent_.iface(i + n, j);
}
inline structured_grid::jface_t structured_grid::cell_t::jface(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Jface index is out of range.");
    auto [i, j] = parent_.compute_cell_coordinates(id_);
    return parent_.jface(i, j + n);
}
inline double structured_grid::cell_t::volume() const {
    return parent_.cell_volumes_[id_];
}

inline bool structured_grid::iface_t::operator==(const iface_t& other) const {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline vector2 structured_grid::iface_t::area() const {
    return parent_.iface_areas_[id_];
}
inline vector2 structured_grid::iface_t::vertex(std::size_t n) const {
    check_precondition(0 <= n && n < 2, "Vertex index is out of range.");
    auto [i, j] = parent_.compute_iface_coordinates(id_);
    return parent_.vertex(i, j + 1 - n);
}

inline bool structured_grid::jface_t::operator==(const jface_t& other) const {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline vector2 structured_grid::jface_t::area() const {
    return parent_.jface_areas_[id_];
}
inline vector2 structured_grid::jface_t::vertex(std::size_t n) const {
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