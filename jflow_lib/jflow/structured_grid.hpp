#ifndef JFLOW_STRUCTURED_GRID_HPP
#define JFLOW_STRUCTURED_GRID_HPP

#include "jflow/common.hpp"
#include <array>
#include <istream>
#include <limits>
#include <ostream>
#include <vector>

namespace jflow {

//----------------------------------------------------------------------------------
// structured_grid
//----------------------------------------------------------------------------------
// A class to represent a 2D structured grid. The figure belows shows how
// nodes/cells/faces are numbered (c => cell, if => iface, jf => jface).
// Numbering is "row-major" in the standard C/C++ fashion (j-index fastest).
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
    class cell_handle;   // Handle object to a cell in the mesh
    class iface_handle;  // Handle object to a constant-i facet in the mesh
    class jface_handle;  // Handle object to a constant-j facet in the mesh
    template <typename T>
    class range2d;  // An object representing a range of mesh elements

  public:
    // Constructors
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
        update_dependent_members();
    }

    // Vertex methods
    auto num_vertex() const -> std::size_t {
        return size_vertex_[0] * size_vertex_[1];
    }
    auto size_vertex(std::size_t dim) const -> std::size_t {
        return size_vertex_[dim];
    }
    auto vertex(size2 coordinates) const -> vector2 {
        return vertices_[compute_vertex_id(coordinates)];
    }
    auto vertex(std::size_t i, std::size_t j) const -> vector2 {
        return vertex(size2{ i, j });
    }
    auto vertices() const -> const std::vector<vector2>& {
        // TODO: Have this return a range like other objects
        return vertices_;
    }

    // Cell methods
    auto num_cell() const -> std::size_t {
        return size_cell(0) * size_cell(1);
    }
    auto size_cell(std::size_t dim) const -> std::size_t {
        return size_cell_[dim];
    }
    auto cell(size2 coordinates) const -> cell_handle;
    auto cell(std::size_t i, std::size_t j) const -> cell_handle;
    auto cells() const -> range2d<cell_handle>;

    // Constant-i face methods
    auto num_iface() const -> std::size_t {
        return size_iface(0) * size_iface(1);
    }
    auto size_iface(std::size_t dim) const -> std::size_t {
        return size_iface_[dim];
    }
    auto iface(size2 coordinates) const -> iface_handle;
    auto iface(std::size_t i, std::size_t j) const -> iface_handle;
    auto ifaces() const -> range2d<iface_handle>;
    auto min_ifaces() const -> range2d<iface_handle>;
    auto max_ifaces() const -> range2d<iface_handle>;
    auto interior_ifaces() const -> range2d<iface_handle>;

    // Constant-j face methods
    auto num_jface() const -> std::size_t {
        return size_jface(0) * size_jface(1);
    }
    auto size_jface(std::size_t dim) const -> std::size_t {
        return size_jface_[dim];
    }
    auto jface(size2 coordinates) const -> jface_handle;
    auto jface(std::size_t i, std::size_t j) const -> jface_handle;
    auto jfaces() const -> range2d<jface_handle>;
    auto min_jfaces() const -> range2d<jface_handle>;
    auto max_jfaces() const -> range2d<jface_handle>;
    auto interior_jfaces() const -> range2d<jface_handle>;

    // Grid mutation
    auto translate(vector2 offset) -> void;

    // Serialization
    static auto read(const std::string& filename) -> structured_grid;
    static auto read(std::istream& in) -> structured_grid;
    auto write(const std::string& filename) const -> void;
    auto write(std::ostream& out) const -> void;

  private:
    // Initializers
    auto update_face_areas() -> void;
    auto update_cell_volumes() -> void;
    auto update_dependent_members() -> void {
        update_face_areas();
        update_cell_volumes();
    }

    // id => i,j calculations
    auto compute_coordinates(std::size_t id, size2 size) const -> size2;
    auto compute_vertex_coordinates(std::size_t id) const -> size2;
    auto compute_cell_coordinates(std::size_t id) const -> size2;
    auto compute_iface_coordinates(std::size_t id) const -> size2;
    auto compute_jface_coordinates(std::size_t id) const -> size2;

    // i,j => id calculations
    auto compute_id(size2 coords, size2 size) const -> std::size_t;
    auto compute_vertex_id(size2 coords) const -> std::size_t;
    auto compute_cell_id(size2 coords) const -> std::size_t;
    auto compute_iface_id(size2 coords) const -> std::size_t;
    auto compute_jface_id(size2 coords) const -> std::size_t;

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
// Helper Classes
//----------------------------------------------------------------------------------
class structured_grid::cell_handle {
  public:
    cell_handle(const structured_grid& parent, std::size_t id)
        : parent_(parent)
        , id_(id) {}
    auto operator==(const cell_handle& other) const -> bool;
    auto id() const -> std::size_t;
    auto iface(std::size_t n) const -> iface_handle;
    auto jface(std::size_t n) const -> jface_handle;
    auto vertex(std::size_t i) const -> vector2;
    auto volume() const -> double;

  private:
    const structured_grid& parent_;  // Parent grid containing this face
    const std::size_t id_;           // Location in the set of cells
};
class structured_grid::iface_handle {
  public:
    iface_handle(const structured_grid& parent, std::size_t id)
        : parent_(parent)
        , id_(id) {}
    auto operator==(const iface_handle& other) const -> bool;
    auto area() const -> vector2;
    auto cell(std::size_t n) const -> cell_handle;
    auto id() const -> std::size_t;
    auto vertex(std::size_t n) const -> vector2;

  private:
    const structured_grid& parent_;  // Parent grid containing this face
    const std::size_t id_;           // Location in the set of i-faces
};
class structured_grid::jface_handle {
  public:
    jface_handle(const structured_grid& parent, std::size_t id)
        : parent_(parent)
        , id_(id) {}
    auto operator==(const jface_handle& other) const -> bool;
    auto area() const -> vector2;
    auto cell(std::size_t n) const -> cell_handle;
    auto id() const -> std::size_t;
    auto vertex(std::size_t n) const -> vector2;

  private:
    const structured_grid& parent_;  // Parent grid containing this face
    const std::size_t id_;           // Location in the set of i-faces
};
template <typename T>
class structured_grid::range2d {
    // This range object allows iterating over a rectangular subblock of elements in
    // the grid. It essentially implements a pair of nested for loops. For example,
    //
    //   for (i = irange[0], i < jrange[1]; ++i)
    //     for (j = jrange[0], j < jrange[1]; ++j)
    //       mutate(grid.element(i,j))
    //
    // Will result in the same post-loop state as
    //
    //   for (auto& e : grid.range2d(irange, jrange))
    //     mutate(e)
    //
    // However, the order in which elements are visited may not be the same.
    // Note that in keeping with C++ convention, irange and jrange are half-open
    // intervals, e.g. [imin, imax)
  public:
    class iterator;
    range2d(const structured_grid& parent, size2 irange, size2 jrange, size2 size)
        : parent_(parent)
        , start_(parent_.compute_id(size2{ irange[0], jrange[0] }, size))
        , end_(parent_.compute_id(size2{ irange[1] - 1, jrange[0] }, size) + size[1])
        , interval_(jrange[1] - jrange[0])
        , offset_(size[1] - (jrange[1] - jrange[0])) {
        check_precondition(irange[0] < irange[1], "irange must be strictly sorted");
        check_precondition(irange[1] <= size[0], "irange exceeds grid size");
        check_precondition(jrange[0] < jrange[1], "jrange must be strictly sorted");
        check_precondition(jrange[1] <= size[1], "jrange exceeds grid size");
    }
    auto begin() -> iterator;
    auto end() -> iterator;

  private:
    const structured_grid& parent_;  // Parent grid
    const std::size_t start_;        // Index of the first element in sequence
    const std::size_t end_;          // Index of the last+next element in sequence
    const std::size_t interval_;     // Number of increments between offsets
    const std::size_t offset_;       // Offset to move from [i,jmax+1] -> [i+1,jmin]
};
template <typename T>
class structured_grid::range2d<typename T>::iterator {
  public:
    iterator(const range2d<T>& range, std::size_t current)
        : range_(range)
        , current_(current)
        , count_until_jump_(range.interval_) {}
    auto operator++() -> iterator&;
    auto operator--() -> iterator&;
    auto operator==(const iterator& other) const -> bool;
    auto operator!=(const iterator& other) const -> bool;
    auto operator*() const -> T;

  private:
    const range2d<T>& range_;       // Parent range
    std::size_t current_;           // Current position in the set of elements
    std::size_t count_until_jump_;  // Countdown until next offset
};

//----------------------------------------------------------------------------------
// structured_grid Implementation (Inline)
//----------------------------------------------------------------------------------
inline auto structured_grid::cell(size2 coordinates) const -> cell_handle {
    return cell_handle(*this, compute_cell_id(coordinates));
}
inline auto structured_grid::cell(std::size_t i, std::size_t j) const -> cell_handle {
    return cell(size2{ i, j });
}
inline auto structured_grid::iface(size2 coordinates) const -> iface_handle {
    return iface_handle(*this, compute_iface_id(coordinates));
}
inline auto structured_grid::iface(std::size_t i, std::size_t j) const -> iface_handle {
    return iface(size2{ i, j });
}
inline auto structured_grid::jface(size2 coordinates) const -> jface_handle {
    return jface_handle(*this, compute_jface_id(coordinates));
}
inline auto structured_grid::jface(std::size_t i, std::size_t j) const -> jface_handle {
    return jface(size2{ i, j });
}
inline auto structured_grid::cells() const -> range2d<cell_handle> {
    auto irange = size2{ 0, size_cell(0) };
    auto jrange = size2{ 0, size_cell(1) };
    return range2d<cell_handle>(*this, irange, jrange, size_cell_);
}
inline auto structured_grid::ifaces() const -> range2d<iface_handle> {
    auto irange = size2{ 0, size_iface(0) };
    auto jrange = size2{ 0, size_iface(1) };
    return range2d<iface_handle>(*this, irange, jrange, size_iface_);
}
inline auto structured_grid::min_ifaces() const -> range2d<iface_handle> {
    auto irange = size2{ 0, 1 };
    auto jrange = size2{ 0, size_iface(1) };
    return range2d<iface_handle>(*this, irange, jrange, size_iface_);
}
inline auto structured_grid::max_ifaces() const -> range2d<iface_handle> {
    auto irange = size2{ size_iface(0) - 1, size_iface(0) };
    auto jrange = size2{ 0, size_iface(1) };
    return range2d<iface_handle>(*this, irange, jrange, size_iface_);
}
inline auto structured_grid::interior_ifaces() const -> range2d<iface_handle> {
    auto irange = size2{ 1, size_iface(0) - 1 };
    auto jrange = size2{ 0, size_iface(1) };
    return range2d<iface_handle>(*this, irange, jrange, size_iface_);
}
inline auto structured_grid::jfaces() const -> range2d<jface_handle> {
    auto irange = size2{ 0, size_jface(0) };
    auto jrange = size2{ 0, size_jface(1) };
    return range2d<jface_handle>(*this, irange, jrange, size_jface_);
}
inline auto structured_grid::min_jfaces() const -> range2d<jface_handle> {
    auto irange = size2{ 0, size_jface(0) };
    auto jrange = size2{ 0, 1 };
    return range2d<jface_handle>(*this, irange, jrange, size_jface_);
}
inline auto structured_grid::max_jfaces() const -> range2d<jface_handle> {
    auto irange = size2{ 0, size_jface(0) };
    auto jrange = size2{ size_jface(1) - 1, size_jface(1) };
    return range2d<jface_handle>(*this, irange, jrange, size_jface_);
}
inline auto structured_grid::interior_jfaces() const -> range2d<jface_handle> {
    auto irange = size2{ 0, size_jface(0) };
    auto jrange = size2{ 1, size_jface(1) - 1 };
    return range2d<jface_handle>(*this, irange, jrange, size_jface_);
}
inline auto structured_grid::compute_coordinates(std::size_t id, size2 size) const -> size2 {
    check_precondition(0 <= id && id < size[0] * size[1], "id is out of range.");
    auto i = id / size[1];
    auto j = id - i * size[1];
    return size2{ i, j };
}
inline auto structured_grid::compute_vertex_coordinates(std::size_t id) const -> size2 {
    return compute_coordinates(id, size_vertex_);
}
inline auto structured_grid::compute_cell_coordinates(std::size_t id) const -> size2 {
    return compute_coordinates(id, size_cell_);
}
inline auto structured_grid::compute_iface_coordinates(std::size_t id) const -> size2 {
    return compute_coordinates(id, size_iface_);
}
inline auto structured_grid::compute_jface_coordinates(std::size_t id) const -> size2 {
    return compute_coordinates(id, size_jface_);
}
inline auto structured_grid::compute_id(size2 coords, size2 size) const -> std::size_t {
    check_precondition(0 <= coords[0] && coords[0] < size[0], "i-index is out of range.");
    check_precondition(0 <= coords[1] && coords[1] < size[1], "j-index is out of range.");
    return coords[0] * size[1] + coords[1];
}
inline auto structured_grid::compute_vertex_id(size2 coords) const -> std::size_t {
    return compute_id(coords, size_vertex_);
}
inline auto structured_grid::compute_cell_id(size2 coords) const -> std::size_t {
    return compute_id(coords, size_cell_);
}
inline auto structured_grid::compute_iface_id(size2 coords) const -> std::size_t {
    return compute_id(coords, size_iface_);
}
inline auto structured_grid::compute_jface_id(size2 coords) const -> std::size_t {
    return compute_id(coords, size_jface_);
}

//----------------------------------------------------------------------------------
// Helper Class Implementation (Inline)
//----------------------------------------------------------------------------------
inline auto structured_grid::cell_handle::operator==(const cell_handle& other) const -> bool {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline auto structured_grid::cell_handle::id() const -> std::size_t {
    return id_;
}
inline auto structured_grid::cell_handle::iface(std::size_t n) const -> iface_handle {
    check_precondition(0 <= n && n < 2, "Iface index is out of range.");
    auto [i, j] = parent_.compute_cell_coordinates(id_);
    return parent_.iface(i + n, j);
}
inline auto structured_grid::cell_handle::jface(std::size_t n) const -> jface_handle {
    check_precondition(0 <= n && n < 2, "Jface index is out of range.");
    auto [i, j] = parent_.compute_cell_coordinates(id_);
    return parent_.jface(i, j + n);
}
inline auto structured_grid::cell_handle::vertex(std::size_t n) const -> vector2 {
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
inline auto structured_grid::cell_handle::volume() const -> double {
    return parent_.cell_volumes_[id_];
}
inline auto structured_grid::iface_handle::operator==(const iface_handle& other) const -> bool {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline auto structured_grid::iface_handle::area() const -> vector2 {
    return parent_.iface_areas_[id_];
}
inline auto structured_grid::iface_handle::cell(std::size_t n) const -> cell_handle {
    check_precondition(n < 2, "Cell index is out of range");
    auto coords = parent_.compute_iface_coordinates(id_);
    coords[0] += (n - 1);
    return cell_handle(parent_, parent_.compute_cell_id(coords));
}
inline auto structured_grid::iface_handle::id() const -> std::size_t {
    return id_;
}
inline auto structured_grid::iface_handle::vertex(std::size_t n) const -> vector2 {
    check_precondition(0 <= n && n < 2, "Vertex index is out of range.");
    auto [i, j] = parent_.compute_iface_coordinates(id_);
    return parent_.vertex(i, j + 1 - n);
}
inline auto structured_grid::jface_handle::operator==(const jface_handle& other) const -> bool {
    return &parent_ == &other.parent_ && id_ == other.id_;
}
inline auto structured_grid::jface_handle::area() const -> vector2 {
    return parent_.jface_areas_[id_];
}
inline auto structured_grid::jface_handle::cell(std::size_t n) const -> cell_handle {
    check_precondition(n < 2, "Cell index is out of range");
    auto coords = parent_.compute_jface_coordinates(id_);
    coords[1] += (n - 1);
    return cell_handle(parent_, parent_.compute_cell_id(coords));
}
inline auto structured_grid::jface_handle::id() const -> std::size_t {
    return id_;
}
inline auto structured_grid::jface_handle::vertex(std::size_t n) const -> vector2 {
    check_precondition(0 <= n && n < 2, "Vertex index out of range.");
    auto [i, j] = parent_.compute_jface_coordinates(id_);
    return parent_.vertex(i + n, j);
}
template <typename T>
inline auto structured_grid::range2d<T>::begin() -> iterator {
    return iterator(*this, start_);
}
template <typename T>
inline auto structured_grid::range2d<T>::end() -> iterator {
    return iterator(*this, end_);
}
template <typename T>
inline auto structured_grid::range2d<T>::iterator::operator++() -> iterator& {
    ++current_;
    --count_until_jump_;
    if (count_until_jump_ == 0) {
        current_ += range_.offset_;
        count_until_jump_ = range_.interval_;
    }
    return *this;
}
template <typename T>
inline auto structured_grid::range2d<T>::iterator::operator--() -> iterator& {
    if (count_until_jump_ == range_.interval_) {
        count_until_jump_ = 0;
        current_ -= range_.offset_;
    }
    ++count_until_jump_;
    --current_;
    return *this;
}
template <typename T>
inline auto structured_grid::range2d<T>::iterator::operator==(const iterator& other) const
    -> bool {
    check_precondition(
        &range_.parent_ == &other.range_.parent_, "Element are from different grids.");
    return current_ == other.current_;
}
template <typename T>
inline auto structured_grid::range2d<T>::iterator::operator!=(const iterator& other) const
    -> bool {
    return !(*this == other);
}
template <typename T>
inline auto structured_grid::range2d<T>::iterator::operator*() const -> T {
    return (T(range_.parent_, current_));
}

//----------------------------------------------------------------------------------
// Design Notes
//----------------------------------------------------------------------------------
// [1] Currently, iface and jface are two different types. This will be a problem
// when we implement iterators that loop over all faces in a memory-efficient order
// because we will want to interleave ifaces and jfaces. I see three options:
//
//   1. Create a base class. This is less than awesome, since objects will be
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
//      inevitable once we get tri faces, etc.
//
// [2] The way we're doing ranges (proxy class that provides iterators) feels
// over-engineered. To properly loop of the interior faces of the mesh, we have to
// be able to skip over any faces that fall on the boundary. In the case of the
// j-faces, this mean that we either:
//   1. Build in the ability to insert regular gaps into the sequence
//   2. Construct some kind of lazy filter a-la the ranges library
// I've opted for approach (1) for now, but in Python this is where I would use a
// generator function. Will have to explore co-routines and the ranges TS when they
// show up in future versions of VS. Also, right now I'm using the more complicated
// range2d iterator even in places where a simple linear iterator would suffice.
// The compiler might be smart enough to figure this out and optize away the unused
// offset logic, but this might be a place where the implementation can be refined.
//----------------------------------------------------------------------------------
}  // namespace jflow

#endif
