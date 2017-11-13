#ifndef JFLOW_COMMON_HPP
#define JFLOW_COMMON_HPP

#include "blaze/Math.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace jflow {

//-----------------------------------------------------------------------
// Vocabulary Types
//-----------------------------------------------------------------------
using size2 = blaze::StaticVector<std::size_t, 2>;
using size3 = blaze::StaticVector<std::size_t, 3>;

using vector2 = blaze::StaticVector<double, 2>;
using vector3 = blaze::StaticVector<double, 3>;
using vector4 = blaze::StaticVector<double, 4>;
using vector5 = blaze::StaticVector<double, 5>;

//-----------------------------------------------------------------------
// Error Handling
//-----------------------------------------------------------------------

// Base error type from which all other library errors derive
struct runtime_error : public std::runtime_error {
    runtime_error(const std::string& what)
        : std::runtime_error(what) {}
};

struct precondition_error : public std::runtime_error {
    precondition_error(const std::string& what)
        : std::runtime_error(what) {}
};

// Use this instead of assert to throw instead of abort on error.
// This makes testing easier (test framework can't handle abort).
// Like assert, reduces to a no-op if NDEBUG is defined.
inline void check_precondition(bool check, const std::string& what) {
#ifndef NDEBUG
    if (!check) {
        throw(precondition_error(what));
    }
#endif
}

//-----------------------------------------------------------------------
// Mathematical Constants
//-----------------------------------------------------------------------
namespace constants {
    extern const double pi;
};

}  // namespace jflow

//-----------------------------------------------------------------------
// Blaze Extensions
//-----------------------------------------------------------------------
// Blaze does not provide a cross product specialization for 2D vectors,
// and does not provide vector norms. Here we add those functions to the
// blaze namespace so they can be discoved via ADL.
namespace blaze {

template <typename VT1, bool TF1, typename VT2, bool TF2>
auto cross2d(const Vector<VT1, TF1>& base_x, const Vector<VT2, TF2>& base_y) {
    const auto& x = ~base_x;
    const auto& y = ~base_y;
    return x[0] * y[1] - x[1] * y[0];
}

template <typename VT, bool TF>
auto norm(const blaze::Vector<VT, TF>& v) {
    return std::sqrt(dot(v, v));
}
}

// Blaze does not provide structured bindings for its fixed-size vectors.
// The code below adds this functionality. Based on:
//   https://blog.tartanllama.xyz/structured-bindings/
namespace std {

template <typename T, std::size_t N>
struct tuple_size<blaze::StaticVector<T, N>> : std::integral_constant<std::size_t, N> {};

template <std::size_t I, typename T, std::size_t N>
decltype(auto) get(blaze::StaticVector<T, N>& x) {
    return x[I];
}

template <std::size_t I, typename T, std::size_t N>
struct tuple_element<I, blaze::StaticVector<T, N>> {
    using type = decltype(std::get<I>(blaze::StaticVector<T, N>()));
};
}

#endif
