#include <stdexcept>
#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace jflow {

    //-----------------------------------------------------------------------
    // Vocabulary Types
    //-----------------------------------------------------------------------
    using vector2 = Eigen::Matrix<double, 2, 1>;
    using vector3 = Eigen::Matrix<double, 3, 1>;
    using vector4 = Eigen::Matrix<double, 4, 1>;
    using vector5 = Eigen::Matrix<double, 4, 1>;

    using matrix22 = Eigen::Matrix<double, 2, 2>;
    using matrix33 = Eigen::Matrix<double, 3, 3>;
    using matrix44 = Eigen::Matrix<double, 4, 4>;
    using matrix55 = Eigen::Matrix<double, 5, 5>;


    //-----------------------------------------------------------------------
    // Error Handling
    //-----------------------------------------------------------------------

    // Base error type from which all other library errors derive
    struct runtime_error : public std::runtime_error {
        runtime_error(const std::string& what) : std::runtime_error(what) {}
    };

    struct precondition_error : public runtime_error {
        precondition_error(const std::string& what) : runtime_error(what) {}
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
    // Miscellaneous Utilities
    //-----------------------------------------------------------------------

    // Eigen does not provide a cross product specialization for 2D vectors,
    // so we provide on here. For symmetry, we define a vector3 overload as
    // well, but since it returns a vector3 it will shout-circuts Eigen's
    // expression-template engine and may not optimize as well..
    inline double cross(const vector2& x, const vector2& y) {
        return x[0]*y[1] - x[1]*y[0];
    }

    // TODO: Fix cross for vector3 inpuths.
    // This overloaad causes problems when arguments are Eigen expression
    // templates; apparently expressions templates aren't smart enough to
    // fail when implicitly converted to type w/ wrong size. Can probably
    // fix this with some SFINAE dark magic.
    //inline vector3 cross(const vector3& x, const vector3& y) {
    //    return x.cross(y); // Delegates to Eigen
    //}

}
