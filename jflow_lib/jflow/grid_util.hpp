#ifndef JFLOW_GRID_UTIL_HPP
#define JFLOW_GRID_UTIL_HPP

#include "jflow/structured_grid.hpp"

namespace jflow {

auto make_cartesian_grid(
    vector2 xrange,  // {xmin, xmax}
    vector2 yrange,  // {ymin, ymax}
    size2 size       // Number of points in x,y
    ) -> structured_grid;

auto make_elliptic_grid(
    double eccentricity,  // Linear eccentricity: distance from orgin to focus [m]
    vector2 mu_range,     // {mu_min, mu_max}: coordinate along the hyperbolas
    vector2 nu_range,     // Linear eccentricity: distance from orgin to focus [m]
    size2 size            // Number of points in mu, nu
    ) -> structured_grid;

auto make_hyperbolic_forebody_grid(
    double length,          // Length of the body [m]
    double base_radius,     // Radius of the body at the base [m]
    double nose_radius,     // Radius of curvature at the nose [m]
    double boundary_angle,  // Angle of farfield boundary at outflow plane [rad]
    size2 size              // Number of point along, normal to the body
    ) -> structured_grid;

}  // namespace jflow

#endif
