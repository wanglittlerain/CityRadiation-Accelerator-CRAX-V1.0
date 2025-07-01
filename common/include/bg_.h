#pragma once
#ifndef _BOOST_GEOMETRY_H_
#define _BOOST_GEOMETRY_H_
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
using bg_point = boost::geometry::model::d2::point_xy<double>;
using bg_polygon = boost::geometry::model::polygon<bg_point, false, false>;
using bg_ring = boost::geometry::model::ring<bg_point, false, false>;
using bg_segment = boost::geometry::model::segment<bg_point>;
using bg_box = boost::geometry::model::box<bg_point>;
using Polys = std::vector<bg_polygon>;

inline bg_polygon scalePoly(const bg_polygon& poly, double scale) {
    if (boost::geometry::is_empty(poly)) return poly;
    auto res{poly};
    bg_point c(0.0, 0.0);
    boost::geometry::centroid(res, c);
    auto fun = [&](bg_point& p) {
        using namespace boost::geometry;
        detail::for_each_dimension<bg_point>([&](auto i) {
            const auto ci{get<i>(c)}; 
            set<i>(p, (get<i>(p) - ci) * scale + ci);
        });
    };
    boost::geometry::for_each_point(res, fun);
    return res;
}

inline bool pointInPolys(const bg_point& p, const Polys& polys) {
    for (auto& poly : polys) {
        if (boost::geometry::covered_by(p, poly)) return true;
        if (boost::geometry::distance(p, poly) < 1e-3) return true;
    }
    return false;
}

inline double polysArea(const Polys& polys) {
    auto area{0.0};
    for (auto& p : polys) {
        area += std::abs(boost::geometry::area(p));
    }
    return area;
}

inline bool intersectArea(const bg_polygon& poly1, const bg_polygon& poly2, double& area) {
    Polys out;
    if (boost::geometry::intersection(poly1, poly2, out) && !out.empty()) {
        area += polysArea(out);
        return true;
    }
    return false;
}

inline double intersectArea(const bg_polygon& poly, const Polys& polys) {
    auto area{0.0};
    for (auto& p : polys) {
        intersectArea(poly, p, area);
    }
    return area;
}
#endif