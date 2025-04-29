#pragma once
#ifndef _BOOST_GEOMETRY_H_
#define _BOOST_GEOMETRY_H_
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
using bg_point = boost::geometry::model::d2::point_xy<double>;
using bg_polygon = boost::geometry::model::polygon<bg_point, false, false>;
using bg_ring = boost::geometry::model::ring<bg_point, false, false>;
using bg_segment = boost::geometry::model::segment<bg_point>;
using bg_mp = boost::geometry::model::multi_point<bg_point>;
using bg_ls = boost::geometry::model::linestring<bg_point>;
using bg_mls = boost::geometry::model::multi_linestring<bg_ls>;
using bg_box = boost::geometry::model::box<bg_point>;
using Polys = std::vector<bg_polygon>;

using bg_pl = boost::geometry::model::d2::point_xy<long>;
using bg_ringl = boost::geometry::model::ring<bg_pl, false, false>;
using bg_polygonl = boost::geometry::model::polygon<bg_pl, false, false>;
using Polys_d2l = std::vector<bg_polygonl>;

static constexpr auto c_d2l{100.0};

inline bg_pl bgp_d2l(const bg_point& p) {
    bg_pl pl = {static_cast<long>(p.x() * c_d2l), static_cast<long>(p.y() * c_d2l)};
    return pl;
}

inline void poly_d2l(const bg_polygon& poly, bg_polygonl& polyl) {
	for (const auto& p : poly.outer()) {
		polyl.outer().emplace_back(bgp_d2l(p));
	}
	for (const auto& in : poly.inners()) {
		bg_ringl inl;
		for (const auto& p : in) {
			inl.emplace_back(bgp_d2l(p));
		}
		polyl.inners().emplace_back(inl);
	}
	if (!boost::geometry::is_valid(polyl)) {
		boost::geometry::correct(polyl);
	}
}

inline void polys_d2l(const Polys& polys, Polys_d2l& ploysl) {
    ploysl.clear();
    ploysl.reserve(polys.size());
    for (auto& poly : polys) {
        bg_polygonl polyl;
        poly_d2l(poly, polyl);
        ploysl.emplace_back(polyl);
    }
}

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
        auto dis{boost::geometry::distance(p, poly)};
        if (dis < 1e-3) {
            return true;
        }
    }
    return false;
}

inline bool intersectArea(const bg_polygon& poly1, const bg_polygon& poly2, double& area) {
    Polys out;
    if (boost::geometry::intersection(poly1, poly2, out) && !out.empty()) {
        for (auto& o : out) {
            area += std::abs(boost::geometry::area(o));
        }
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