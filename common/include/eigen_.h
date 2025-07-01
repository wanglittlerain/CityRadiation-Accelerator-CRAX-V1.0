#pragma once
#ifndef _EIGEN_
#define _EIGEN_
#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include "const_.h"
using ep3 = Eigen::Vector3d;
using epts = std::vector<ep3>;
using ebbox = Eigen::Matrix<double, 6, 1>;
using ep2 = Eigen::Vector2d;
using e_rmat = Eigen::Matrix3d;
struct e_seg {
    std::array<ep2, 3> seg{};
    double len{};

    e_seg() = default;

    e_seg(const ep2& b, const ep2& e) {
        seg[0] = b;
        seg[1] = e;
        direction();
    }

    void direction() {
        seg[2] = seg[1] - seg[0];
        len = seg[2].norm();
    }

    constexpr ep2& operator[](int pos) noexcept {
        return seg[pos];
    }

    constexpr const ep2& operator[](int pos) const noexcept {
        return seg[pos];
    }
};

inline bool e_p_equal(const auto& p1, const auto& p2, int axis = 3, double prec = 1e-6) {
    for (int i = 0; i < axis; ++i) {
        if (std::abs(p1(i) - p2(i)) > prec) return false;
    }
    return true;
}

inline bool e_dlte(double v1, double v2, double prec = c_5e_2) {
    return v1 < v2 || std::abs(v1 - v2) <= prec;
}

inline double e_cross(const ep2& v1, const ep2& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}

inline bool e_inpos(const e_seg& sab, const e_seg& scd, e_seg& seg) {
    if (sab.len < scd.len) return e_inpos(scd, sab, seg);
    auto flag = [lab = sab.len](const ep2& ax, const ep2& bx) {
        auto lax{ax.norm()};
        auto lbx{bx.norm()};
        auto f{0};
        if (e_dlte(lax, lab) && e_dlte(lbx, lab)) {
            f = -1;
        } else if (lax > lbx) {
            f = 1;
        }
        return f;
    };

    auto f1{flag(ep2(sab[0] - scd[0]), ep2(sab[1] - scd[0]))};
    auto f2{flag(ep2(sab[0] - scd[1]), ep2(sab[1] - scd[1]))};
    if (f1 != -1 && f2 != -1) {
        return false;
    }
    f1 == -1 ? seg[0] = scd[0] : seg[0] = sab[f1];
    f2 == -1 ? seg[1] = scd[1] : seg[1] = sab[f2];
    seg.direction();
    return true;
}

inline int e_collinear(const e_seg& sab, const e_seg& scd, e_seg& res) {
    ep2 ac = sab[0] - scd[0];
    if (std::abs(e_cross(sab[2], ac)) < c_1e_2) {
        ep2 ad = sab[0] - scd[1];
        if (std::abs(e_cross(sab[2], ad)) < c_1e_2) {
            if (!e_inpos(sab, scd, res)) return -1;
            if (res.len < c_5e_2) return -1;
            return std::abs(sab.len - res.len) < c_5e_2 && std::abs(scd.len - res.len) < c_5e_2;
        }
    }
    return -1;
}

inline void e_push_points(epts& ps, ep3& p1, ep3& p2) {
    ps.reserve(4);
    ps.emplace_back(p1);
    double z = p1(2);
    p1(2) = p2(2);
    ps.emplace_back(p1);
    ps.emplace_back(p2);
    p2(2) = z;
    ps.emplace_back(p2);    
}

#endif