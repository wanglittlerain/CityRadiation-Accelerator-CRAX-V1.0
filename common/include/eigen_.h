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
using e_seg = std::array<ep2, 2>;
using e_rmat = Eigen::Matrix3d;

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

inline bool e_inpos(const e_seg& pab, const e_seg& pcd, double lab, e_seg& seg) {
    auto flag = [lab](const ep2& ax, const ep2& bx) {
        auto lax{ax.norm()};
        auto lbx{bx.norm()};
        auto f{0};
        if (e_dlte(lax, lab) && e_dlte(lbx, lab)) {
            f = -1;
        } else {
            if (lax > lbx) {
                f = 1;
            }
        }
        return f;
    };

    ep2 ac = pab[0] - pcd[0];
    ep2 bc = pab[1] - pcd[0];
    auto f1{flag(ac, bc)};

    ep2 ad = pab[0] - pcd[1];
    ep2 bd = pab[1] - pcd[1];
    auto f2{flag(ad, bd)};

    if (f1 != -1 && f2 != -1) {
        return false;
    }
    if (f1 == -1) {
        seg[0] = pcd[0];
    } else {
        seg[0] = pab[f1];
    } 
    if (f2 == -1) {
        seg[1] = pcd[1];
    } else {
        seg[1] = pab[f2];
    }
    return true;
}

inline int e_collinear(const e_seg& pab, const e_seg& pcd, e_seg& res) {
    ep2 ab = pab[0] - pab[1];
    ep2 ac = pab[0] - pcd[0];
    if (std::abs(e_cross(ab, ac)) < c_1e_2) {
        ep2 ad = pab[0] - pcd[1];
        if (std::abs(e_cross(ab, ad)) < c_1e_2) {
            auto lab{ab.norm()};
            ep2 cd = pcd[0] - pcd[1];
            auto lcd{cd.norm()};
            if (lab > lcd) {
                if (!e_inpos(pab, pcd, lab, res)) return -1;
            } else {
                if (!e_inpos(pcd, pab, lcd, res)) return -1;
            }
            ep2 dir = res[0] - res[1];
            auto dis{dir.norm()};
            if (dis < c_5e_2) {
                return -1;
            }
            if (std::abs(lab - dis) < c_5e_2 && std::abs(lcd - dis) < c_5e_2) {
                return 1;
            }
            return 0;
        }
    }
    return -1;
}

inline void e_push_points(epts& ps, ep3& p1, ep3& p2) {
    ps.emplace_back(p1);
    double z = p1(2);
    p1(2) = p2(2);
    ps.emplace_back(p1);
    ps.emplace_back(p2);
    p2(2) = z;
    ps.emplace_back(p2);    
}

#endif