#pragma once
#ifndef _Generate_Polyhedron_
#define _Generate_Polyhedron_
#include "base.h"
#include "concave.h"
struct InstUnit {
    double minH{};
    e_seg seg;
};
struct SegmentInsts {
    bool integral{};
    double minH{};
    std::vector<InstUnit> units;
};
using Insts = std::map<int, SegmentInsts>;

struct InnerSeg {
    double minh{};
    ep2 norm{};
    e_seg seg;
    inline void swap() {
        if (seg[2].dot(norm) < 0.0) {
            std::swap(seg[0], seg[1]);
        }
    }
    inline bool operator < (const InnerSeg& rhs) const {
        ep2 dir{seg[0] - rhs.seg[0]};
        return dir.dot(norm) < 0.0;
    }
};
using Isegs = std::vector<InnerSeg>;

class GeneratePolyhedron {
public:
    bool stretch(Polyhedron& phd, const auto& points, const auto& poly, 
        const Insts& insts, double maxz, double minz) const {
        auto& out{poly.outer()};
        if (out.empty() || std::abs(maxz - minz) < c_5e_2 || minz > maxz) {
            return false;
        }
        auto size{static_cast<int>(points.size())};
        auto collinear = [&](const auto& p, int& self, int other) {
            for (int i = 0; i < size; ++i) {
                auto j{(i + 1) % size};
                if (std::abs(cxd::handedness(points[i], p, points[j])) < c_1e_2) { //3p1l
                    if (other == i) {
                        self = j;
                    } else if (other == j) {
                        self = i;
                    }
                    return self == i || self == j;
                }
            }
            return false;
        };

        auto realId = [&](const auto& pre, const auto& cur) {
            auto id{c_virtual_srf};
            if (points.empty()) return id;
            auto pid{-1};
            auto cid{-1};
            for (int i = 0; i < size; ++i) {
                auto& pt{points[i]};
                if (pid == -1 && e_p_equal(pt, pre, 2, c_5e_2)) {
                    pid = i;
                }
                if (cid == -1 && e_p_equal(pt, cur, 2, c_5e_2)) {
                    cid = i;
                }
                if (pid >= 0 && cid >= 0) break;
            }
            if (pid == -1 && cid == -1) return id;
            if (pid == -1 && !collinear(pre, pid, cid)) return id;
            if (cid == -1 && !collinear(cur, cid, pid)) return id;
            auto dis{cid - pid};
            auto mid{cid};
            if (dis < 0) {
                dis *= -1;
                mid = pid;
            }
            if (dis == 1) {
                id = mid + 1;
            } else if (dis == size - 1) {
                id = 1;
            }
            return id;
        };

        auto& srfs = phd.srfs;
        auto add_srf = [&](const auto& pb, const auto& pe, double minh, bool gap = false) {
            Surface srf;
            srf.id = static_cast<int>(srfs.size());
            srf._h = maxz;
            srf._gap = gap;
            ep3 p1(pb.x(), pb.y(), maxz);
            ep3 p2(pe.x(), pe.y(), minh);
            e_push_points(srf._pts, p1, p2);
            srf.gseg = bg_segment(bg_point{pb.x(), pb.y()}, bg_point{pe.x(), pe.y()});
            srfs.emplace_back(std::move(srf));
        };

        auto gap2srf = [&](const ep2& pb, const ep2& pe, double minh) {
            if (e_dlte(maxz, minh)) return;
            add_srf(pb, pe, minh, true);
        };

        auto pre = out.back();
        auto nth_srf = [&](const auto& cur) {
            int id{realId(pre, cur)};
            if (id <= c_virtual_srf) return true;
            double mh{minz};
            auto it = insts.find(id - 1);
            if (it != insts.end()) {
                if (it->second.integral) {
                    if (e_dlte(maxz, it->second.minH)) {
                        return true;
                    } else if (it->second.minH > minz) {
                        mh = it->second.minH;
                    }
                } else {
                    Isegs segs;
                    ep2 pb{pre.x(), pre.y()};
                    ep2 pe{cur.x(), cur.y()};
                    e_seg ls{pb, pe};
                    auto area{ls.len * (maxz - mh)};
                    for (const auto& u : it->second.units) {
                        if (e_dlte(u.minH, mh)) continue;
                        e_seg res{};
                        if (e_inpos(ls, u.seg, res)) {
                            double h = maxz < u.minH ? maxz : u.minH;
                            if (e_dlte(ls.len, res.len) && e_dlte(maxz, h)) return true;
                            area -= res.len * (h - mh);
                            if (area < c_5e_2) return true;
                            InnerSeg iseg;
                            iseg.norm = ls[2];
                            iseg.minh = h;
                            iseg.seg = res;
                            iseg.swap();
                            segs.emplace_back(std::move(iseg));
                        }
                    }
                    if (!segs.empty()) {
                        std::sort(segs.begin(), segs.end());
                        for (const auto& iseg : segs) {
                            if (!e_p_equal(pb, iseg.seg[0], 2, c_5e_2)) {
                                gap2srf(pb, iseg.seg[0], mh);
                            }
                            gap2srf(iseg.seg[0], iseg.seg[1], iseg.minh);
                            pb = iseg.seg[1];
                        }
                        if (!e_p_equal(pb, pe, 2, c_5e_2)) {
                            gap2srf(pb, pe, mh);
                        }
                        return true;
                    }
                }
            }
            add_srf(pre, cur, mh);
            return true;
        };

        srfs.resize(1);
        auto& roof = srfs[c_roof_Id];
        roof.id = c_roof_Id;
        roof._h = maxz;
        roof._pts.reserve(out.size());
        for (const auto& it : out) {
            srfs[c_roof_Id]._pts.emplace_back(it.x(), it.y(), maxz);
            if (!nth_srf(it)) return false;
            pre = it;
        }
        return true;
    }
};
static const GeneratePolyhedron g_GP;
#endif