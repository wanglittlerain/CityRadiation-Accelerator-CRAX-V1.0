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
    e_seg seg{};
    inline void swap() {
        ep2 dir{seg[0] - seg[1]};
        if (dir.dot(norm) < 0.0) {
            std::swap(seg[0], seg[1]);
        }
    }
    inline bool operator < (const InnerSeg& rhs) const {
        ep2 dir{rhs.seg[0] - seg[0]};
        return dir.dot(norm) < 0.0;
    }
};
using Isegs = std::vector<InnerSeg>;

class GeneratePolyhedron {
public:
    bool stretch(Polyhedron& phd, const auto& points, const auto& poly, 
        const Insts& insts, double maxz, double minz) const {
        const auto& out{poly.outer()};
        if (out.empty() || std::abs(maxz - minz) < c_5e_2 || minz > maxz) {
            return false;
        }

        auto p3l1 = [](const auto& p0, const auto& p1, const auto& p2) {
            double hd = cxd::handedness(p0, p1, p2);
            if (std::abs(hd) < c_1e_2) return true;
            return false;
        };
        auto size{static_cast<int>(points.size())};
        auto collinear = [&](const auto& p, int& self, int other) {
            auto pid{-1};
            auto nid{-1};
            for (int i = 0; i < size; ++i) {
                auto j{(i + 1) % size};
                auto& p0{points[i]};
                auto& p2{points[j]};
                if (p3l1(p0, p, p2)) {
                    pid = i;
                    nid = j;
                    break;
                }
            }
            if (pid < 0) return false;
            if (other == pid) {
                self = nid;
            } else if (other == nid) {
                self = pid;
            } else {
                return false;
            }
            return true;
        };

        auto real = [&](const auto& pre, const auto& cur, int& realId) {
            realId = c_virtual_srf;
            if (points.empty()) return false;
            auto fp{false};
            auto fc{false};
            auto pid{-1};
            auto cid{-1};
            for (int i = 0; i < size; ++i) {
                auto& pt{points[i]};
                if (!fp && e_p_equal(pt, pre, 2, c_5e_2)) {
                    pid = i;
                    fp = true;
                }
                if (!fc && e_p_equal(pt, cur, 2, c_5e_2)) {
                    cid = i;
                    fc = true;
                }
                if (fp && fc) break;
            }
            if (!fp && !fc) return false;
            if (!fp && !collinear(pre, pid, cid)) return false;
            if (!fc && !collinear(cur, cid, pid)) return false;
            auto dis{cid - pid};
            auto mid{cid};
            if (dis < 0) {
                dis *= -1;
                mid = pid;
            }
            if (dis == 1) {
                realId = mid + 1;
            } else if (dis == size - 1) {
                realId = 1;
            }
            return true;
        };

        auto& srfs = phd.srfs;
        auto add_srf = [&](int realId, const auto& pb, const auto& pe, double minh) {
            Surface srf;
            srf.id = realId;
            srf._h = maxz;
            ep3 p1(pb.x(), pb.y(), maxz);
            ep3 p2(pe.x(), pe.y(), minh);
            e_push_points(srf._pts, p1, p2);
            bg_point sp1{pb.x(), pb.y()};
            bg_point sp2{pe.x(), pe.y()};
            srf.gseg = bg_segment(sp1, sp2);
            srfs.emplace(realId, std::move(srf));
        };

        auto gap2srf = [&](const ep2& pb, const ep2& pe, double minh, int& sid) {
            if (e_dlte(maxz, minh)) return;
            add_srf(sid, pb, pe, minh);
            sid += c_split_realId;
        };

        auto& roof = srfs[c_roof_Id];
        roof.id = c_roof_Id;
        roof._h = maxz;

        int idx{};
        auto ptPre = out.back();
        auto nth_srf = [&](const auto& ptCur) {
            int realId{};
            double mh{minz};
            real(ptPre, ptCur, realId);
            auto it = insts.end();
            if (realId < 0) {
                realId = --idx;
            } else {
                it = insts.find(realId - 1);
                if (it != insts.end()) {
                    if (it->second.integral) {
                        if (e_dlte(maxz, it->second.minH)) {
                            realId = --idx;
                        } else if (it->second.minH > minz) {
                            mh = it->second.minH;
                        }
                    }
                }
            }
            if (realId <= c_virtual_srf) {
                return true;
            }
            auto area{0.0};
            if (it != insts.end() && !it->second.integral) {
                Isegs segs;
                ep2 pb{ptPre.x(), ptPre.y()};
                ep2 pe{ptCur.x(), ptCur.y()};
                e_seg ls{pb, pe};
                ep2 l{pb - pe};
                auto len{l.norm()};
                area += len * (maxz - mh);
                for (const auto& u : it->second.units) {
                    if (e_dlte(u.minH, mh)) continue;
                    ep2 ul = u.seg[0] - u.seg[1];
                    auto ulen{ul.norm()};
                    e_seg res{};
                    bool find{};
                    if (ulen < len) {
                        find = e_inpos(ls, u.seg, len, res);
                    } else {
                        find = e_inpos(u.seg, ls, ulen, res);
                    }
                    if (find) {
                        double h = maxz < u.minH ? maxz : u.minH;
                        ep2 dir = res[0] - res[1];
                        auto dis{dir.norm()};
                        if (e_dlte(len, dis) && e_dlte(maxz, h)) {
                            return true;
                        }
                        area -= dis * (h - mh);
                        if (area < c_5e_2) {
                            return true;
                        }
                        InnerSeg iseg;
                        iseg.norm = l;
                        iseg.minh = h;
                        iseg.seg = res;
                        iseg.swap();
                        segs.emplace_back(std::move(iseg));
                    }
                }
                if (!segs.empty()) {
                    std::sort(segs.begin(), segs.end());
                    int sid{realId + c_split_realId};
                    for (const auto& iseg : segs) {
                        auto& seg = iseg.seg;
                        if (!e_p_equal(pb, seg[0], 2, c_5e_2)) {
                            gap2srf(pb, seg[0], mh, sid);
                        }
                        gap2srf(seg[0], seg[1], iseg.minh, sid);
                        pb = seg[1];
                    }
                    if (!e_p_equal(pb, pe, 2, c_5e_2)) {
                        gap2srf(pb, pe, mh, sid);
                    }
                    return true;
                }
            }
            add_srf(realId, ptPre, ptCur, mh);
            return true;
        };
        for (const auto& it : out) {
            roof._pts.emplace_back(it.x(), it.y(), maxz);
            if (!nth_srf(it)) return false;
            ptPre = it;
        }
        return true;
    }
};
static const GeneratePolyhedron sGP;
#endif