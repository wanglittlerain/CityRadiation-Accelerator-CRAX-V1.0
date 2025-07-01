#include "shadow.h"
#include "tp.h"
void Shadow::calculate(GPhds& phds, Solar& s, const ShadowSet& cset) const {
    if (cset.day <= 1) {
        cfun(phds, s.h, s.v, cset, false);
        return;
    }
    auto len{s._step - s._start};
    if (len < 0) return;
    if (len < 24) {
        cfun(phds, s.h, s.v, cset, true, len == 0);
        return;
    }
    auto day{len / 24};
    if (day < cset.day) {
        auto pos{len % 24};
        for (auto& phd : phds) {
            if (!cset.have(phd._buildId)) continue;
            for (auto& srf : phd.srfs) {
                srf.calcsr(s.v, pos);
            }
        }
        return;
    }
    if (day == cset.day) {
        s._start = s._step;
        cfun(phds, s.h, s.v, cset, true, true);
    }
}

void Shadow::cfun(GPhds& phds, double sh, const ep3& sv, const ShadowSet& cset, bool save, bool clear) const {
    init(phds, sh, sv, save, clear);
    if (sh <= 0.0) return;
    auto tansh{std::tan(sh)};
    auto mt{cset.diffuse && cset.mt};
    auto type{cset.type()};
    auto f_nbhds = [&](Polyhedron& phd) {
        if (!mt) return;
        for (const auto& nphd : phds) {
            if (phd._id == nphd._id) continue;
            if (phd._buildId == nphd._buildId) {
                add(phd.nbhds, nphd.srfs);
                continue;
            }
            auto dis{ndis(phd, nphd, tansh)};
            checkAdd(phd, nphd, sv, dis, cset.diffuse);
        }
    };
    auto f = [&](int i, int begin, int end) {
        for (i = begin; i < end; ++i) {
            auto& phd = phds.at(i);
            if (!cset.have(phd._buildId)) continue;
            f_nbhds(phd);
            for (auto& srf : phd.srfs) {
                if (srf._scalc) continue;
                if (phd.nbhds.empty()) {
                    srf._shdg = c_shadow_none;
                    if (save) {
                        srf.setsr(0.0, 0.0, save);
                    }
                    continue;
                }
                forward(srf, phd.box(2), phd.nbhds, sv, type, save);
            }
        }
    };
    
    if (mt) {
        asyncF(f, phds.size());
    } else {
        neighborhoods(tansh, sv, phds, cset);
        f(0, 0, static_cast<int>(phds.size()));
    }
}

void Shadow::diffuse(GPhds& phds, SrfMeshs& smeshs, const ShadowSet& cset, int begin, int end) const {
    auto constexpr DPhi{0.5 * c_pi / c_phi};       // Altitude step size
    auto constexpr DTheta{2.0 * c_pi / c_theta};   // Azimuth step size
    auto constexpr PhiMin{0.5 * DPhi};           // Minimum altitude
    ep3 sv;
    auto phi{0.0};
    auto cosPhi{0.0};
    auto f_phi = [&](int iphi) {
        phi = PhiMin + iphi * DPhi;
        sv(2) = -std::sin(phi);
        cosPhi = std::cos(phi);
    };
    auto f_dr = [&](int iphi, bool dome, bool last = false) {
        cfun(phds, phi, sv, cset, false);
        if (cset.separate) {
            for (auto& [_, meshs] : smeshs) {
                for (auto& mesh : meshs) {
                    if (mesh.srf->_cos_s < 0.0) {
                        mesh.shadow(0, false, false, true);
                        mesh._dr.add(cosPhi, mesh.srf->_cos_s, mesh.m.sr, iphi, dome);
                    }
                    if (last) {
                        mesh._dr.calc(mesh.srf->_sin_g, mesh.srf->_normi(2));
                    }
                }
            }
            return;
        }
        for (auto& phd : phds) {
            for (auto& srf : phd.srfs) {
                if (srf._cos_s < 0.0) {
                    srf._dr.add(cosPhi, srf._cos_s, srf._sr, iphi, dome);
                }
                if (last) {
                    srf._dr.calc(srf._sin_g, srf._normi(2));
                }
            }
        }
    };
    auto loop = [&](int iphi, int itheta, bool last = false) {
        auto theta{itheta * DTheta};
        auto cosTheta{std::cos(theta)};
        auto sinTheta{std::sin(theta)};
        sv(0) = cosPhi * sinTheta;
        sv(1) = cosPhi * cosTheta;
        f_dr(iphi, true);
        sv(2) *= -1.0;
        f_dr(iphi, false, last);
        sv(2) *= -1.0;
    };

    if (cset.mt) {
        for (int i = begin; i < end; ++i) {
            auto iphi{i / c_theta};
            auto itheta{i % c_theta};
            f_phi(iphi);
            loop(iphi, itheta);
        }
    } else {
        for (int iphi = 0; iphi < c_phi; ++iphi) { // Loop over patch altitude values
            f_phi(iphi);
            for (int itheta = 0; itheta < c_theta; ++itheta) { // Loop over patch azimuth values
                auto last{iphi == c_phi - 1 && itheta == c_theta - 1};
                loop(iphi, itheta, last);
            }
        }
        if (cset.separate) return;
        for (auto& [_, meshs] : smeshs) {
            for (auto& mesh : meshs) {
                mesh._dr = mesh.srf->_dr;
            }
        }
    }
}

void Shadow::init(GPhds& phds, const double sh, const ep3& sv, bool save, bool clear) const {
    for (auto& phd : phds) {
        phd.nbhds.clear();
        for (auto& srf : phd.srfs) {
            srf.shadowInit(sh, sv, save, clear);
        }
    }
}

double Shadow::ndis(const Polyhedron& phd1, const Polyhedron& phd2, double tansh) const {
    ep3 pt{phd1._center - phd2._center};
    pt(2) = 0.0;
    auto dis{pt.norm() - phd1._dis - phd2._dis};
    return dis * tansh;
}

bool Shadow::check(const ebbox& b1, const ebbox& b2, const ep3& sv) const {
    auto cf = [&](int i) {
        int j{i + 3};
        if (sv(i) < -c_1e_8) return b1(i) < b2(j);
        else if (sv(i) > c_1e_8) return b2(i) < b1(j);
        return false;
    };
    return cf(0) && cf(1);
}

void Shadow::add(NBHDs& res, const Surfaces& nsrfs) const {
    for (auto& srf : nsrfs) {
        if (srf._shdg < c_shadow_full) {
            res.emplace_back(&srf);
        }
    }
}

void Shadow::checkAdd(Polyhedron& phd, const Polyhedron& nphd, const ep3& sv, double dis, bool dfuse) const {
    auto relative{nphd.box(5) - phd.box(2)};
    if (dfuse && phd.box(2) > 0.0) {
        relative = nphd.box(5);
    }
    if (dis < relative && (dfuse || check(phd.box, nphd.box, sv))) {
        add(phd.nbhds, nphd.srfs);
    }
}

void Shadow::neighborhoods(double tansh, const ep3& sv, GPhds& phds, const ShadowSet& cset) const {
    for (auto it = phds.begin(); it != phds.end(); ++it) {
        auto have{cset.have(it->_buildId)};
        auto jt{it};
        for (++jt; jt != phds.end(); ++jt) {
            if (it->_buildId == jt->_buildId) {
                if (have) {
                    add(it->nbhds, jt->srfs);
                    add(jt->nbhds, it->srfs);
                }
                continue;
            }
            auto jhave{cset.have(jt->_buildId)};
            if (!have && !jhave) continue;
            auto dis{ndis(*it, *jt, tansh)};
            if (have) {
                checkAdd(*it, *jt, sv, dis, cset.diffuse);
            }
            if (jhave) {
                checkAdd(*jt, *it, sv, dis, cset.diffuse);
            } 
        }
    }
}

void Shadow::forward(Surface& srf, double te, const NBHDs& nbhds, const ep3& sv, int type, bool save) const {
    auto has{false};
    bg_polygon poly;
    te = std::min(0.0, te);
    for (auto n : nbhds) {
        if (n->_h < te) continue;
        if (!_projection(poly, srf, sv, n->_pts, n->_h, te)) continue;
        auto sarea{std::abs(boost::geometry::area(poly))};
        if (sarea < 1e-1) continue;
        auto area{0.0};
        if (!intersectArea(srf._poly, poly, area)) continue;
        if (area < 1e-1 || area > std::min(sarea, srf._area)) continue;
        has = true;
        if (difference(srf._lights, srf._poly, poly)) break;
    }
    if (!has) {
        srf._shdg = c_shadow_none;
        srf.setsr(0.0, 0.0, save);
        return;
    }
    srfcalc(srf, type, save);
}

bool Shadow::_projection(bg_polygon& poly, const Surface& srf, 
    const ep3& sv, const epts& pts, double fh, double te) const {
    auto terrain = [&](double& h) {
        if (h < fh - c_1e_2) {
            h = te;
        }
    };
    auto insert = [&](ep3 p0, ep3 p1) {
        terrain(p0(2));
        terrain(p1(2));
        ep3 dir{p0 - p1};
        auto ddot{srf._normi.dot(dir)};
        auto t{(srf._d - p0.dot(srf._normi)) / ddot};
        ep3 sp{p0 + t * dir};
        ep3 pt{srf._rmat * sp};
        poly.outer().emplace_back(base_tdp(pt));
    };
    ep3 sp;
    auto proj = [&](ep3 pt) {
        terrain(pt(2));
        auto t{srf._d - pt.dot(srf._normi)};
        if (t < 0.0) {
            t /= srf._cos_s;
            sp = pt + t * sv;
            return true;
        }
        return false;
    };
    poly.clear();
    auto front{proj(pts.back())};
    for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
        auto before{proj(pts[i])};
        if (front != before) {
            if (i == 0) {
                insert(pts.back(), pts[i]);
            } else {
                insert(pts[i - 1], pts[i]);
            }
        }
        if (before) {
            ep3 pt{srf._rmat * sp};
            poly.outer().emplace_back(base_tdp(pt)); 
        }
        front = before;
    }
    return poly.outer().size() > 2;
}

void Shadow::srfcalc(Surface& srf, int type, bool save) const {
    if (type <= 0 || srf._scalc) return;
    auto larea{polysArea(srf._lights)};
    auto sr{1.0 - larea / srf._area};
    if (type == 1) {
        srf.setsr(sr, sr, save);
        return;
    }
    auto wsr{0.0};
    if (!srf.window.polys.empty()) {
        auto swinArea{0.0};
        for (const auto& winPoly : srf.window.polys) {
            swinArea += substract(winPoly, srf._lights);
        }
        wsr = div_s(swinArea, srf.window.area);
        sr = div_s(std::max(srf._area - larea - swinArea, 0.0), srf._opaque);
    } else if (srf.window.full()) {
        wsr = sr;
        sr = 0.0;
    }
    srf.setsr(sr, wsr, save);
}

bool Shadow::difference(Polys& res, const bg_polygon& totalPoly, const bg_polygon& poly) const {
    auto df = [&](const bg_polygon& p1) {
        Polys out;
        boost::geometry::difference(p1, poly, out);
        for (auto& p : out) {
            if (boost::geometry::area(p) < 1e-1) continue;
            res.emplace_back(std::move(p));
        }
    };
    if (res.empty()) {
        df(totalPoly);
    } else {
        Polys tmp;
        tmp.swap(res);
        for (const auto& p : tmp) {
            df(p);
        }
    }
    return res.empty();
}

double Shadow::substract(const bg_polygon& totalPoly, const Polys& intersectPolys) const {
    if (intersectPolys.empty()) return std::abs(boost::geometry::area(totalPoly));
    Polys res;
    for (auto& poly : intersectPolys) {
        if (difference(res, totalPoly, poly)) break;
    }
    return polysArea(res);
}
