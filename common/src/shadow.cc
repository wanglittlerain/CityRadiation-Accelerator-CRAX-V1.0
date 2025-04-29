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
            if (!cset.should(phd._buildId)) continue;
            for (auto& [_, srf] : phd.srfs) {
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
    auto f_nbhds = [&](Polyhedron& phd) {
        if (!cset.should(phd._buildId)) return false;
        if (!cset.mt) return true;
        for (const auto& nphd : phds) {
            if (phd._id == nphd._id) continue;
            if (phd._buildId == nphd._buildId) {
                add(phd.nbhds, nphd.srfs);
                continue;
            }
            auto dis{ndis(phd, nphd, tansh)};
            if (dis < nphd.box(5) && (cset.diffuse || check(phd.box, nphd.box, sv))) {
                add(phd.nbhds, nphd.srfs);
            }
        }
        return true;
    };
    auto f = [&](int begin, int end) {
        for (int i = begin; i < end; ++i) {
            auto& phd = phds.at(i);
            if (!f_nbhds(phd)) continue;
            auto noNbhd{phd.nbhds.empty()};
            for (auto& [_, srf] : phd.srfs) {
                if (srf._scalc) continue;
                if (noNbhd) {
                    srf._shdg = c_shadow_none;
                    if (save) {
                        srf.setsr(0.0, 0.0, save);
                    }
                    continue;
                }
                if (cset.mode) {
                    forward(srf, phd.nbhds, sv, cset.ctype, save);
                } else {
                    backward(srf, phd.nbhds, -sv, cset.ctype, save);
                }
            }
        }
    };
    if (cset.mt) {
        asyncF(f, phds.size());
    } else {
        neighborhoods(tansh, sv, phds, cset);
        f(0, static_cast<int>(phds.size()));
    }
}

void Shadow::diffuse(GPhds& phds, SrfMeshs& smeshs, const ShadowSet& cset, bool alone) const {
    auto constexpr NPhi = 6;                      // Number of altitude angle steps for sky integration
    auto constexpr NTheta = 24;                   // Number of azimuth angle steps for sky integration
    auto constexpr DPhi = 0.5 * c_pi / NPhi;       // Altitude step size
    auto constexpr DTheta = 2.0 * c_pi / NTheta;   // Azimuth step size
    auto constexpr PhiMin = 0.5 * DPhi;           // Minimum altitude
    auto fd_fun = [](DiffuseR& fd, double cos_Phi, double sunCosTheta, double sr, int iPhi, bool dome) {
        auto woShdg{cos_Phi * sunCosTheta * -1.0};
        auto withShdg{woShdg * (1.0 - sr)};
        if (dome) {
            fd.withShdgIsoSky += withShdg;
            fd.rdome += woShdg;
            if (iPhi == 0) {
                fd.withShdgHoriz += withShdg;
                fd.rhorizon += woShdg;
            }
        } else {
            fd.withShdgGround += withShdg;
            fd.rdomeG += woShdg;
        }
    };
    alone &= !smeshs.empty();
    if (alone) {
        cset.dset(true, 0);
    } else {
        cset.dset(true, 2);
    }
    ep3 sv;
    for (int IPhi = 0; IPhi < NPhi; ++IPhi) { // Loop over patch altitude values
        auto Phi = PhiMin + IPhi * DPhi;
        sv(2) = -std::sin(Phi);
        auto cos_Phi{std::cos(Phi)};
        for (int ITheta = 0; ITheta < NTheta; ++ITheta) { // Loop over patch azimuth values
            auto Theta = ITheta * DTheta;
            auto cos_Theta{std::cos(Theta)};
            auto sin_Theta{std::sin(Theta)};
            sv(0) = cos_Phi * sin_Theta;
            sv(1) = cos_Phi * cos_Theta;
            
            auto r_fun = [&](bool dome) {
                cfun(phds, Phi, sv, cset, false);
                if (alone) {
                    for (auto& [_, meshs] : smeshs) {
                        for (auto& mesh : meshs) {
                            auto sunCosTheta{mesh.srf->_cos_s};
                            if (sunCosTheta >= 0.0) continue;
                            mesh.shadow_(0, false, false, true);
                            fd_fun(mesh._dr, cos_Phi, sunCosTheta, mesh.m.sr, IPhi, dome);
                        }
                    }
                } else {
                    for (auto& phd : phds) {
                        for (auto&[_, srf] : phd.srfs) {
                            if (srf.id <= c_virtual_srf || srf._cos_s >= 0.0) continue;
                            fd_fun(srf._dr, cos_Phi, srf._cos_s, srf._sr, IPhi, dome);
                        }
                    }
                }
            };
            r_fun(true);
            sv(2) *= -1.0;
            r_fun(false);
            sv(2) *= -1.0;
        }
    }

    if (!alone) {
        for (auto& phd : phds) {
            for (auto&[_, srf] : phd.srfs) {
                srf._dr.calc();
            }
        }
    }
    for (auto& [_, meshs] : smeshs) {
        for (auto& mesh : meshs) {
            if (alone) {
                mesh._dr.calc();
                continue;
            }
            mesh._dr.rdome = mesh.srf->_dr.rdome;
            mesh._dr.rhorizon = mesh.srf->_dr.rhorizon;
            mesh._dr.rdomeG = mesh.srf->_dr.rdomeG;
        }
    }
    cset.dset(false, 1);
}

void Shadow::init(GPhds& phds, const double sh, const ep3& sv, bool save, bool clear) const {
    for (auto& phd : phds) {
        phd.nbhds.clear();
        for (auto& [_, srf] : phd.srfs) {
            srf.shadowInit(sh, sv, save, clear);
        }
    }
}

double Shadow::ndis(const Polyhedron& phd1, const Polyhedron& phd2, const double tansh) const {
    ep3 pt{phd1._center - phd2._center};
    pt(2) = 0.0;
    auto dis{pt.norm() - phd1._dis - phd2._dis};
    dis *= tansh;
    return dis;
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
    for (auto& [_, srf] : nsrfs) {
        if (srf._shdg < c_shadow_full && srf.id > c_virtual_srf) {
            res.emplace_back(&srf);
        }
    }
}

void Shadow::neighborhoods(double tansh, const ep3& sv, GPhds& phds, const ShadowSet& cset) const {
    for (auto it = phds.begin(); it != phds.end(); ++it) {
        auto need{cset.should(it->_buildId)};
        auto jt{it};
        for (++jt; jt != phds.end(); ++jt) {
            if (it->_buildId == jt->_buildId) {
                if (need) {
                    add(it->nbhds, jt->srfs);
                    add(jt->nbhds, it->srfs);
                }
                continue;
            }
            auto jneed{cset.should(jt->_buildId)};
            if (need || jneed) {
                auto dis{ndis(*it, *jt, tansh)};
                if (need && dis < jt->box(5) && (cset.diffuse || check(it->box, jt->box, sv))) {
                    add(it->nbhds, jt->srfs);
                }
                if (jneed && dis < it->box(5) && (cset.diffuse || check(jt->box, it->box, sv))) {
                    add(jt->nbhds, it->srfs);
                }
            }
        }
    }
}

void Shadow::forward(Surface& srf, const NBHDs& nbhds, const ep3& sv, int type, bool save) const {
    auto has{false};
    bg_polygon poly;
    auto fun = [&](const epts& ps, double fh) {
        if (!_projection(poly, srf, sv, ps, fh)) return false;
        auto sarea{std::abs(boost::geometry::area(poly))};
        if (sarea < 1e-1) return false;
        auto area{0.0};
        if (!intersectArea(srf._poly, poly, area)) return false;
        if (area < 1e-1 || area > std::min(sarea, srf._area)) return false;
        has = true;
        if (difference(srf._lights, srf._poly, poly)) return true;
        return false;
    };
    for (auto n : nbhds) {
        if (fun(n->_pts, n->_h)) break;
    }
    if (!has) {
        srf._shdg = c_shadow_none;
        srf.setsr(0.0, 0.0, save);
        return;
    }
    srfcalc(srf, type, save);
}

bool Shadow::_projection(bg_polygon& poly, const Surface& srf, const ep3& sv, 
    const epts& pts, double fh) const {
    auto terrain = [&](double& h) {
        if (h < fh - c_1e_2) {
            h = 0.0;
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
        poly.outer().emplace_back(base_td(pt.x()), base_td(pt.y()));
    };
    poly.clear();
    ep3 sp;
    auto proj = [&](ep3 pt) {
        terrain(pt(2));
        auto before{false};
        auto t{(srf._d - pt.dot(srf._normi)) / srf._cos_s};
        if (t > c_1e_8) {
            before = true;
            sp = pt + t * sv;
        }
        return before;        
    };
    auto front{proj(pts.back())};
    int i{};
    for (const auto& it : pts) {
        auto before{proj(it)};
        if (front != before) {
            if (i == 0) {
                insert(pts.back(), it);
            } else {
                insert(pts[i - 1], it);
            }
        }
        if (before) {
            ep3 pt{srf._rmat * sp};
            poly.outer().emplace_back(base_td(pt.x()), base_td(pt.y())); 
        }
        front = before;
        ++i;
    }
    return !poly.outer().empty();   
}

void Shadow::srfcalc(Surface& srf, int type, bool save) const {
    if (type <= 0 || srf._scalc) return;
    auto larea{srf.lightArea()};
    auto asr{1.0 - larea / srf._area};
    if (type == 2) {
        srf.setsr(asr, asr, save);
        return;
    }
    auto sr{asr};
    auto wsr{0.0};
    if (!srf.window.polys.empty()) {
        auto swinArea{0.0};
        for (const auto& winPoly : srf.window.polys) {
            swinArea += substract(winPoly, srf._lights);
        }
        if (srf.window.area > c_1e_8) {
            wsr = swinArea / srf.window.area;
        }
        if (srf._opaque < c_1e_8) {
            sr = 0.0;
        } else {
            auto osa{srf._area - larea - swinArea};
            if (osa < c_1e_8) {
                osa = 0.0;
            }
            sr = osa / srf._opaque;
        }
    } else {
        if (srf.window.wwr > 1.0 - c_1e_8) {
            wsr = asr;
            sr = 0.0;
        }
    }
    srf.setsr(sr, wsr, save);
}

bool Shadow::difference(Polys& res, const bg_polygon& totalPoly, const bg_polygon& poly) const {
    auto df = [](Polys& res, const bg_polygon& p1, const bg_polygon& p2) {
        Polys out;
        boost::geometry::difference(p1, p2, out);
        for (const auto& pot : out) {
            if (boost::geometry::area(pot) < 1e-1) continue;
            res.emplace_back(pot);
        }
    };
    if (res.empty()) {
        df(res, totalPoly, poly);
    } else {
        Polys tmp;
        tmp.swap(res);
        for (const auto& tpoly : tmp) {
            df(res, tpoly, poly);
        }
    }
    return res.empty();
}

double Shadow::substract(const bg_polygon& totalPoly, const Polys& intersectPolys) const {
    if (intersectPolys.empty()) return std::abs(boost::geometry::area(totalPoly));
    Polys res;
    for (const auto& poly : intersectPolys) {
        if (difference(res, totalPoly, poly)) break;
    }
    auto area{0.0};
    for (const auto& poly : res) {
        area += std::abs(boost::geometry::area(poly));
    }
    return area;
}

bool Shadow::projection_(const Surface* nsrf, const ep3& sv, const ep3& p) const {
    if (nsrf == nullptr) return false;
    auto t{(nsrf->_d - nsrf->_normi.dot(p)) / nsrf->_cos_s};
    ep3 ap{t * sv};
    if (std::abs(ap[0]) < c_5e_2 && std::abs(ap[1]) < c_5e_2 && std::abs(ap[2]) < c_5e_2) return false;
    if (ap.dot(sv) < c_5e_2) return false;
    ep3 sp = p + ap;
    ep3 pt{nsrf->_rmat * sp};
    bg_point bp(pt(0), pt(1));
    if (boost::geometry::covered_by(bp, nsrf->_poly)) {
        return true;
    }
    return false;
}

void Shadow::backward(Surface& srf, const NBHDs& nbhds, const ep3& sv, int type, bool save) const {
    auto sarea{0.0};
    auto wsarea{0.0};
    for (auto& m : srf.meshs) {
        for (auto n : nbhds) {
            ep3 c = srf._rmat_inv * ep3(m.c_rot.x(), m.c_rot.y(), srf._rotz); //m.c
            if (!projection_(n, sv, c)) {
                m.sr = 1.0;
                continue;
            }
            m.sr = 0.0;
            if (type > 0) {
                if (m.wall) {
                    sarea += m.area;
                } else {
                    wsarea += m.area;
                }
            }
            break;
        }
    }
    if (type > 0) {
        auto fun = [](double sarea, double tarea, double& r){
            if (tarea > c_1e_8) {
                r = sarea / tarea;
            }
        };
        auto sr{0.0};
        auto wsr{0.0};
        fun(sarea, srf._opaque, sr);
        fun(wsarea, srf.window.area, wsr);
        if (type == 1) {
            srf.setsr(sr, wsr, save);
        } else {
            fun(sarea + wsarea, srf._area, sr);
            srf.setsr(sr, sr, save);
        }
    }
}
