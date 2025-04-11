#include "shadow.h"
void Shadow::calc_(Solar& s, GPhds& phds, const ShadowSet& cset) const {
    if (cset.day <= 1) {
        cfun(s.h, s.v, phds, cset, false);
        return;
    }
    auto len{s._step - s._start};
    if (len < 0) return;
    if (len < 24) {
        if (len == 0) {
            clear(phds);
        }
        cfun(s.h, s.v, phds, cset, true);
        return;
    }
    auto day{len / 24};
    if (day < cset.day) {
        auto pos{len % 24};
        for (auto& phd : phds) {
            for (auto& [_, srf] : phd.srfs) {
                srf.calcsr(pos);
            }
        }
        return;
    }
    if (day == cset.day) {
        s._start = s._step;
        clear(phds);
        cfun(s.h, s.v, phds, cset, true);
    }
}

void Shadow::clear(GPhds& phds) const {
    for (auto& phd : phds) {
        for (auto& [_, srf] : phd.srfs) {
            srf.clearsr();
        }
    }
}

void Shadow::cfun(double sh, const ep3& sv, GPhds& phds, const ShadowSet& cset, bool save) const {
    init(sv, phds, save);
    neighbours(sh, sv, phds, cset);
    calculate(sv, phds, cset.ctype, cset.mode, save);
}

void Shadow::diffuse(GPhds& phds, SrfMeshs& smeshs, const ShadowSet& cset, bool alone) const {
    auto constexpr NPhi = 6;                      // Number of altitude angle steps for sky integration
    auto constexpr NTheta = 24;                   // Number of azimuth angle steps for sky integration
    auto constexpr DPhi = 0.5 * c_pi / NPhi;       // Altitude step size
    auto constexpr DTheta = 2.0 * c_pi / NTheta;   // Azimuth step size
    //auto constexpr DThetaDPhi = DTheta * DPhi;    // Product of DTheta and DPhi
    auto constexpr PhiMin = 0.5 * DPhi;           // Minimum altitude
    auto fd_fun = [](DiffuseR& fd, double cos_Phi, double sunCosTheta, double sr, int iPhi, bool dome) {
        auto woShdg{cos_Phi * sunCosTheta * -1.0};
        auto withShdg{woShdg * (1.0 - sr)};
        if (dome) {
            fd.withShdgIsoSky += withShdg;
            fd.rdome += woShdg;
            // Horizon region
            if (iPhi == 0) {
                fd.withShdgHoriz += withShdg;
                fd.rhorizon += woShdg;
            }
        } else {
            fd.withShdgGround += withShdg;
            fd.rdomeG += woShdg;
        }
    };
    if (alone && !smeshs.empty()) {
        cset.dset(true, 0);
    } else {
        cset.dset(true, 2);
    }
    ep3 sv;
    for (int IPhi = 0; IPhi < NPhi; ++IPhi) { // Loop over patch altitude values
        auto Phi = PhiMin + IPhi * DPhi;
        auto sh{Phi};
        sv(2) = -std::sin(Phi);
        auto cos_Phi{std::cos(Phi)};
        for (int ITheta = 0; ITheta < NTheta; ++ITheta) { // Loop over patch azimuth values
            auto Theta = ITheta * DTheta;
            auto cos_Theta{std::cos(Theta)};
            auto sin_Theta{std::sin(Theta)};
            sv(0) = cos_Phi * sin_Theta;
            sv(1) = cos_Phi * cos_Theta;
            
            auto srffun = [&](bool dome) {
                for (auto& phd : phds) {
                    for (auto&[_, srf] : phd.srfs) {
                        if (srf.id <= c_virtual_srf) continue;
                        auto sunCosTheta{srf._normi.dot(sv)};
                        if (sunCosTheta >= 0.0) continue;
                        fd_fun(srf._dr, cos_Phi, sunCosTheta, srf._sr, IPhi, dome);
                    }
                }         
            };
            auto meshfun = [&](bool dome) {
                for (auto& [_, meshs] : smeshs) {
                    for (auto& mesh : meshs) {
                        auto sunCosTheta{mesh.srf->_normi.dot(sv)};
                        if (sunCosTheta >= 0.0) continue;
                        mesh.shadow_(0, false, false, true);
                        fd_fun(mesh._dr, cos_Phi, sunCosTheta, mesh.m.sr, IPhi, dome);
                    }
                }
            };
            auto rfun = [&](bool dome) {
                cfun(sh, sv, phds, cset, false);
                if (alone && !smeshs.empty()) {
                    meshfun(dome);
                } else {
                    srffun(dome);
                }
            };

            rfun(true);
            sv(2) *= -1.0;
            rfun(false);
            sv(2) *= -1.0;
        }
    }

    auto f_cfun = [&]() {
        for (auto& phd : phds) {
            for (auto&[_, srf] : phd.srfs) {
                srf._dr.calc();
            }
        }
    };
    auto c_cfun = [&](bool alone) {
        for (auto& [_, meshs] : smeshs) {
            for (auto& mesh : meshs) {
                if (alone) {
                    mesh._dr.calc();
                } else {
                    mesh._dr.rdome = mesh.srf->_dr.rdome;
                    mesh._dr.rhorizon = mesh.srf->_dr.rhorizon;
                    mesh._dr.rdomeG = mesh.srf->_dr.rdomeG;
                }
            }
        }
    };

    if (alone && !smeshs.empty()) {
        c_cfun(alone);
    } else {
        f_cfun();
        c_cfun(alone);
    }
    cset.dset(false, 1);
}

void Shadow::init(const ep3& sv, GPhds& phds, bool save) const {
    for (auto& phd : phds) {
        phd.neighbours.clear();
        for (auto& [_, srf] : phd.srfs) {
            srf.shadowInit(sv, save);
        }
    }
}

void Shadow::neighbours(double sh, const ep3& sv, GPhds& phds, const ShadowSet& cset) const {
    auto tansh{std::tan(sh)};
    for (auto it = phds.begin(); it != phds.end(); ++it) {
        if (!cset.should(it->_buildId)) continue;
        auto jt{it};
        for (++jt; jt != phds.end(); ++jt) {
            if (it->_buildId == jt->_buildId) {
                add(it->neighbours, *jt);
                add(jt->neighbours, *it);
                continue;
            }
            auto dis{0.0};
            if (ndis(*it, *jt, sh, dis)) {
                if (check(it->box, jt->box, tansh, sv, dis, cset.diffuse)) {
                    add(it->neighbours, *jt);
                }
                if (!cset.should(jt->_buildId)) continue;
                if (check(jt->box, it->box, tansh, sv, dis, cset.diffuse)) {
                    add(jt->neighbours, *it);
                }
            }         
        }
    }
}

bool Shadow::ndis(const Polyhedron& phd1, const Polyhedron& phd2, double sh, double& dis) const {
    if (sh < c_1e_8) return false;
    auto hmax{std::max(phd1.box(5), phd2.box(5))};
    if (hmax < c_1e_8) return false; //underground
    ep3 pt{phd1._center - phd2._center};
    dis = std::sqrt(pt(1) * pt(1) + pt(0) * pt(0)) - phd1._dis - phd2._dis;
    auto threshold{hmax / sh};
    if (dis > threshold) return false;
    return true;
}

bool Shadow::check(const ebbox& b1, const ebbox& b2, const double tansh, 
    const ep3& sv, const double dis, bool noSolar) const {
    double h2{b2(5) / tansh};
    if (h2 < dis) return false;
    if (noSolar) return true;
    auto f = [&](int i) {
        int j{i + 3};
        if (sv(i) < -c_1e_8) return b1(i) < b2(j);
        else if (sv(i) > c_1e_8) return b2(i) < b1(j);
        return false;
    };
    return f(0) && f(1);
}

void Shadow::add(std::vector<const Surface*>& res, const Polyhedron& other) const {
    for (const auto& [_, srf] : other.srfs) {
        if (srf._shdg < c_shadow_full && srf.id > c_virtual_srf) {
            res.emplace_back(&srf);
        }
    }
}

void Shadow::calculate(const ep3& sv, GPhds& phds, int type, bool mode, bool save) const {
    for (auto& phd : phds) {
        auto noflag{phd.neighbours.empty()};
        for (auto& [_, srf] : phd.srfs) {
            if (srf._scalc) continue;
            if (noflag) {
                srf._shdg = c_shadow_none;
                if (save) {
                    srf.setsr(0.0, 0.0, save);
                }
                continue;
            }
            if (mode) {
                forward(srf, phd, sv, type, save);
            } else {
                backward(srf, phd, -sv, type, save);
            }
        }
    }
}

void Shadow::forward(Surface& srf, Polyhedron& phd, const ep3& sv, int type, bool save) const {
    auto ndot{srf._norm.dot(sv)};
    auto num{0};
    bg_polygon poly;
    auto fun = [&](const epts& ps, double fh) {
        if (!_projection(poly, srf, sv, ndot, ps, fh)) return false;
        auto sara{std::abs(boost::geometry::area(poly))};
        if (sara < 1e-1) return false;
        Polys out;
        if (!boost::geometry::intersection(srf._poly, poly, out) || out.empty()) return false;
        auto area{0.0};
        for (auto& o : out) {
            area += std::abs(boost::geometry::area(o));
        }
        if (area > sara || area > srf._area + c_1e_2) return false;
        ++num;
        if (difference(srf._lights, srf._poly, poly)) return true;
        return false;
    };
    for (const auto& it : phd.neighbours) {
        if (fun(it->_pts, it->_h)) break;
    }
    if (num == 0) {
        srf._shdg = c_shadow_none;
        srf.setsr(0.0, 0.0, save);
        return;
    }
    srfcalc(srf, type, save);
}

bool Shadow::_projection(bg_polygon& poly, const Surface& srf, const ep3& sv, 
    double ndot, const epts& pts, double fh) const {
    auto terrain = [&](double h) {
        if (h < fh - c_1e_2) {
            h = 0.0;
        }
        return h;
    };
    auto insert = [&](const ep3& p0, const ep3& p1) {
        ep3 dir{p0 - p1};
        dir(2) = terrain(p0(2)) - terrain(p1(2));
        const auto ddot{srf._norm.dot(dir)};
        const auto t{(srf._d - p0.dot(srf._norm)) / ddot};
        const ep3 sp{p0 + t * dir};
        const ep3 pt{srf._rmat * sp};
        poly.outer().emplace_back(pt.x(), pt.y());
    };
    poly.clear();
    ep3 sp;
    auto proj = [&](ep3 pt) {
        pt(2) = terrain(pt(2));
        auto before{false};
        const auto t{(srf._d - pt.dot(srf._norm)) / ndot};
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
            const ep3 pt{srf._rmat * sp};
            poly.outer().emplace_back(pt.x(), pt.y()); 
        }
        front = before;
        ++i;
    }
    return !poly.outer().empty();   
}

void Shadow::srfcalc(Surface& srf, int type, bool save) const {
    if (type <= 0) return;
    if (srf._scalc) return;
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
    auto f = [](Polys& res, const bg_polygon& p1, const bg_polygon& p2) {
        Polys out;
        boost::geometry::difference(p1, p2, out);
        for (const auto& pot : out) {
            if (boost::geometry::area(pot) < c_5e_2) continue;
            res.emplace_back(pot);
        }
    };
    if (res.empty()) {
        f(res, totalPoly, poly);
    } else {
        Polys tmp;
        tmp.swap(res);
        for (const auto& tpoly : tmp) {
            f(res, tpoly, poly);
        }
    }
    return res.empty();
}

double Shadow::substract(const bg_polygon& totalPoly, const Polys& intersectPolys) const {
    if (intersectPolys.empty()) {
        return std::abs(boost::geometry::area(totalPoly));
    }
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
    auto t{(nsrf->_d - nsrf->_norm.dot(p)) / nsrf->_norm.dot(sv)};
    ep3 ap{t * sv};
    if (std::abs(ap[0]) < c_5e_2 && std::abs(ap[1]) < c_5e_2 && std::abs(ap[2]) < c_5e_2) {
        return false;
    }
    if (ap.dot(sv) < c_5e_2) {
        return false;
    }
    ep3 sp = p + ap;
    const ep3 pt{nsrf->_rmat * sp};
    bg_point bp(pt(0), pt(1));
    if (boost::geometry::covered_by(bp, nsrf->_poly)) {
        return true;
    }
    return false;
}

void Shadow::backward(Surface& srf, Polyhedron& phd, const ep3& sv, int type, bool save) const {
    auto sarea{0.0};
    auto wsarea{0.0};
    for (auto& m : srf.meshs) {
        for (auto n : phd.neighbours) {
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
            } else {
                r = 0.0;
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
