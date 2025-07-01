#include "base.h"
#include "tp.h"
int surfaceDir(const ep3& norm) {
    static constexpr std::array<std::array<int, 2>, 4> cDir = {{
        {0, 1}, {1, 1}, {1, -1}, {0, -1}
    }};//ensw
    static constexpr auto r2{0.70710678};
    for (int i = 0; i < 4; ++i) {
        auto& arr = cDir[i];
        auto d{norm(arr[0]) * arr[1]};
        if (d >= r2 && d <= 1.0) {
            return i;
        }
    }
    return -1;
}

void Solar::sv(int step) {
    _step = step;
    auto d{step / 24 + 1};
    auto hour{(step + 1) % 24 - 0.5};
    if (day != d) {
        day = d;
        auto declination{deg2rad(23.45 * std::sin((284.0 + day) * c_declination))};
        sinD = std::sin(declination);
        cosD = std::cos(declination);
        auto W{deg2rad(day)};
        e = (-0.0002786409 + 0.1227715 * std::cos(W + 1.498311)
            - 0.1654575 * std::cos(2 * W - 1.261546) 
            - 0.00535383 * std::cos(3 * W - 1.1571)) * 0.25 + longitude;
    }
    auto w{deg2rad(e + 15.0 * hour)};
    auto sinh{sinlat * sinD + coslat * cosD * std::cos(w)};
    h = std::asin(sinh);
    auto cosh{std::cos(h)};
    auto cosa{1.0};
    auto azimuth{0.0};
    if (std::abs(w) > c_1e_8) {
        cosa = div_s(sinh * sinlat - sinD, cosh * coslat);
        azimuth = std::acos(cosa);
        if (w < 0.0) {
            azimuth *= -1.0;
        }
    }
    v << cosh * std::sin(azimuth), cosh * cosa, -sinh;
}

bool Surface::init(const ep3& center) {
    if (_gcalc) return true;
    if (!equation(center)) return false;
    rotation();
    polygon();
    if (!direction()) return false;
    _gcalc = true;
    return true;
}

bool Surface::equation(const ep3& c) {
    auto normv = [&](int i0 = 0, int i1 = 1, int i2 = 2) {
        ep3 p10{_pts[i1] - _pts[i0]};
        ep3 p20{_pts[i2] - _pts[i0]};
        ep3 _norm = p10.cross(p20);
        ep3 l{0.5 * (_pts[i0] + _pts[i2]) - c};
        auto ret{0};
        if (_norm.dot(l) < -c_1e_8) {
            _norm *= -1.0;
            ret = 1;
        }
        _normi = _norm.normalized();
        if (_normi.norm() < c_5e_2) return -1;
        return ret;
    };

    auto size{static_cast<int>(_pts.size())};
    if (size < 3) return false;
    auto ret{normv()};
    if (ret < 0 && size > 3) {
        ret = normv(size - 3, size - 2, size - 1);
        if (ret < 0) return false;
    }
    if (ret == 1) {
        std::reverse(_pts.begin(), _pts.end());
    }
    _d = _normi.dot(_pts[0]);
    return true;
}

void Surface::rotation() {
    const static ep3 c_ground(0.0, 0.0, 1.0);
    ep3 rotAxis{_normi.cross(c_ground)};
    _sin_g = rotAxis.norm();
    if (_sin_g < c_1e_8) {
        _rmat = Eigen::Matrix3d::Identity();
        return;
    }
    rotAxis.normalize();
    auto angle{std::acos(_normi(2))};
    Eigen::AngleAxisd raa(angle, rotAxis);
    _rmat = raa.matrix();
}

e_rmat Surface::rmat_inv() const {
    e_rmat rmat_inv = _rmat.transpose();
    return rmat_inv;
}

void Surface::polygon() {
    auto rzset{false};
    auto addRing = [&](bg_ring& ring, const epts& pts) {
        for (auto& pt : pts) {
            ep3 p{_rmat * pt};
            ring.emplace_back(base_tdp(p));
            if (!rzset) {
                _rotz = p(2);
                rzset = true;
            }
        }
    };
    addRing(_poly.outer(), _pts);
    if (!inners.empty()) {
        _poly.inners().resize(inners.size());
        for (int i = 0; i < static_cast<int>(inners.size()); ++i) {
            addRing(_poly.inners()[i], inners[i]);
        }
    }
    if (!boost::geometry::is_valid(_poly)) {
        boost::geometry::correct(_poly);
    }
    _area = std::abs(boost::geometry::area(_poly));
}

bool Surface::direction() {
    if (id == c_roof_Id) {
        _dir = 4;
        bg_point bc(0.0, 0.0);
        boost::geometry::centroid(_poly, bc);
        _center = {bc.x(), bc.y(), _h};
    } else {
        _dir = surfaceDir(_normi);
        _center = (_pts[0] + _pts[2]) * 0.5;
    }
    return _dir >= 0;
}

void Surface::wpoly(int floor) { //only for g_GP::stretch
    if (!window.polys.empty()) return;
    if (window.wwr > c_1e_8 && window.wwr < 1 - c_1e_8) {
        auto r{std::sqrt(window.wwr)};
        auto add = [&](const bg_polygon& poly) {
            auto wp{scalePoly(poly, r)};
            if (!boost::geometry::is_valid(wp)) {
                boost::geometry::correct(wp);
            }
            window.area += std::abs(boost::geometry::area(wp));
            window.polys.emplace_back(wp);
        };
        auto bottom = [&](ep3& minp, ep3& maxp) {
            auto first{false};
            for (const auto& pt : _pts) {
                if (std::abs(_h - pt(2)) > c_5e_2) {
                    if (!first) {
                        minp = pt;
                        first = true;
                    } else {
                        maxp = pt;
                        break;
                    }
                }
            }
        };
        
        if (floor == 1 || id == c_roof_Id || _gap) {
            add(_poly);
        } else { //srf only have 4 points
            ep3 minp(0, 0, 0);
            ep3 maxp(0, 0, 0);
            bottom(minp, maxp);
            auto h{(_h - minp(2)) / floor};
            if (h < c_5e_2) return;
            ep3 pt1{_rmat * minp};
            ep3 pt2{_rmat * maxp};
            bg_polygon poly;
            poly.outer().resize(4);
            for (int i = 0; i < floor; ++i) {
                minp(2) += h;
                ep3 pt{_rmat * minp};
                poly.outer()[0] = base_tdp(pt);
                poly.outer()[1] = base_tdp(pt1);
                poly.outer()[2] = base_tdp(pt2);
                maxp(2) += h;
                pt2 = _rmat * maxp;
                poly.outer()[3] = base_tdp(pt2);
                pt1 = pt;
                add(poly);
            }
        }
    }
    if (window.full()) {
        window.area = _area;
    }
    _opaque = _area - window.area;
}

void Surface::reset() {
    _scalc = false;
    _shdg = c_shadow_init;
    _sr = 0.0;
    _wsr = 0.0;
    _lights.clear();
}

void Surface::setsr(double sr, double wsr, bool save) {
    base_correct_sr(_sr, sr);
    base_correct_sr(_wsr, wsr);
    _scalc = true;
    if (save) {
        daysrs.emplace_back(std::pair{_sr, _wsr});
    }
}

void Surface::calcsr(const ep3& sv, int pos) {
    _cos_s = _normi.dot(sv);
    std::tie(_sr, _wsr) = daysrs.at(pos);
    _shdg = c_shadow_other_set;
}

void Surface::shadowInit(double sh, const ep3& sv, bool save, bool clear) {
    if (clear) {
        daysrs.clear();
    }
    reset();
    _cos_s = _normi.dot(sv);
    if (sh <= 0.0 || _cos_s > 0.0) {
        _shdg = c_shadow_full;
        setsr(1.0, 1.0, save);
        return;
    }
    if (std::abs(_cos_s) < 0.01) {
        _shdg = c_shadow_none;
        setsr(0.0, 0.0, save);
    }
}

bool Polyhedron::init(int id, int buildId) {
    _id = id;
    _buildId = buildId;
    ep3 minp(0, 0, 0);
    ep3 maxp(0, 0, 0);
    auto f{true};
    auto extremum = [&](const epts& pts) {
        for (const auto& pt : pts) {
            if (f) {
                minp = pt;
                maxp = pt;
                f = false;
                continue;
            }
            for (int i = 0; i < 3; ++i) {
                if (minp(i) > pt(i)) {
                    minp(i) = pt(i);
                } else if (maxp(i) < pt(i)) {
                    maxp(i) = pt(i);
                }
            }
        }        
    };
    for (const auto& srf : srfs) {
        extremum(srf._pts);
    }
    box << minp, maxp;
    _center = 0.5 * (minp + maxp);
    for (auto& srf : srfs) {
        if (!srf.init(_center)) return false;
        for (const auto& pt : srf._pts) {
            ep3 p{_center - pt};
            auto dis{p(0) * p(0) + p(1) * p(1)};
            if (dis > _dis) {
                _dis = dis;
            }
        }
    }
    _dis = std::sqrt(_dis);
    return true;
}

void base_draw_meshs(std::vector<SrfMesh>& bmeshs, int phd, const Surface& srf, 
    const bg_polygon& poly, double du, double dv, bool wall) {
    auto generate = [&](const bg_polygon& grid) {
        Polys outs;
        if (!boost::geometry::intersection(poly, grid, outs) || outs.empty()) return;
        for (auto& out : outs) {
            if (!boost::geometry::is_valid(out)) {
                boost::geometry::correct(out);
            }
            SrfMesh sm;
            sm.phd = phd;
            sm.srf = &srf;
            auto& mesh = sm.m;
            std::swap(mesh.poly, out);
            mesh.area = boost::geometry::area(mesh.poly) * 1e-3;
            boost::geometry::centroid(mesh.poly, mesh.c_rot);
            mesh.wall = wall;
            sm.setDir();
            bmeshs.emplace_back(std::move(sm));
        }
    };
    bg_box box;
    boost::geometry::envelope(poly, box);
    auto umin{box.min_corner().x()};
    auto vmin{box.min_corner().y()};
    auto umax{box.max_corner().x()};
    auto vmax{box.max_corner().y()};
    auto u_div{static_cast<int>(std::ceil((umax - umin) / du))};
    auto v_div{static_cast<int>(std::ceil((vmax - vmin) / dv))};
    bg_polygon grid;
    auto& outer = grid.outer();
    outer.resize(4);
    for (int i = 0; i < u_div; ++i) {
        auto u1{umin + i * du};
        auto u2{u1 + du};
        for (int j = 0; j < v_div; ++j) {
            auto v1{vmin + j * dv};
            auto v2{v1 + dv};
            outer[0] = bg_point(u1, v1);
            outer[1] = bg_point(u2, v1);
            outer[2] = bg_point(u2, v2);
            outer[3] = bg_point(u1, v2);
            generate(grid);
        }
    }
}

void base_draw(std::vector<SrfMesh>& bmeshs, int phd, const Surface& srf, double du, double dv) {
    auto wallPoly{srf._poly};
    for (const auto& wp : srf.window.polys) {
        base_draw_meshs(bmeshs, phd, srf, wp, du, dv, false);
        wallPoly.inners().emplace_back(wp.outer());
    }
    if (!boost::geometry::is_valid(wallPoly)) {
        boost::geometry::correct(wallPoly);
    }
    base_draw_meshs(bmeshs, phd, srf, wallPoly, du, dv, true);
}

void base_generate_meshs(SrfMeshs& smeshs, const GPhds& phds, const OneMany& b2pids, 
    const ShadowSet& cset, double droof, double dwall) {
    for (auto id : cset.ids) {
        smeshs.emplace(id, std::vector<SrfMesh>{});
    }
    auto fun = [&](int i, int begin, int end) {
        for (i = begin; i < end; ++i) {
            auto bid = cset.ids.at(i);
            auto& pids = b2pids.at(bid);
            auto& bmeshs = smeshs.at(bid);
            for (auto id : pids) {
                auto& phd = phds.at(id);
                for (auto& srf : phd.srfs) {
                    if (srf.id == c_roof_Id) {
                        base_draw(bmeshs, id, srf, droof, droof);
                    } else {
                        base_draw(bmeshs, id, srf, dwall, dwall);
                    }
                }
            }
        }
    };
    asyncF(fun, cset.ids.size());
}
