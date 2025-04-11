#include "base.h"
const static ep3 c_ground(0.0, 0.0, 1.0);
int wallDir(const ep3& norm) {
    const static std::array cDir = {
        ep3(1.0, 0.0, 0.0), ep3(0.0, 1.0, 0.0), 
        ep3(0.0, -1.0, 0.0), ep3(-1.0, 0.0, 0.0)
    };//ensw
    static const auto r2{std::pow(2.0, 0.5) * 0.5};
    for (int i = 0; i < 4; ++i) {
        auto d{norm.dot(cDir[i])};
        if (d >= r2 && d <= 1.0) {
            return i;
        }
    }
    return -1;
}

inline double div_s(double x, double y) {
    if (y > -c_1e_8 && y < c_1e_8) {
        return 0.0;
    }
    return x / y;
}

void Solar::sv(int d, double hour, int step) {
    _step = step;
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
        if (w < 0) {
            azimuth *= -1.0;
        }
    }
    v << cosh * std::sin(azimuth), cosh * cosa, -sinh;
}

void Solar::init(double lat, double lon, double ls, int start) {
    lat = deg2rad(lat);
    sinlat = std::sin(lat);
    coslat = std::cos(lat);
    longitude = lon - ls - 180.0;
    _start = start; 
}

void base_gen_grids(Polys& grids, const bg_polygon& poly, double du, double dv) {
    bg_box box;
    boost::geometry::envelope(poly, box);
    auto umin{box.min_corner().x()};
    auto vmin{box.min_corner().y()};
    auto umax{box.max_corner().x()};
    auto vmax{box.max_corner().y()};
    auto u_div{static_cast<int>(std::ceil((umax - umin) / du))};
    auto v_div{static_cast<int>(std::ceil((vmax - vmin) / dv))};
    grids.reserve(u_div * v_div);
    for (int i = 0; i < u_div; ++i) {
        auto u1{umin + i * du};
        auto u2{u1 + du};
        for (int j = 0; j < v_div; ++j) {
            auto v1{vmin + j * dv};
            auto v2{v1 + dv};
            bg_point p1{u1, v1};
            bg_point p2{u2, v1};
            bg_point p3{u2, v2};
            bg_point p4{u1, v2};
            bg_polygon mp;
            mp.outer().emplace_back(p1);
            mp.outer().emplace_back(p2);
            mp.outer().emplace_back(p3);
            mp.outer().emplace_back(p4);
            grids.emplace_back(mp);
        }
    }
}

void base_gen_meshs(Meshs& meshs, const bg_polygon& poly, double du, double dv, bool wall) {
    Polys grids;
    base_gen_grids(grids, poly, du, dv);
    for (const auto& grid : grids) {
        Polys outs;
        if (!boost::geometry::intersection(poly, grid, outs) || outs.empty()) continue;
        for (auto& out : outs) {
            if (!boost::geometry::is_valid(out)) {
                boost::geometry::correct(out);
            }
            Mesh mesh;
            mesh.poly.outer().swap(out.outer());
            mesh.poly.inners().swap(out.inners());
            mesh.area = boost::geometry::area(mesh.poly);
            boost::geometry::centroid(mesh.poly, mesh.c_rot);
            //mesh.c = mat_inv * ep3(mesh.c_rot.x(), mesh.c_rot.y(), rz);
            mesh.wall = wall;
            meshs.emplace_back(mesh);
        }
    }    
}

bool Surface::init(const ep3& center) {
    if (id <= c_virtual_srf || _gcalc) return true;
    if (!equation(center)) {
        std::cout << "init equation error " << id << std::endl; 
        return false;
    }
    rotation();
    setPoly();
    setDir();
    _gcalc = true;
    return true;
}

bool Surface::equation(const ep3& c) {
    auto get_norm = [&](bool reverse = false) {
        auto size{static_cast<int>(_pts.size())};
        if (size < 3) return -1;
        auto i0{0};
        auto i1{1};
        auto i2{2};
        if (reverse) {
            i0 = size - 1;
            i1 = size - 2;
            i2 = size - 3;
        }
        ep3 p10{_pts[i1] - _pts[i0]};
        ep3 p20{_pts[i2] - _pts[i0]};
        _norm = p10.cross(p20);
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
    auto ret{get_norm()};
    if (ret < 0) {
        ret = get_norm(true);
        if (ret < 0) return false;
    }
    if (ret == 1) {
        std::reverse(_pts.begin(), _pts.end());
    }
    _d = _norm.dot(_pts[0]);
    return true;
}

void Surface::rotation() {
    _gdot = c_ground.dot(_normi);
    ep3 rotAxis{_normi.cross(c_ground)};
    if (rotAxis.norm() < c_1e_8) {
        _rmat = Eigen::Matrix3d::Identity();
        _rmat_inv = Eigen::Matrix3d::Identity();
        return;
    }
    rotAxis.normalize();
    auto angle{std::acos(_gdot)};
    Eigen::AngleAxisd raa(angle, rotAxis);
    _rmat = raa.matrix();
    _rmat_inv = _rmat.transpose();
    _sinsita = std::sqrt(1 - std::pow(_gdot, 2));
}

void Surface::setPoly() {
    auto rzset{false};
    auto addRing = [&](bg_ring& ring, const epts& pts) {
        for (const auto& pt : pts) {
            ep3 p{_rmat * pt};
            ring.emplace_back(p.x(), p.y());
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

void Surface::setDir() {
    if (id == c_roof_Id) {
        _dir = 4;
        bg_point bc(0.0, 0.0);
        boost::geometry::centroid(_poly, bc);
        _center = {bc.x(), bc.y(), _h};
    } else {
        _dir = wallDir(_normi);
        _center = (_pts[0] + _pts[2]) * 0.5;
    }
}

bool Surface::wpoly(int floor) {
    //only for sGP::stretch
    if (id <= c_virtual_srf) return true;
    if (!window.polys.empty()) return true;
    if (window.wwr > c_1e_8 && window.wwr < 1 - c_1e_8) {
        auto r{std::sqrt(window.wwr)};
        auto add = [&](const bg_polygon& poly) {
            auto&& wp{scalePoly(poly, r)};
            if (!boost::geometry::is_valid(wp)) {
                boost::geometry::correct(wp);
            }
            window.area += std::abs(boost::geometry::area(wp));
            window.polys.emplace_back(wp);
        };

        auto bottom = [&](ep3& minp, ep3& maxp, double& minz) {
            auto first{false};
            for (const auto& pt : _pts) {
                if (std::abs(_h - pt(2)) > c_5e_2) {
                    if (!first) {
                        minp = pt;
                        minz = pt(2);
                        first = true;
                    } else {
                        maxp = pt;
                        break;
                    }
                }
            }
        };
        if (floor == 1 || id == c_roof_Id || id >= c_split_realId) {
            add(_poly);
        } else {
            //srf only have 4 points
            ep3 minp(0, 0, 0);
            ep3 maxp(0, 0, 0);
            auto minz{0.0};
            bottom(minp, maxp, minz);
            auto h{(_h - minz) / floor};
            if (h < c_5e_2) return true;
            ep3 pt1{_rmat * minp};
            ep3 pt2{_rmat * maxp};
            for (int i = 0; i < floor; ++i) {
                bg_polygon poly;
                minp(2) += h;
                ep3 pt{_rmat * minp};
                poly.outer().emplace_back(pt.x(), pt.y());
                poly.outer().emplace_back(pt1.x(), pt1.y());
                poly.outer().emplace_back(pt2.x(), pt2.y());
                maxp(2) += h;
                pt2 = _rmat * maxp;
                poly.outer().emplace_back(pt2.x(), pt2.y());
                pt1 = pt;
                add(poly);
            }
        }
    }
    _opaque = _area - window.area;
    return true; 
}

void Surface::light2Json(json& j, bool all) const {
    auto ring2Json = [&](const bg_ring& ring, json& jr) {
        for (auto& bp : ring) {
            const ep3 ep{bp.x(), bp.y(), _rotz};
            ep3 realP = _rmat_inv * ep;
            json jp;
            for (int i = 0; i < 3; ++i) {
                jp.emplace_back(realP(i));
            }
            jr[jr.size()].swap(jp);
        }
        if (!jr.empty()) {
            jr.emplace_back(jr[0]);
        }
    };

    auto poly2Json = [&] (const bg_polygon& poly, json& jp) {
        auto& jpo = jp["o"];
        ring2Json(poly.outer(), jpo);
        if (!poly.inners().empty()) {
            auto& jpis = jp["ins"];
            for (size_t i = 0; i < poly.inners().size(); ++i) {
                auto& jpi = jpis[i];
                ring2Json(poly.inners()[i], jpi);
            }
        }
    };

    if (all) {
        json jp;
        poly2Json(_poly, jp);
        j[j.size()].swap(jp);
    } else if (!_lights.empty()) {
        for (const auto& poly : _lights) {
            json jp;
            poly2Json(poly, jp);
            j[j.size()].swap(jp);
        }
    }
}

void Surface::reset() {
    _scalc = false;
    _shdg = c_shadow_init;
    _sr = 0.0;
    _wsr = 0.0;
    _lights.clear();
}

double Surface::lightArea() const {
    double larea{};
    for (const auto& lpoly : _lights) {
        larea += std::abs(boost::geometry::area(lpoly));
    }
    return larea;
}

void Surface::setsr(double sr, double wsr, bool save) {
    _sr = sr;
    _wsr = wsr;
    base_correct_sr(_sr);
    base_correct_sr(_wsr);
    _scalc = true;
    if (save) {
        daysrs.emplace_back(std::pair{_sr, _wsr});
    }
}

bool Surface::calcsr(int pos) {
    if (pos < 0 || pos >= static_cast<int>(daysrs.size())) return false;
    auto [sr, wsr] = daysrs[pos];
    _sr = sr;
    _wsr = wsr;
    _shdg = c_shadow_other_set;
    return true;
}

bool Surface::shadowInit(const ep3& sv, bool save) {
    if (id <= c_virtual_srf)  return false;
    reset();
    auto dot{_norm.dot(sv)};
    if (dot > c_1e_8) {
        _shdg = c_shadow_full;
        setsr(1.0, 1.0, save);
        return true;
    }
    if (std::abs(dot) < 0.01) {
        _shdg = c_shadow_none;
        setsr(0.0, 0.0, save);
    }
    return false;
}

void Surface::draw(double du, double dv) {
    meshs.clear();
    auto wallPoly{_poly};
    for (const auto& wp : window.polys) {
        base_gen_meshs(meshs, wp, du, dv, false);
        wallPoly.inners().emplace_back(wp.outer());
    }
    if (!boost::geometry::is_valid(wallPoly)) {
        boost::geometry::correct(wallPoly);
    }
    base_gen_meshs(meshs, wallPoly, du, dv, true);
}

void Surface::draw(Meshs& meshs_, double du, double dv) const {
    auto wallPoly{_poly};
    for (const auto& wp : window.polys) {
        base_gen_meshs(meshs_, wp, du, dv, false);
        wallPoly.inners().emplace_back(wp.outer());
    }
    if (!boost::geometry::is_valid(wallPoly)) {
        boost::geometry::correct(wallPoly);
    }
    base_gen_meshs(meshs_, wallPoly, du, dv, true);
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
    for (const auto& [_, srf] : srfs) {
        extremum(srf._pts);
    }
    box << minp, maxp;
    _center = 0.5 * (minp + maxp);

    _dis = 0.0;
    for (auto& [_, srf] : srfs) {
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

void base_generateSrfMeshs(SrfMeshs& smeshs, const GPhds& phds, const OneMany& b2pids, 
    const ShadowSet& cset, double droof, double dwall) {
    Meshs meshs;
    for (auto bid : cset.ids) {
        auto& bmeshs = smeshs[bid];
        auto pit = b2pids.find(bid);
        if (pit == b2pids.end()) continue;
        for (auto id : pit->second) {
            if (id < 0 || id >= static_cast<int>(phds.size())) continue;
            auto& phd = phds.at(id);
            for (auto& [_, srf] : phd.srfs) {
                meshs.clear();
                if (srf.id == c_roof_Id) {
                    srf.draw(meshs, droof, droof);
                } else {
                    srf.draw(meshs, dwall, dwall);
                }
                for (auto& m : meshs) {
                    SrfMesh sm;
                    sm.phd = id;
                    sm.dir = srf._dir;
                    if (!m.wall) {
                        sm.dir += 4;
                    }
                    if (srf.id == c_roof_Id) {
                        sm.dir = 8;
                    }
                    sm.m = m;
                    sm.m.area *= 1e-3;
                    sm.srf = &srf;
                    bmeshs.emplace_back(sm);
                }
            }
        }
    }    
}