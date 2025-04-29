#pragma once
#ifndef _BASE_H_
#define _BASE_H_
#include "env_.h"
#include <bitset>
inline double deg2rad(double deg) {
    static constexpr auto c_d2r{c_pi / 180.0};
    return deg * c_d2r;
}

inline void base_correct_sr(double& sr) {
    if (sr > 1.0) {
        sr = 1.0;
    } else if (sr < 0.0) {
        sr = 0.0;
    }
}

inline double base_td(double val) {
    return std::round(val * 1e4) * 1e-4;
}

struct Solar {
    ep3 v;
    double h{};

    double longitude{};
    double sinlat{};
    double coslat{};
    double e{};
    double sinD{};
    double cosD{};

    int day{-1};
    int _start{};
    int _step{};

    inline void init(double lat, double lon, double ls, int start = 0) {
        lat = deg2rad(lat);
        sinlat = std::sin(lat);
        coslat = std::cos(lat);
        longitude = lon - ls - 180.0;
        _start = start; 
    }
    void sv(int d, double hour, int step);
};

using NameIds = std::map<std::string, int>;
struct ShadowSet {
    NameIds nameIds;
    std::vector<int> ids;
    int day{};
    bool all{true};
    bool mode{true};
    bool mt{true};
    mutable bool diffuse{};
    mutable int ctype{1}; //0(no calc srf) 1(calc and separate) 2(calc as one wall)

    inline bool should(int id) const {
        if (!all && !std::binary_search(ids.begin(), ids.end(), id)) return false;
        return true;
    }

    inline void dset(bool diff, int type) const {
        diffuse = diff;
        ctype = type;
    }
};

struct Subsurface {
    Polys polys; //no holes
    double wwr{};
    double area{};

    inline bool exist() const {
        return wwr > c_1e_8 || !polys.empty();
    }
};

struct WWRatio {
    std::array<double, 5> rs{};
    inline void set(Subsurface& sf, int dir) const {
        sf.wwr = rs[dir];
        base_correct_sr(sf.wwr);
    }
};

struct DiffuseR {
    double withShdgIsoSky{};
    double rdome{};
    double withShdgHoriz{};
    double rhorizon{};
    double withShdgGround{};
    double rdomeG{};

    void calc() {
        auto fun = [](double a, double b) {
            if (std::abs(b) > c_1e_8) {
                return a / b;
            } else {
                return a / (b + c_1e_8);
            }
        };
        rdome = fun(withShdgIsoSky, rdome);
        rhorizon = fun(withShdgHoriz, rhorizon);
        rdomeG = fun(withShdgGround, rdomeG);
    }
};

struct Mesh {
    //ep3 c; //if shadow backward use center  To save space
    bg_point c_rot{};
    bg_polygon poly; //no holes
    double area{};
    double sr{};
    bool wall{};

    inline std::string poly2str() const {
        std::string pstr{};
        for (const auto& p : poly.outer()) {
            if (!pstr.empty()) {
                pstr += " ";
            }
            auto val{std::format("{} {}", p.x(), p.y())};
            pstr += val;
        }
        return pstr;
    }
};
using Meshs = std::vector<Mesh>;
void base_gen_meshs(Meshs& meshs, const bg_polygon& poly, double du, double dv, bool wall);

struct Surface {
    int id{};
    int _dir{-1};
    double _d{};
    double _rotz{};
    double _h{};
    double _area{};
    double _opaque{};
    double _sin_g{};
    bool _gcalc{};

    bool _scalc{};
    int _shdg{};
    double _cos_s{};
    mutable double _sr{}; //shadow ratio -> (mutable) wall total radiation 
    mutable double _wsr{};

    epts _pts;
    std::vector<epts> inners;
    Subsurface window;
    ep3 _normi;
    ep3 _center{};
    e_rmat _rmat;
    e_rmat _rmat_inv;
    bg_polygon _poly;
    Polys _lights;
    DiffuseR _dr;
    bg_segment gseg;
    Meshs meshs;
    std::vector<std::pair<double, double>> daysrs;

    bool init(const ep3& center);
    bool wpoly(int floor);

    void reset();
    void setsr(double sr, double wsr, bool save);
    void shadowInit(double sh, const ep3& sv, bool save, bool clear);

    double lightArea() const;

    bool calcsr(const ep3& sv, int pos);
    void draw(double du, double dv);
    void draw(Meshs& meshs_, double du, double dv) const;
private:
    bool equation(const ep3& c);
    void rotation();
    void polygon();
    void direction();
};
using Surfaces = std::map<int, Surface>;
using NBHDs = std::vector<const Surface*>;

struct Polyhedron {
    int _id{};
    int _buildId{};
    int floor{1};
    Surfaces srfs;
    ebbox box;
    ep3 _center;
    double _dis{};
    NBHDs nbhds;
    bool init(int id, int buildId);
};
using GPhds = std::vector<Polyhedron>;

struct SrfMesh {
    int phd{};
    int dir{};
    Mesh m;
    const Surface* srf{};
    DiffuseR _dr{};
    std::bitset<24> daysrs;

    inline void shadow_(int pos, bool clear, bool save, bool intersect = false) {
        if (srf->_shdg == c_shadow_other_set) {
            m.sr = daysrs[pos];
            return;
        }
        if (clear && save) {
            daysrs.set();
        }
        auto s{false};
        if (srf->_shdg == c_shadow_full) {
            m.sr = 1.0;
            s = true;
        } else if (srf->_shdg == c_shadow_none) {
            m.sr = 0.0;
        } else {
            if (!intersect || m.poly.outer().empty()) {
                if (!pointInPolys(m.c_rot, srf->_lights)) {
                    m.sr = 1.0;
                    s = true;
                } else {
                    m.sr = 0.0;
                }
            } else {
                auto larea{intersectArea(m.poly, srf->_lights)};
                m.sr = 1.0 - larea / m.area;
                base_correct_sr(m.sr);
            }
        }
        if (save && !s) {
            daysrs.reset(pos);
        }
    }
};
using SrfMeshs = std::map<int, std::vector<SrfMesh>>;
using OneMany = std::map<int, std::vector<int>>;

void base_generateSrfMeshs(SrfMeshs&, const GPhds&, const OneMany&, const ShadowSet&, double, double);
#endif