#pragma once
#ifndef _BASE_H_
#define _BASE_H_
#include "env_.h"
#include <bitset>
inline double deg2rad(double deg) {
    static constexpr auto c_d2r{c_pi / 180.0};
    return deg * c_d2r;
}

inline void base_correct_sr(double& sr, double v) {
    sr = std::min(1.0, std::max(0.0, v));
}

inline double base_td(double val) {
    return std::round(val * 1e4) * 1e-4;
}

inline bg_point base_tdp(const ep3& p) {
    return bg_point{base_td(p.x()), base_td(p.y())};
}

inline double div_s(double x, double y) {
    return (y > -c_1e_8 && y < c_1e_8) ? 0.0 : x / y;
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

    Solar(double lat, double lon, double ls, int start = 0) {
        lat = deg2rad(lat);
        sinlat = std::sin(lat);
        coslat = std::cos(lat);
        longitude = lon - ls - 180.0;
        _start = start; 
    }

    void sv(int step);
};

struct ShadowSet {
    std::vector<int> ids;
    int day{};
    bool all{true};
    bool diffuse{};
    bool separate{};
    bool mt{true};

    inline bool have(int id) const {
        return all || std::binary_search(ids.begin(), ids.end(), id);
    }

    inline void dset(bool df = false, bool mesh = false) {
        diffuse = df;
        if (diffuse) {
            separate &= mesh;
        }
    }

    inline int type() const {
        return diffuse ? (separate ? 0 : 1) : 2; //0(no calc srf) 1(as one wall) 2(wall + window)
    }
};

struct Subsurface {
    Polys polys; //no holes
    double wwr{};
    double area{};

    inline bool exist() const {
        return wwr > c_1e_8 || !polys.empty();
    }

    inline bool full() const {
        return wwr > 1.0 - c_1e_8;
    }
};

struct DiffuseR {
    double shdgIsoSky{}; //-> reuse reflection_factor
    double rdome{}; //woShdgIsoSky rdome *= (1.0 + cos_g)
    double shdgHoriz{};
    double rhorizon{}; //woShdgHoriz rhorizon *= sin_g
    double shdgGround{};
    double rdomeG{}; //woShdgGround rdomeg *= 0.5 * (1.0 - cos_g)

    inline void add(double cosPhi, double sunCosTheta, double sr, int iphi, bool dome) {
        auto woShdg{cosPhi * sunCosTheta};
        auto withShdg{woShdg * (1.0 - sr)};
        if (dome) {
            shdgIsoSky += withShdg;
            rdome += woShdg;
            if (iphi == 0) {
                shdgHoriz += withShdg;
                rhorizon += woShdg;
            }
        } else {
            shdgGround += withShdg;
            rdomeG += woShdg;
        }        
    }

    inline DiffuseR& operator+=(const DiffuseR& rhs) {
        shdgIsoSky += rhs.shdgIsoSky;
        rdome += rhs.rdome;
        shdgHoriz += rhs.shdgHoriz;
        rhorizon += rhs.rhorizon;
        shdgGround += rhs.shdgGround;
        rdomeG += rhs.rdomeG;
        return *this;
    }

    inline void calc(double sin_g, double cos_g) {//reuse
        rdome = div_s(shdgIsoSky * (1.0 + cos_g), rdome);
        rhorizon = div_s(shdgHoriz * sin_g, rhorizon);
        rdomeG = div_s(shdgGround * 0.5 * (1.0 - cos_g), rdomeG);
        shdgIsoSky = 1.0 - 0.5 * rdome - rdomeG;
    }

    inline double factor() const {
        return shdgIsoSky;
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
        std::string str{};
        for (const auto& p : poly.outer()) {
            str += std::format("{} {} ", p.x(), p.y());
        }
        if (!str.empty()) {
            str.pop_back();
        }
        return str;
    }
};

struct Surface {
    int id{};
    int _dir{-1};
    double _d{};
    double _rotz{};
    double _h{};
    double _area{};
    double _opaque{};
    double _amwall{};
    double _amwin{};
    double _sin_g{};
    bool _gap{};
    bool _gcalc{};
    bool _scalc{};
    int _shdg{};
    mutable double _cos_s{};
    mutable double _sr{}; //shadow ratio -> (mutable) wall baseIntensity
    mutable double _wsr{};

    epts _pts;
    std::vector<epts> inners;
    Subsurface window;
    ep3 _normi{};
    ep3 _center{};
    e_rmat _rmat;
    bg_polygon _poly;
    Polys _lights;
    DiffuseR _dr;
    bg_segment gseg;
    std::vector<std::pair<double, double>> daysrs;

    bool init(const ep3& center);
    e_rmat rmat_inv() const;

    void wpoly(int floor);

    void reset();
    void setsr(double sr, double wsr, bool save);
    void shadowInit(double sh, const ep3& sv, bool save, bool clear);
    void calcsr(const ep3& sv, int pos);
private:
    bool equation(const ep3& c);
    void rotation();
    void polygon();
    bool direction();
};
using Surfaces = std::vector<Surface>;
using NBHDs = std::vector<const Surface*>;

struct Polyhedron {
    int _id{};
    int _buildId{};
    int floor{1};
    Surfaces srfs;
    ebbox box; //te=box(2) h=box(5)
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

    inline void setDir() {
        if (srf == nullptr) return;
        dir = (m.wall && srf->id != c_roof_Id) ? srf->_dir : srf->_dir + 4;
    }

    inline void shadow(int pos, bool clear, bool save, bool intersect = false) {
        if (srf->_shdg == c_shadow_other_set) {
            m.sr = daysrs[pos];
            return;
        }
        if (clear && save) {
            daysrs.set();
        }
        auto f = [&](bool s) {
            m.sr = s;
            if (!s && save) {
                daysrs.reset(pos);
            }
        };
        if (srf->_shdg == c_shadow_init) {
            if (!intersect || m.poly.outer().empty()) {
                f(!pointInPolys(m.c_rot, srf->_lights));
            } else {
                auto larea{intersectArea(m.poly, srf->_lights)};
                base_correct_sr(m.sr, 1.0 - larea / m.area);
            } 
        } else {
            f(srf->_shdg == c_shadow_full);
        }
    }
};

using SrfMeshs = std::map<int, std::vector<SrfMesh>>;
using OneMany = std::map<int, std::vector<int>>;
void base_generate_meshs(SrfMeshs&, const GPhds&, const OneMany&, const ShadowSet&, double, double);
#endif
