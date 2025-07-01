#pragma once
#ifndef _CONFIG_H_
#define _CONFIG_H_
#include "base.h"
static constexpr std::string_view c_srf_name{"srf{}"};
struct Location {
    double latitude{};
    double longitude{};
    double ls{}; //timeZone * 15
};

struct Weather {
    double dnr{};
    double dhr{};
    std::string date{};
};
using Weathers = std::vector<Weather>;

struct GeoConf {
    int floor{};
    double te{};
    double h{};
    std::string name{};
    std::vector<e_seg> segs;

    inline void pushSegs(const bg_polygon& poly) {
        if (poly.outer().empty()) return;
        segs.reserve(poly.outer().size());
        auto pre{poly.outer().back()};
        for (const auto& cur : poly.outer()) {
            segs.emplace_back(e_seg{ep2{pre.x(), pre.y()}, ep2{cur.x(), cur.y()}});
            pre = cur;
        }
    }
};

struct WWRatio {
    std::array<double, 5> rs{};

    inline void set(Subsurface& sf, int dir) const {
        base_correct_sr(sf.wwr, rs.at(dir));
    }
};
using WWRs = std::map<int, WWRatio>;

struct Materials{
    double G_win{}; //1.0-conf.G_win
    double r_roof{};
    double r_wall{};
};
using BMaterials = std::map<int, Materials>;

struct BuildingFlag {
    std::string name{};
    double te{};
};
using BFlags = std::map<int, BuildingFlag>;
using NameIds = std::map<std::string, int>;

class Config {
public:
    bool readWeather(const std::filesystem::path&, Location&, Weathers&) const;

    bool readGeoFile(const std::filesystem::path&, std::vector<GeoConf>&, 
        int&, NameIds&, ShadowSet&, WWRs&, double, double, bool) const;

    bool readGeoData(const std::filesystem::path&, const std::filesystem::path&,
        std::vector<GeoConf>&, NameIds&, ShadowSet&, WWRs&, double, double, double, bool) const;
    
    bool confs2Phds(std::vector<GeoConf>&, GPhds&, OneMany&) const;

    bool readMaterials(const std::filesystem::path&, const NameIds&, BMaterials&) const;

    int meshSrf(SrfMesh&, const Surface&, const ep3&, const ep3&, double, double) const;
    void meshType(SrfMesh&, std::string_view typeName) const;
    bool readMeshs_(std::vector<SrfMesh>&, const std::vector<std::string_view>&,
         const GPhds&, const std::vector<int>&, double, double) const;
    bool _readMeshs(std::vector<SrfMesh>&, const std::vector<std::string_view>&, const GPhds&) const;
    bool readMeshs(SrfMeshs&, const std::filesystem::path&, const GPhds&, 
        const OneMany&, const NameIds&, double, double) const;

    void writeMeshs(const SrfMeshs&, const std::filesystem::path&, const BFlags&) const;
};
static const Config g_config;
#endif