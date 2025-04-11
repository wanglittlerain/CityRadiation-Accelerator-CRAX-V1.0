#pragma once
#ifndef _CONFIG_H_
#define _CONFIG_H_
#include "base.h"
static constexpr std::array c_orientations{"east", "north", "south", "west", "top"};
static constexpr std::array c_srfType{"walls", "windows", "roofs"};
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
    bg_polygon poly;
    std::vector<e_seg> segs;
    std::vector<int> ins;

    inline void pushSegs() {
        for (int i = 0; i < static_cast<int>(poly.outer().size()); ++i) {
            bg_point pre{};
            if (i == 0) {
                pre = poly.outer().back();
            } else {
                pre = poly.outer()[i - 1];
            }
            auto& cur = poly.outer()[i];
            e_seg seg{ep2{pre.x(), pre.y()}, ep2{cur.x(), cur.y()}};
            segs.emplace_back(seg);
        }  
    }
};
using WWRs = std::map<int, WWRatio>;

struct Intersects {
    int i1{};
    int i2{};
    int p1{};
    int p2{};
    bool integral{};
    double minH{};
    e_seg seg;
};

struct Materials{
    double G_win{};
    double r_roof{};
    double r_wall{};
};

struct BuildingFlag {
    std::string name{};
    double te{};
};
using BFlags = std::map<int, BuildingFlag>;
using BMaterials = std::map<int, Materials>;

class ReadConfig {
public:
    bool readWeather(const std::filesystem::path&, Location&, Weathers&) const;

    bool readGeoFile(const std::filesystem::path&, std::vector<GeoConf>&, 
        int&, ShadowSet&, WWRs&, double, bool) const;

    bool readGeoData(const std::filesystem::path&, const std::filesystem::path&,
        std::vector<GeoConf>&, ShadowSet&, WWRs&, double, double, bool) const;
    
    bool confs2Phds(std::vector<GeoConf>&, GPhds&, OneMany&) const;

    bool readMaterials(const std::filesystem::path&, const NameIds&, BMaterials&) const;

    int meshSrf(SrfMesh&, const Surface&, const ep3&, const ep3&, double, double) const;
    bool readMeshs_(std::vector<SrfMesh>&, const std::vector<std::string>&,
         const GPhds&, const std::vector<int>&, double, double) const;
    bool _readMeshs(std::vector<SrfMesh>&, const std::vector<std::string>&, const GPhds&) const;
    bool readMeshs(SrfMeshs&, const std::filesystem::path&, const GPhds&, 
        const OneMany&, const ShadowSet&, double, double) const;

    void writeMeshs(const SrfMeshs&, const std::filesystem::path&, const BFlags&) const;
};
static const ReadConfig sRC;
#endif