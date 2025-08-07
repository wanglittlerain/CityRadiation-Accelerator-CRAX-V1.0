#include "env_.h"
#include "config.h"
void init(GPhds& phds, const BFlags& bflags, const WWRs& wwrs) {
    for (auto& phd : phds) {
        if (!bflags.contains(phd._buildId)) continue;
        auto it = wwrs.find(phd._buildId);
        for (auto& srf : phd.srfs) {
            if (it != wwrs.end()) {
                it->second.set(srf.window, srf._dir);
                srf.wpoly(phd.floor);
            }
        }
    }
}

int main(int argc, char** argv) {
    Elapsed rt;
    if (argc < 2) {
        std::cout << "CRAX mesh req err" << std::endl;
        return 1;
    }
    std::string req = argv[1];
    if (req.empty()) {
        std::cout << "CRAX mesh req empty" << std::endl;
        return 1;
    }
    std::filesystem::path file{req};
    if (!std::filesystem::exists(file)) {
        std::cout << "CRAX mesh req err " << req << std::endl;
        return 1;
    }

    std::ifstream ifs(file);
    json jconf = json::parse(ifs);
    auto jval = [&jconf]<typename T>(std::string_view key, T&& def) {
        T val{def};
        if (jconf.contains(key)) {
            val = jconf[key];
        }
        return val;
    };

    bool sensor = jconf["calculate-sensor-data"];
    if (!sensor) {
        std::cout << "CRAX mesh err calculate-sensor-data false" << std::endl;
        return 1;
    }
    bool input = jconf["using-cea-sensor"];
    if (input) {
        std::cout << "CRAX mesh err using-cea-sensor true" << std::endl;
        return 1;
    }

    auto fileDir{file.parent_path()};
    std::string sdir = jconf["CRAX_input_folder"];
    fileDir = sdir;

    NameIds nameIds;
    const json& jbs = jconf["buildings"];
    for (auto& b : jbs) {
        nameIds.emplace(std::string(b), 0);
    }

    std::string zoneGeo = jval("zone_geo", std::string("zone_building_geometry.csv"));
    auto geo1{fileDir / zoneGeo};
    double sdis1 = jval("zone-geometry", 0.0);
    std::string surrGeo = jval("surrounding-geo", std::string("surroundings_building_geometry.csv"));
    double sdis2 = jval("surrounding-geometry", 0.0);
    auto geo2{fileDir / surrGeo};
    bool neglect = jconf["neglect-adjacent-buildings"];
    GPhds phds;
    ShadowSet cset;
    BFlags bflags;
    OneMany b2pids;
    {
        std::vector<GeoConf> confs;
        WWRs wwrs;
        if (!g_config.readGeoData(geo1, geo2, confs, nameIds, cset, wwrs, sdis1, sdis2, 0.0, neglect)) {
            std::cout << "readGeoData err" << std::endl;
            return 1;
        }
        g_config.confs2Phds(confs, phds, b2pids);
        for (auto id : cset.ids) {
            const auto& conf = confs[id];
            auto& info = bflags[id];
            info.name = conf.name;
            info.te = conf.te;
        }
        init(phds, bflags, wwrs);
    }

    SrfMeshs meshs;
    std::string output = jconf["CRAX_result_folder"];
    auto droof{jval("roof-grid", 1.0)};
    auto dwall{jval("walls-grid", 1.0)};
    base_generate_meshs(meshs, phds, b2pids, cset, droof, dwall);

    fileDir = output;
    g_config.writeMeshs(meshs, fileDir, bflags);
    rt.view("running time: ");
    return 0;
}
