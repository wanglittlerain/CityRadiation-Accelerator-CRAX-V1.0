#include "env_.h"
#include "config.h"
void initPhdSrfs(GPhds& phds, const BFlags& bflags, const WWRs& wwrs) {
    for (auto& phd : phds) {
        auto pit = bflags.find(phd._buildId);
        if (pit == bflags.end()) continue;
        auto it = wwrs.find(phd._buildId);
        for (auto& [_, srf] : phd.srfs) {
            if (srf.id <= c_virtual_srf) continue;
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
    auto jval_fun = [&jconf](auto& val, std::string_view key) {
        if (jconf.contains(key)) {
            val = jconf[key];
        }
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

    ShadowSet cset;
    const json& jbs = jconf["buildings"];
    for (auto& b : jbs) {
        cset.nameIds.emplace(std::string(b), 0);
    }

    std::string zoneGeo{"zone_building_geometry.csv"};
    jval_fun(zoneGeo, "zone-geo");
    auto geo1{fileDir / zoneGeo};
    double sdis1{};
    jval_fun(sdis1, "zone-geometry");
    std::string surrGeo{"surroundings_building_geometry.csv"};
    jval_fun(surrGeo, "surrounding-geo");
    double sdis2{};
    jval_fun(sdis2, "surrounding-geometry");
    auto geo2{fileDir / surrGeo};
    bool neglect = jconf["neglect-adjacent-buildings"];
    GPhds phds;
    WWRs wwrs;
    BFlags bflags;
    OneMany b2pids;
    {
        std::vector<GeoConf> confs;
        if (!sRC.readGeoData(geo1, geo2, confs, cset, wwrs, sdis1, sdis2, neglect)) {
            std::cout << "readGeoData err" << std::endl;
            return 1;
        }
        sRC.confs2Phds(confs, phds, b2pids);
        for (auto id : cset.ids) {
            const auto& conf = confs[id];
            auto& info = bflags[id];
            info.name = conf.name;
            info.te = conf.te;
        }
    }

    initPhdSrfs(phds, bflags, wwrs);

    SrfMeshs meshs;
    std::string output = jconf["CRAX_result_folder"];
    auto droof{1.0};
    jval_fun(droof, "roof-grid");
    auto dwall{1.0};
    jval_fun(dwall, "walls-grid");
    base_generateSrfMeshs(meshs, phds, b2pids, cset, droof, dwall);

    fileDir = output;
    sRC.writeMeshs(meshs, fileDir, bflags);
    rt.view("running time: ");
    return 0;
}
