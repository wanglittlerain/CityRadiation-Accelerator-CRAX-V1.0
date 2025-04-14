#include "tp.h"
#include "rad_func.h"
#include "arrow_env.h"
#include "shadow.h"
#include "config.h"
static constexpr std::string_view R_title{"Date,roofs_top_kW,walls_east_kW,walls_north_kW,walls_south_kW,walls_west_kW,windows_east_kW,windows_north_kW,windows_south_kW,windows_west_kW,roofs_top_m2,walls_east_m2,walls_north_m2,walls_south_m2,walls_west_m2,windows_east_m2,windows_north_m2,windows_south_m2,windows_west_m2\n"};
static constexpr std::string_view R_row{"{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"};
static constexpr std::string_view R_name{"{}_radiation.csv"};

struct BuildingAreas {
    int id{};
    std::string name{};
    std::array<double, 10> as{};
};
using BAreas = std::vector<BuildingAreas>;

BFlags bflags;
BMaterials bms;
using Radiations = std::map<int, std::vector<std::array<double, 10>>>;
Radiations radRes;

using MeshRadiations = std::map<int, std::vector<std::vector<double>>>;
MeshRadiations meshRads;
std::vector<int> meshIds;

struct NeighborReflectSrfs { //neighbor_reflect_radiation
    std::map<int, std::tuple<double, double>> nfs; //uid,h,area
    double tarea{};
};
struct NeighborReflectBuilding {
    std::set<int> nfs;
    double tarea{};
};
std::map<int, NeighborReflectSrfs> srf_nrfs; //wall neighbour
std::map<int, NeighborReflectBuilding> building_nrfs; //roof neighbour
using KdVs = std::map<int, double>;

static constexpr int uniqueKey{10000};
int uniqueId(int pid, int sid) {
    return pid * uniqueKey + sid;
}

void phd_srf(int uid, int& pid, int& sid) {
    pid = uid / uniqueKey;
    sid = uid % uniqueKey;
}

inline void srf_neighbors(const Surface& srf, int pid, double tansh, const GPhds& phds) {
    if (srf.id == c_roof_Id) return;
    auto uid{uniqueId(pid, srf.id)};
    auto& fns = srf_nrfs[uid];
    for (auto& phd : phds) {
        if (phd._id == pid) continue;
        for (auto& [_1, nsrf] : phd.srfs) {
            if (nsrf.id <= c_virtual_srf) continue;
            if (nsrf.id == c_roof_Id && srf._h <= nsrf._h) continue;
            auto threshold{nsrf._h / tansh};
            ep3 dir_v{srf._center - nsrf._center};
            auto dis{dir_v.norm()};
            if (dis > threshold) continue;

            auto ndot{srf._normi.dot(nsrf._normi)};
            if (ndot >= 0.0) continue;

            auto ddot{srf._normi.dot(dir_v)};
            if (ddot >= 0.0) continue;

            auto nuid{uniqueId(phd._id, nsrf.id)};
            auto it = fns.nfs.emplace(nuid, std::make_tuple(nsrf._h, nsrf._area));
            if (it.second) {
                fns.tarea += nsrf._area;
            }
        }
    }
}

inline void roof_neighbors(const GPhds& phds, const OneMany& b2pids) {
    for (const auto& [bid, pids] : b2pids) {
        auto& bns = building_nrfs[bid];
        for (auto pid : pids) {
            auto& phd = phds.at(pid);
            auto bh{phd.box(5)};
            for (auto& [_, srf] : phd.srfs) {
                if (srf.id <= c_roof_Id) continue;
                auto uid{uniqueId(pid, srf.id)};
                const auto nit = srf_nrfs.find(uid);
                if (nit == srf_nrfs.end()) continue;
                for (auto& [nid, val] : nit->second.nfs) {
                    auto[nh, narea] = val;
                    if (bh >= nh) continue;
                    auto it = bns.nfs.emplace(nid);
                    if (it.second) {
                        bns.tarea += narea;
                    }
                }
            }
        }        
    }
}

void initPhdSrfs(GPhds& phds, const WWRs& wwrs, const OneMany& b2pids, double phi) {
    phi = deg2rad(phi);
    auto tansh{std::tan(phi)};
    for (auto& phd : phds) {
        auto it = wwrs.find(phd._buildId);
        for (auto& [_, srf] : phd.srfs) {
            if (srf.id <= c_virtual_srf) continue;
            if (it != wwrs.end()) {
                it->second.set(srf.window, srf._dir);
                srf.wpoly(phd.floor);
            }
            srf_neighbors(srf, phd._id, tansh, phds);
        }
    }
    roof_neighbors(phds, b2pids);
}

inline void refRadiation(const Surface& srf, int uid, const Materials& m, KdVs& nrrs) {
    auto nrr{0.0};
    if (srf.id == c_roof_Id) {
        nrr = srf._area * srf._tr * m.r_roof;
    } else {
        nrr = (srf._opaque * srf._tr * m.r_wall + srf.window.area * srf._wtr * (1 - m.G_win));
    }
    nrrs[uid] = nrr;
}

inline double nrefRadiation(int uid, const KdVs& nrrs, bool roof) {
    auto tnrr{0.0};
    auto tarea{0.0};
    auto nrr_fun = [&](int nid) {
        auto nit = nrrs.find(nid);
        if (nit == nrrs.end()) return;
        tnrr += nit->second;
    };

    if (roof) {
        const auto it = building_nrfs.find(uid);
        if (it != building_nrfs.end()) {
            tarea = it->second.tarea;
            for (auto nid : it->second.nfs) {
                nrr_fun(nid);
            }
        }
    } else {
        const auto it = srf_nrfs.find(uid);
        if (it == srf_nrfs.end()) return tnrr;
        tarea = it->second.tarea;
        for (auto& [nid, _] : it->second.nfs) {
            nrr_fun(nid);
        }
    }
    if (tarea > c_1e_8) {
        return tnrr / tarea;
    }
    return 0.0;
}

inline void correct_cosi(double& cosi, double sh) {
    if (sh <= 0.0 || cosi >= 0.0) {
        cosi = 0.0;
    } else {
        cosi *= -1.0;
    }
}

bool srfRadiation(const Surface& srf, KdVs& nrrs, int phdId, const Solar& solar, 
    const StepCoef& sc, const Materials& m) {
    if (srf.id <= c_virtual_srf) return false;
    if (srf._dir < 0) return false;
    auto cosi{solar.v.dot(srf._normi)};
    correct_cosi(cosi, solar.h);
    auto wall_dr{cosi * sc.dnr * (1 - srf._sr)};
    auto dr{diffuse_(solar.h, sc, cosi, srf._gdot, srf._sinsita, 
        srf._dr.rdome, srf._dr.rhorizon, srf._sr)};
    auto grr{ground_reflect_(solar.h, sc.c4, srf._gdot, srf._dr.rdomeG)};
    srf._tr = wall_dr + dr + grr;

    srf._wtr = 0.0;
    if (srf.window.exist()) {
        auto window_dr{cosi * sc.dnr * (1 - srf._wsr)};
        auto dr1{diffuse_(solar.h, sc, cosi, srf._gdot, srf._sinsita, 
            srf._dr.rdome, srf._dr.rhorizon, srf._wsr)};
        srf._wtr = window_dr + dr1 + grr;
    }
    auto uid{uniqueId(phdId, srf.id)};
    refRadiation(srf, uid, m, nrrs);
    return true;
}

bool srfAllRadiation(const Surface& srf, const KdVs& nrrs, int phdId, int buildId, int step, double avgrr) {
    if (srf.id <= c_virtual_srf) return false;
    if (srf._dir < 0) return false;
    auto uid{uniqueId(phdId, srf.id)};
    if (srf.id != c_roof_Id) {
        avgrr = nrefRadiation(uid, nrrs, false);
    }
    auto factor{reflection_factor(srf._dr.rdome, srf._dr.rdomeG, srf._gdot)};
    auto nrr{factor * avgrr};
    auto& rs = radRes[buildId];
    auto& r = rs[step];
    if (srf.id == c_roof_Id) {
        r[8] += srf._opaque * (srf._tr + nrr) * 1e-3;
        r[9] += srf.window.area * (srf._wtr + nrr) * 1e-3;
    } else {
        r[srf._dir] += srf._opaque * (srf._tr + nrr) * 1e-3;;
        r[srf._dir + 4] += srf.window.area * (srf._wtr + nrr) * 1e-3;
    }
    return true;
}

void phdsRadiation(const GPhds& phds, KdVs& nrrs, const Solar& solar, int step, const StepCoef& sc, bool calc) {
    for (const auto& phd : phds) {
        auto mit = bms.find(phd._buildId);
        if (mit == bms.end()) {
            mit = bms.begin();
            if (mit == bms.end()) continue;
        }
        for (auto& [_, srf] : phd.srfs) {
            srfRadiation(srf, nrrs, phd._id, solar, sc, mit->second);
        }
    }
    if (!calc) return;
    int buildId{-1};
    auto build_avgrr{0.0};
    for (const auto& phd : phds) {
        if (bflags.find(phd._buildId) == bflags.end()) continue;
        if (phd._buildId != buildId) {
            buildId = phd._buildId;
            build_avgrr = nrefRadiation(buildId, nrrs, true);
        }
        for (auto& [_, srf] : phd.srfs) {
            srfAllRadiation(srf, nrrs, phd._id, buildId, step, build_avgrr);
        }
    }
}

void meshsRadiation(const GPhds& phds, SrfMeshs& smeshs, KdVs& nrrs, 
    const Solar& solar, int step, const StepCoef& sc, bool save) {
    phdsRadiation(phds, nrrs, solar, step, sc, smeshs.empty());
    auto pos{(solar._step - solar._start) % 24};
    auto clear{solar._step == solar._start};
    for (auto& [buildId, meshs] : smeshs) {
        auto& rs = meshRads[buildId];
        for (auto& mesh : meshs) {
            mesh.shadow_(pos, clear, save);
            auto cosi{solar.v.dot(mesh.srf->_normi)};
            correct_cosi(cosi, solar.h);
            auto srf_dr{cosi * sc.dnr * (1.0 - mesh.m.sr)};
            auto dr{diffuse_(solar.h, sc, cosi, mesh.srf->_gdot, 
                mesh.srf->_sinsita, mesh._dr.rdome, mesh._dr.rhorizon, mesh.m.sr)};
            auto grr{ground_reflect_(solar.h, sc.c4, mesh.srf->_gdot, mesh._dr.rdomeG)};

            auto avgrr{0.0};
            if (mesh.srf->id == c_roof_Id) {
                avgrr = nrefRadiation(buildId, nrrs, true);
            } else {
                auto uid{uniqueId(mesh.phd, mesh.srf->id)};
                avgrr = nrefRadiation(uid, nrrs, false);
            }
            auto factor{reflection_factor(mesh._dr.rdome, mesh._dr.rdomeG, mesh.srf->_gdot)};
            auto nrr{factor * avgrr};
            auto t{srf_dr + dr + grr + nrr};
            rs[step].emplace_back(t);
        }
    } 
}

inline void td(double& val) {
    val = std::round(val * 1e3) * 1e-3;
}

inline double td_(double val) {
    val = std::round(val * 1e3) * 1e-3;
    return val;
}

inline float tf_(double val) {
    float v = static_cast<float>(std::round(val * 1e3) * 1e-3);
    return v;
}

void buildingArea(BAreas& bas, const GPhds& phds) {
    bas.reserve(bflags.size());
    std::array<double, 10> as{};
    int preId{-1};
    std::string name{};
    auto add = [&]() {
        BuildingAreas ba;
        ba.id = preId;
        for (int i = 0; i < static_cast<int>(as.size()); ++i) {
            ba.as[i] = td_(as[i]);
        }
        ba.name = name;
        bas.emplace_back(std::move(ba));
        as = {};
    };

    for (auto& phd : phds) {
        auto it = bflags.find(phd._buildId);
        if (it != bflags.end()) {
            if (preId >= 0 && preId != phd._buildId) {
                add();
            }
            for (auto& [_, srf] : phd.srfs) {
                if (srf.id <= c_virtual_srf) continue;
                if (srf._dir < 0) continue;
                if (srf.id == c_roof_Id) {
                    as[8] += srf._opaque;
                    as[9] += srf.window.area;
                } else {
                    as[srf._dir] += srf._opaque;
                    as[srf._dir + 4] += srf.window.area;
                }
            }
            if (preId != phd._buildId) {
                preId = phd._buildId;
                name = it->second.name;
            }
        }
    }
    if (preId >= 0) {
        add();
    }
}

void writePhdsRadiation(const std::filesystem::path& fileDir, const BAreas& bas, 
    const Weathers& ws, size_t hc) {
    std::filesystem::create_directory(fileDir);
    auto fun = [&](int th, int begin, int end) {
        th = begin;
        for (; th < end; ++th) {
            auto& ba = bas[th];
            auto name{std::format(R_name, ba.name)};
            auto file{fileDir / name};
            std::ofstream ofs(file);
            ofs << R_title;
            const auto& r = radRes[ba.id];
            auto idx{0};
            for (auto& a : r) {
                auto& w = ws[idx];
                auto line{std::format(R_row, w.date, td_(a[8]), td_(a[0]), td_(a[1]), td_(a[2]), 
                    td_(a[3]), td_(a[4]), td_(a[5]), td_(a[6]), td_(a[7]), ba.as[8], ba.as[0], 
                    ba.as[1], ba.as[2], ba.as[3], ba.as[4], ba.as[5], ba.as[6], ba.as[7])};
                ofs << line;
                ++idx;
            }
            ofs << std::endl;  
        }
    };
    concurrencyF(fun, bas.size(), hc);
}

void writeMeshsRadiation(const std::filesystem::path& fileDir, const SrfMeshs& meshs, 
    const BAreas& bas, const Weathers& ws, size_t hc) {
    static constexpr std::string_view mname{"{}_insolation_Whm2.feather"};
    std::filesystem::create_directory(fileDir);
    auto findBa = [&](int id) {
        int left{};
        int right{static_cast<int>(bas.size())};
        while (left <= right) {
            auto mid{(left + right) / 2};
            auto bid{bas[mid].id};
            if (bid == id) {
                return mid;
            } else if (bid < id) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        return -1;
    };


    auto fun = [&](int th, int begin, int end) {
        th = begin;
        for (; th < end; ++th) {
            auto id = meshIds[th];
            const auto it = meshs.find(id);
            if (it == meshs.end() || it->second.empty()) continue;
            auto bidx{findBa(id)};
            if (bidx < 0) continue;
            auto& ba = bas[bidx];
            auto name{std::format(R_name, ba.name)};
            auto file{fileDir / name};
            std::ofstream ofs(file);
            ofs << R_title;
            
            auto meshName{std::format(mname, ba.name)};
            auto meshFile{fileDir / meshName};
            std::vector<std::shared_ptr<arrow::Field>> fields;
            fields.reserve(it->second.size());
            for (int i = 0; i < static_cast<int>(it->second.size()); ++i) {
                fields.emplace_back(arrow::field(std::format(c_srf_name, i), arrow::float32()));
            }
            std::vector<arrow::FloatBuilder> builders;
            builders.resize(fields.size());
            const auto& rs = meshRads[id];
            std::string line{};
            for (int i = 0; i < static_cast<int>(rs.size()); ++i) {
                std::array<double, 9> a{};
                auto& rsu = rs[i];
                for (int j = 0; j < static_cast<int>(rsu.size()); ++j) {
                    auto rsuj = rsu[j];
                    builders[j].Append((tf_(rsuj))).ok();
                    auto& mesh = it->second[j];
                    a[mesh.dir] += rsuj * mesh.m.area;
                }
                for (auto& val : a) {
                    td(val);
                }
                auto& w = ws[i];
                line = std::format(R_row, w.date, a[8], a[0], a[1], a[2], a[3], a[4], 
                    a[5], a[6], a[7], ba.as[8], ba.as[0], ba.as[1], ba.as[2], ba.as[3], 
                    ba.as[4], ba.as[5], ba.as[6], ba.as[7]
                );
                ofs << line;
            }
            ofs << std::endl;
            auto writeF = [&]() {
                std::vector<std::shared_ptr<arrow::Array>> arrays;
                arrays.resize(builders.size());
                for (int i = 0; i < static_cast<int>(builders.size()); ++i) {
                    builders[i].Finish(&arrays[i]).ok();
                }
                auto schema = arrow::schema(fields);
                auto table = arrow::Table::Make(schema, arrays);
                std::shared_ptr<arrow::io::FileOutputStream> cfs;
                ARROW_ASSIGN_OR_RAISE(cfs, arrow::io::FileOutputStream::Open(meshFile.string(), false));
                arrow::ipc::feather::WriteProperties prop;
                prop.Defaults();
                prop.compression = arrow::Compression::ZSTD;
                ARROW_RETURN_NOT_OK(arrow::ipc::feather::WriteTable(*table, cfs.get(), prop));
                ARROW_RETURN_NOT_OK(cfs->Close());
                return arrow::Status::OK();
            };
            auto st = writeF();
            st.ok();
        }
    };
    concurrencyF(fun, meshIds.size(), hc);
}

inline void initThreadMeshs(SrfMeshs& tms, const GPhds& tps) {
    for (auto& [_, meshs] : tms) {
        for (auto& mesh : meshs) {
            auto& pit = tps.at(mesh.phd);
            auto fit = pit.srfs.find(mesh.srf->id);
            mesh.srf = &(fit->second);
        }
    }
}

int main(int argc, char** argv) {
    Elapsed rt;
    if (argc < 2) {
        std::cout << "CRAX req err" << std::endl;
        return 1;
    }
    std::string req = argv[1];
    if (req.empty()) {
        std::cout << "CRAX req empty" << std::endl;
        return 1;
    }
    std::filesystem::path file{req};
    if (!std::filesystem::exists(file)) {
        std::cout << "CRAX req err " << req << std::endl;
    }
    Location _lo;
    Weathers _ws;
    ShadowSet cset;
    GPhds phds;
    SrfMeshs meshs;
    auto fileDir{file.parent_path()};
    auto albedo{0.0};
    size_t size{};
    size_t hc{std::thread::hardware_concurrency()};
    {
        std::ifstream ifs(file);
        json jconf = json::parse(ifs);
        auto jval_fun = [&jconf](auto& val, std::string_view key) {
            if (jconf.contains(key)) {
                val = jconf[key];
            }
        };
        albedo = jconf["albedo"];

        std::string sdir = jconf["CRAX_input_folder"];
        fileDir = sdir;
        std::string wconf{"weather.epw"};
        jval_fun(wconf, "weather");
        auto wfile{fileDir / wconf};
        if (!sRC.readWeather(wfile, _lo, _ws)) {
            std::cout << "readWeather err " << wfile << std::endl;
            return 1;
        }
        size = _ws.size();

        cset.day = jconf["update-shadow-day"];
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
        std::vector<GeoConf> confs;
        WWRs wwrs;
        OneMany b2pids;
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
            auto& r = radRes[id];
            r.resize(size);
        }

        NameIds bnameIds;
        auto id{0};
        for (const auto& conf : confs) {
            bnameIds.emplace(conf.name, id);
            ++id;
        }
        std::string mat{"building_materials.csv"};
        jval_fun(mat, "materials");
        auto matFile{fileDir / mat};
        sRC.readMaterials(matFile, bnameIds, bms);

        auto sh{15.0};
        jval_fun(sh, "neighbor-filter-altitude-angle");
        initPhdSrfs(phds, wwrs, b2pids, sh);
    
        std::string output = jconf["CRAX_result_folder"];
        bool sensor = jconf["calculate-sensor-data"];
        if (sensor) {
            bool input = jconf["using-cea-sensor"];
            if (input) {
                auto meshDir{std::filesystem::path{output}};
                auto pnorm{c_1e_2};
                jval_fun(pnorm, "pnorm");
                auto pinner{1e-3};
                jval_fun(pinner, "pinner");
                sRC.readMeshs(meshs, meshDir, phds, b2pids, cset, pnorm, pinner);
            } else {
                auto droof{1.0};
                jval_fun(droof, "roof-grid");
                auto dwall{1.0};
                jval_fun(dwall, "walls-grid");
                base_generateSrfMeshs(meshs, phds, b2pids, cset, droof, dwall);
                fileDir = output;
                sRC.writeMeshs(meshs, fileDir, bflags);
            }
            for (const auto& [buildId, mesh] : meshs) {
                meshIds.emplace_back(buildId);
                auto& r = meshRads[buildId];
                std::vector<double> tmp;
                tmp.reserve(mesh.size());
                r.resize(size, tmp);
            }
        }
        bool alone{};
        jval_fun(alone, "mesh-clipping-view-factor");
        sShdg.diffuse(phds, meshs, cset, alone);
        fileDir = output;     
    }
    {
        std::vector<GPhds> thPhds;
        std::vector<SrfMeshs> thSmeshs;
        thPhds.resize(hc, phds);
        thSmeshs.resize(hc, meshs);
        auto fun = [&](int th, int begin, int end) {
            Solar solar;
            solar.init(_lo.latitude, _lo.longitude, _lo.ls, begin);
            auto& tps = thPhds[th];
            auto& tms = thSmeshs[th];
            initThreadMeshs(tms, tps);
            KdVs nrrs;
            StepCoef coef{};
            for (int step = begin; step < end; ++step) {
                auto day{step / 24 + 1};
                auto hour{(step + 1) % 24 - 0.5};
                solar.sv(day, hour, step);
                sShdg.calc_(solar, tps, cset);
                const auto& w = _ws[step];
                coef.dhr = w.dhr;
                coef.dnr = w.dnr;
                dr_coef(solar.h, albedo, coef);
                meshsRadiation(tps, tms, nrrs, solar, step, coef, cset.day > 1);
            }
        };
        concurrencyF(fun, size, hc);
        rt.view("model running time: ");
    }

    BAreas bas;
    buildingArea(bas, phds);
    if (meshs.empty()) {
        writePhdsRadiation(fileDir, bas, _ws, hc);
    } else {
        writeMeshsRadiation(fileDir, meshs, bas, _ws, hc);
    }
    rt.view("file writing time: ");
    return 0;
}
