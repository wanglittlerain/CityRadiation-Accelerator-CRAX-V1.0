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
using BAs = std::vector<BuildingAreas>;

BFlags bflags;
std::map<int, std::vector<std::array<double, 10>>> brs; //buildings hourly radiation
std::map<int, std::vector<std::vector<double>>> mris; //buildings hourly meshs radiation intensity

struct NbhdReflectSrfs { //neighborhood_reflect_radiation
    std::map<int, std::tuple<double, double>> nfs; //uid,h,area
    double tarea{};
};

struct NbhdReflectBuilding {
    std::set<int> nfs;
    double tarea{};
};
std::map<int, NbhdReflectSrfs> srf_nrfs; //wall neighborhood
std::map<int, NbhdReflectBuilding> building_nrfs; //roof neighborhood

struct NReflect {
    std::map<int, double> rads;

    inline void radiation(const Surface& srf, int uid) {
        rads[uid] = srf._amwall * srf._sr + srf._amwin * srf._wsr;
    }

    inline double avg_intensity(int uid, bool roof) const {
        auto tnrr{0.0};
        auto tarea{0.0};
        auto add = [&](int nid) {
            auto nit = rads.find(nid);
            if (nit != rads.end()) {
                tnrr += nit->second;
            }
        };
        auto rfun = [&](const auto& nrfs) {
            const auto it = nrfs.find(uid);
            if (it != nrfs.end()) {
                tarea = it->second.tarea;
                for (auto& nf : it->second.nfs) {
                    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(nf)>, int>) {
                        add(nf);
                    } else {
                        add(nf.first);
                    }
                }
                if (tarea > c_1e_8) return tnrr / tarea;
            }
            return 0.0;
        };

        if (roof) return rfun(building_nrfs);
        return rfun(srf_nrfs);
    }
};

static constexpr int uniqueKey{1000};
inline int _uid(int pid, int sid) {
    return pid * uniqueKey + sid;
}

void srf_neighborhoods(const Surface& srf, int pid, double te, double tansh, const GPhds& phds) {
    if (srf.id == c_roof_Id) return;
    auto& fns = srf_nrfs[_uid(pid, srf.id)];
    auto h{srf._h - te};
    for (auto& phd : phds) {
        if (phd._id == pid) continue;
        for (auto& nsrf : phd.srfs) {
            auto relative{nsrf._h - te};
            if (nsrf.id == c_roof_Id && h <= relative) continue;
            ep3 dir_v{srf._center - nsrf._center};
            if (srf._normi.dot(dir_v) >= 0.0) continue;
            dir_v(2) = 0.0;
            if (relative < dir_v.norm() * tansh) continue;
            if (srf._normi.dot(nsrf._normi) >= 0.0) continue;
            auto nuid{_uid(phd._id, nsrf.id)};
            auto it = fns.nfs.emplace(nuid, std::make_tuple(relative, nsrf._area));
            if (it.second) {
                fns.tarea += nsrf._area;
            }
        }
    }
}

void roof_neighborhoods(const GPhds& phds, const OneMany& b2pids) {
    for (const auto& [bid, pids] : b2pids) {
        auto& bns = building_nrfs[bid];
        for (auto pid : pids) {
            auto& phd = phds.at(pid);
            auto bh{phd.box(5) - phd.box(2)};
            for (auto& srf : phd.srfs) {
                if (srf.id <= c_roof_Id) continue;
                const auto nit = srf_nrfs.find(_uid(pid, srf.id));
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

void init_srfs(GPhds& phds, const WWRs& wwrs, const BMaterials& bms, 
    const OneMany& b2pids, double phi) {
    phi = deg2rad(phi);
    auto tansh{std::tan(phi)};
    for (auto& phd : phds) {
        auto mit = bms.find(phd._buildId);
        if (mit == bms.end()) {
            mit = bms.begin();
        }
        auto it = wwrs.find(phd._buildId);
        for (auto& srf : phd.srfs) {
            if (it != wwrs.end()) {
                it->second.set(srf.window, srf._dir);
                srf.wpoly(phd.floor);
            }
            if (srf.id == c_roof_Id) {
                srf._amwall = srf._opaque * mit->second.r_roof;
            } else {
                srf._amwall = srf._opaque * mit->second.r_wall;
            }
            srf._amwin = srf.window.area * mit->second.G_win;
            srf_neighborhoods(srf, phd._id, phd.box(2), tansh, phds);
        }
    }
    roof_neighborhoods(phds, b2pids); 
}

inline void correct_cosi(double& cosi, double sh) {
    if (sh <= 0.0 || cosi >= 0.0) {
        cosi = 0.0;
    } else {
        cosi *= -1.0;
    }
}

inline double radi(double sh, const StepCoef& sc, double cosi, const DiffuseR& dr, double sr) {
    return direct_ground_diffuse(sh, sc, cosi, dr.rdome, dr.rhorizon, dr.rdomeG, sr);
}

inline void baseIntensity(const Surface& srf, NReflect& nr, int phdId, double sh, const StepCoef& sc) {
    auto& cosi{srf._cos_s};
    correct_cosi(cosi, sh);
    srf._sr = radi(sh, sc, cosi, srf._dr, srf._sr);
    auto wtr{0.0};
    if (srf.window.exist()) {
        wtr = radi(sh, sc, cosi, srf._dr, srf._wsr);
    }
    srf._wsr = wtr;
    nr.radiation(srf, _uid(phdId, srf.id));
}

inline void arrayAdd(auto& arr, int id, int dir, double v1, double v2) {
    if (id == c_roof_Id) {
        arr[8] += v1;
        arr[9] += v2;
    } else {
        arr[dir] += v1;
        arr[dir + 4] += v2;
    }
}

void phdsRadiation(const GPhds& phds, const NReflect& nr, int step) {
    int buildId{-1};
    auto bnri{0.0};
    for (const auto& phd : phds) {
        if (!bflags.contains(phd._buildId)) continue;
        if (phd._buildId != buildId) {
            buildId = phd._buildId;
            bnri = nr.avg_intensity(buildId, true);
        }
        for (auto& srf : phd.srfs) {
            auto nri{bnri};
            if (srf.id != c_roof_Id) {
                nri = nr.avg_intensity(_uid(phd._id, srf.id), false);
            }
            auto nrr{nri * srf._dr.factor()};
            auto& r = brs.at(buildId).at(step);
            auto t{srf._opaque * (srf._sr + nrr) * 1e-3};
            auto wt{srf.window.area * (srf._wsr + nrr) * 1e-3};
            arrayAdd(r, srf.id, srf._dir, t, wt);
        }
    }
}

void radiations(const GPhds& phds, SrfMeshs& smeshs, NReflect& nr, 
    const Solar& solar, const StepCoef& sc, bool save) {
    for (const auto& phd : phds) {
        for (auto& srf : phd.srfs) {
            baseIntensity(srf, nr, phd._id, solar.h, sc);
        }
    }
    if (smeshs.empty()) {
        phdsRadiation(phds, nr, solar._step);
        return;
    }
    auto pos{(solar._step - solar._start) % 24};
    auto clear{solar._step == solar._start};
    for (auto& [buildId, meshs] : smeshs) {
        auto bnri{nr.avg_intensity(buildId, true)};
        auto last_uid{-1};
        auto last_nri{0.0};
        auto& mrs = mris.at(buildId).at(solar._step);
        mrs.reserve(meshs.size());
        for (auto& mesh : meshs) {
            mesh.shadow(pos, clear, save);
            auto t{radi(solar.h, sc, mesh.srf->_cos_s, mesh._dr, mesh.m.sr)};
            auto nri{bnri};
            if (mesh.srf->id != c_roof_Id) {
                auto uid{_uid(mesh.phd, mesh.srf->id)};
                if (last_uid != uid) {
                    last_uid = uid;
                    last_nri = nr.avg_intensity(uid, false);
                }
                nri = last_nri;
            }
            t += nri * mesh._dr.factor();
            mrs.emplace_back(t);
        }
    } 
}

inline double td_(double val) {
    return std::round(val * 1e3) * 1e-3;
}

inline float tf_(double val) {
    return static_cast<float>(std::round(val * 1e2) * 1e-2);
}

void buildingArea(BAs& bas, const GPhds& phds) {
    bas.reserve(bflags.size());
    BuildingAreas ba;
    ba.id = -1;
    auto add = [&]() {
        if (ba.id < 0) return;
        for (int i = 0; i < static_cast<int>(ba.as.size()); ++i) {
            ba.as[i] = td_(ba.as[i]);
        }
        bas.emplace_back(std::move(ba));
        ba.as = {};
    };

    for (auto& phd : phds) {
        auto it = bflags.find(phd._buildId);
        if (it != bflags.end()) {
            if (ba.id != phd._buildId) {
                add();
                ba.id = phd._buildId;
                ba.name.swap(it->second.name);
            }
            for (auto& srf : phd.srfs) {
                arrayAdd(ba.as, srf.id, srf._dir, srf._opaque, srf.window.area);
            }
        }
    }
    add();
}

inline std::string R_line(std::string_view date, const auto& rs, const auto& as) {
    return std::format(R_row, date, td_(rs[8]), td_(rs[0]), td_(rs[1]), 
        td_(rs[2]), td_(rs[3]), td_(rs[4]), td_(rs[5]), td_(rs[6]), td_(rs[7]), 
        as[8], as[0], as[1], as[2], as[3], as[4], as[5], as[6], as[7]);
}

void write(const std::filesystem::path& folder, const GPhds& phds, 
    const SrfMeshs& meshs, const Weathers& ws, bool mt) {
    std::filesystem::create_directory(folder);
    BAs bas;
    buildingArea(bas, phds);
    auto f1 = [&](int i, int begin, int end) {
        for (i = begin; i < end; ++i) {
            auto& ba = bas[i];
            auto file{folder / std::format(R_name, ba.name)};
            std::ofstream ofs(file);
            ofs << R_title;
            const auto& r = brs[ba.id];
            for (auto idx = 0; idx < static_cast<int>(r.size()); ++idx) {
                ofs << R_line(ws[idx].date, r[idx], ba.as);
            }
            ofs << std::endl;
        }
    };

    auto f2 = [&](int k, int begin, int end) {
        for (k = begin; k < end; ++k) {
            auto& ba = bas[k];
            auto it = meshs.find(ba.id);
            if (it == meshs.end() || it->second.empty()) continue;
            auto file{folder / std::format(R_name, ba.name)};
            std::ofstream ofs(file);
            ofs << R_title;
            std::vector<arrow::FloatBuilder> builders(it->second.size());
            const auto& rs = mris[ba.id];
            for (int i = 0; i < static_cast<int>(rs.size()); ++i) {
                std::array<double, 10> a{};
                auto& rsu = rs[i];
                for (int j = 0; j < static_cast<int>(rsu.size()); ++j) {
                    auto rsuj = rsu[j];
                    builders[j].Append((tf_(rsuj))).ok();
                    auto& mesh = it->second[j];
                    a[mesh.dir] += rsuj * mesh.m.area;
                }
                ofs << R_line(ws[i].date, a, ba.as);
            }
            ofs << std::endl;

            auto writeF = [&]() {
                std::vector<std::shared_ptr<arrow::Array>> arrays(builders.size());
                for (int i = 0; i < static_cast<int>(builders.size()); ++i) {
                    builders[i].Finish(&arrays[i]).ok();
                }
                std::vector<std::shared_ptr<arrow::Field>> fields;
                fields.reserve(it->second.size());
                for (int i = 0; i < static_cast<int>(it->second.size()); ++i) {
                    fields.emplace_back(arrow::field(std::format(c_srf_name, i), arrow::float32()));
                }
                auto schema = arrow::schema(fields);
                auto table = arrow::Table::Make(schema, arrays);
                std::shared_ptr<arrow::io::FileOutputStream> cfs;
                static constexpr std::string_view mname{"{}_insolation_Whm2.feather"};
                auto mfile{folder / std::format(mname, ba.name)};
                ARROW_ASSIGN_OR_RAISE(cfs, arrow::io::FileOutputStream::Open(mfile.string(), false));
                arrow::ipc::feather::WriteProperties prop;
                prop.Defaults();
                prop.compression = arrow::Compression::ZSTD;
                ARROW_RETURN_NOT_OK(arrow::ipc::feather::WriteTable(*table, cfs.get(), prop));
                ARROW_RETURN_NOT_OK(cfs->Close());
                return arrow::Status::OK();
            };
            writeF().ok();
        }
    };

    auto fun = [](const auto& f, size_t size, bool mt) {
        if (mt) {
            asyncF(f, size);
            return;
        }
        f(0, 0, static_cast<int>(size));
    };

    if (meshs.empty()) {
        fun(f1, bas.size(), mt);
        return;
    }
    fun(f2, bas.size(), mt);
}

inline void init_meshs(SrfMeshs& tms, const GPhds& tps) {
    auto phd{-1};
    auto id{-1};
    const Surface* msrf{};
    for (auto& [_, meshs] : tms) {
        for (auto& mesh : meshs) {
            if (phd != mesh.phd || id != mesh.srf->id) {
                phd = mesh.phd;
                id = mesh.srf->id;
                msrf = &(tps.at(phd).srfs.at(id));
            }
            mesh.srf = msrf;
        }
    }
}

int main(int argc, char** argv) {
    Elapsed rt;
    auto err = [](std::string_view msg) {
        std::cout << msg;
        return 1;
    };

    if (argc < 2) return err("argc < 2");
    std::filesystem::path file{argv[1]};
    if (!std::filesystem::exists(file)) {
        return err("no json");
    }
    Location _lo;
    Weathers _ws;
    ShadowSet cset;
    GPhds phds;
    SrfMeshs meshs;
    auto folder{file.parent_path()};
    auto albedo{0.0};
    size_t size{};
    {
        std::ifstream ifs(file);
        json jconf = json::parse(ifs);
        auto jval = [&jconf](std::string_view key, auto&& def) {
            auto val{def};
            if (jconf.contains(key)) {
                val = jconf[key];
            }
            return val;
        };

        albedo = jconf["albedo"];
        std::string sdir = jconf["CRAX_input_folder"];
        folder = sdir;
        auto wconf{jval("weather", std::string("weather.epw"))};
        auto wfile{folder / wconf};
        if (!g_config.readWeather(wfile, _lo, _ws)) {
            return err("weather err");
        }
        size = _ws.size();
        cset.mt = jval("multiprocessing", true);
        cset.separate = jval("mesh-clipping-view-factor", false);
        cset.day = jconf["update-shadow-day"];
        NameIds names;
        const json& jbs = jconf["buildings"];
        for (auto& b : jbs) {
            names.emplace(std::string(b), 0);
        }

        auto zoneGeo{jval("zone-geo", std::string("zone_building_geometry.csv"))};
        auto geo1{folder / zoneGeo};
        auto sdis1{jval("zone-geometry", 0.0)};
        auto surrGeo{jval("surrounding-geo", std::string("surroundings_building_geometry.csv"))};
        auto sdis2{jval("surrounding-geometry", 0.0)};
        auto geo2{folder / surrGeo};
        bool neglect = jconf["neglect-adjacent-buildings"];
        auto te{jval("te", 0.0)};
        std::vector<GeoConf> confs;
        WWRs wwrs;
        OneMany b2pids;
        if (!g_config.readGeoData(geo1, geo2, confs, names, cset, wwrs, sdis1, sdis2, te, neglect)) {
            return err("geo err");
        }
        g_config.confs2Phds(confs, phds, b2pids);
        auto mat{jval("materials", std::string("building_materials.csv"))};
        auto mfile{folder / mat};
        BMaterials bms;
        if (!g_config.readMaterials(mfile, names, bms)) {
            return err("materials err");
        }
        auto sh{jval("neighbor-filter-altitude-angle", 15.0)};
        init_srfs(phds, wwrs, bms, b2pids, sh);
    
        bool sensor = jconf["calculate-sensor-data"];
        for (auto id : cset.ids) {
            auto& conf = confs[id];
            auto& info = bflags[id];
            info.name.swap(conf.name);
            info.te = conf.te;
            if (!sensor) {
                brs[id].resize(size);
            }
        }
        sdir = jconf["CRAX_result_folder"];
        folder = sdir;
        if (sensor) {
            bool input = jconf["using-cea-sensor"];
            if (input) {
                auto pnorm{jval("pnorm", 1e-2)};
                auto pinner{jval("pinner", 1e-3)};
                g_config.readMeshs(meshs, folder, phds, b2pids, names, pnorm, pinner);
            } else {
                auto droof{jval("roof-grid", 1.0)};
                auto dwall{jval("walls-grid", 1.0)};
                base_generate_meshs(meshs, phds, b2pids, cset, droof, dwall);
                g_config.writeMeshs(meshs, folder, bflags);
            }
            for (const auto& [id, _] : meshs) {
                auto& rs = mris[id];
                rs.resize(size);
            }
            if (meshs.empty()) return err("sensor empty");
        }
    }
    {
        auto claculate = [&](GPhds& tps, SrfMeshs& tms, int begin, int end) {
            Solar solar(_lo.latitude, _lo.longitude, _lo.ls, begin);
            NReflect nr;
            StepCoef coef{};
            for (int step = begin; step < end; ++step) {
                solar.sv(step);
                g_shadow.calculate(tps, solar, cset);
                const auto& w = _ws[step];
                coef.dhr = w.dhr;
                coef.dnr = w.dnr;
                dr_coef(coef, solar.h, -solar.v(2), albedo);
                radiations(tps, tms, nr, solar, coef, cset.day > 1);
            }
        };

        cset.dset(true, !meshs.empty());
        if (cset.mt && _ghc > 1) {
            auto hc{static_cast<int>(_ghc) - 1};
            std::vector<GPhds> tphds(hc, phds);
            std::vector<SrfMeshs> tmeshs(hc, meshs);;
            auto df = [&](int i, int begin, int end) {
                if (i < hc) {
                    init_meshs(tmeshs[i], tphds[i]);
                    g_shadow.diffuse(tphds[i], tmeshs[i], cset, begin, end);
                    return;
                }
                g_shadow.diffuse(phds, meshs, cset, begin, end);
            };
            asyncF(df, c_diffuse_step);

            auto merge = [&]() {
                std::vector<DiffuseR*> drs(hc);
                if (!cset.separate) {
                    for (int pi = 0; pi < static_cast<int>(phds.size()); ++pi) {
                        for (auto& srf : phds[pi].srfs) {
                            for (int i = 0; i < hc; ++i) {
                                auto& idr = tphds[i][pi].srfs[srf.id]._dr;
                                drs[i] = &idr;
                                srf._dr += idr;
                            }
                            srf._dr.calc(srf._sin_g, srf._normi(2));
                            for (int i = 0; i < hc; ++i) {
                                drs[i]->operator=(srf._dr);
                            }
                        }
                    }
                }
                for (auto& [buildId, ms] : meshs) {
                    for (int mi = 0; mi < static_cast<int>(ms.size()); ++mi) {
                        auto& mesh = ms[mi];
                        for (int i = 0; i < hc; ++i) {
                            auto& idr = tmeshs[i][buildId][mi]._dr;
                            if (cset.separate) {
                                mesh._dr += idr;
                                drs[i] = &idr;
                            } else {
                                idr = mesh.srf->_dr;
                            }
                        }
                        if (cset.separate) {
                            mesh._dr.calc(mesh.srf->_sin_g, mesh.srf->_normi(2));
                            for (int i = 0; i < hc; ++i) {
                                drs[i]->operator=(mesh._dr);
                            }
                        } else {
                            mesh._dr = mesh.srf->_dr;
                        }
                    }
                }
            };
            merge();
            cset.dset();
            auto fun = [&](int i, int begin, int end) {
                if (i < hc) {
                    claculate(tphds[i], tmeshs[i], begin, end);
                } else {
                    claculate(phds, meshs, begin, end);
                }
            };
            asyncF(fun, size);
        } else {
            g_shadow.diffuse(phds, meshs, cset);
            cset.dset();
            claculate(phds, meshs, 0, static_cast<int>(size));
        }
        rt.view("model running time: ");
    }
    write(folder, phds, meshs, _ws, cset.mt);
    rt.view("file writing time: ");
    return 0;
}
