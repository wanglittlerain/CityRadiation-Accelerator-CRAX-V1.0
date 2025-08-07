#include "config.h"
#include "gen_phd.h"
#include "str_func.h"
static constexpr std::array c_orientations{"east", "north", "south", "west", "top"};
static constexpr std::array c_srfType{"walls", "windows", "roofs"};

double stod_s(const std::string_view sv) {
    return sv.empty() ? 0.0 : std::stod(std::string(sv));
}

int stoi_s(const std::string_view sv) {
    return sv.empty() ? 0 : std::stoi(std::string(sv));
}

bool Config::readWeather(const std::filesystem::path& file, Location& lo, Weathers& ws) const {
    static constexpr std::string_view dateFormat{"{}/{}/{} {}:00:00"};
    if (!std::filesystem::exists(file)) return false;
    ws.reserve(8784);
    std::ifstream f(file);
    std::string l;
    auto first{true};
    while (std::getline(f, l)) {
        if (sfunc::trimr(l)) continue;
        auto&& us = sfunc::split(l, ",");
        if (first) {
            if (us.size() < 9) return false;
            lo.latitude = stod_s(us[6]);
            lo.longitude = stod_s(us[7]);
            lo.ls = stod_s(us[8]) * 15.0;
            sfunc::skip(f, 7);
            first = false;
            continue;
        }
        if (us.size() < 16) return false;
        Weather w;
        w.dnr = stod_s(us[14]);
        w.dhr = stod_s(us[15]);
        auto y{stoi_s(us[0])};
        auto m{stoi_s(us[1])};
        auto d{stoi_s(us[2])};
        auto h{stoi_s(us[3]) - 1};
        w.date = std::format(dateFormat, y, m, d, h);
        ws.emplace_back(std::move(w));
    }
    auto size = static_cast<int>(ws.size());
    return size == 8760 || size == 8784;
}

bool Config::readGeoFile(const std::filesystem::path& file, std::vector<GeoConf>& confs, int& id, 
    NameIds& nameIds, ShadowSet& cset, WWRs& wwrs, double sdis, double te, bool zoneGeo) const {
    if (!std::filesystem::exists(file)) return !zoneGeo;
    std::ifstream f(file);
    std::string l;
    std::getline(f, l);
    while (std::getline(f, l)) {
        if (sfunc::trimr(l)) continue;
        auto&& us = sfunc::split(l, "\"");
        if (us[0].empty()) continue;
        GeoConf conf{};
        conf.name = us[0];
        conf.name.pop_back();
        if (conf.name.empty()) continue;
        auto&& us2 = sfunc::split(us[2], ",");
        conf.floor = stoi_s(us2[1]);
        conf.h = stod_s(us2[2]);
        conf.te = stod_s(us2[3]) + te;
        conf.h += conf.te;
        auto it = nameIds.find(conf.name);
        if (zoneGeo && it != nameIds.end()) {
            it->second = id;
            cset.ids.emplace_back(id);
        } else {
            nameIds.emplace(conf.name, -1 - id); //when id == 0 it->second < 0
        }
        if (us.size() > 3) {
            auto&& us3 = sfunc::split(us[3], ",");
            if (us3.size() == 5) {
                WWRatio& wwr = wwrs[id];
                for (int i = 0; i < 5; ++i) {
                    wwr.rs[i] = stod_s(us3[i]);
                }
            }
        }
        auto&& us1 = sfunc::split(us[1], "), (");
        auto size = static_cast<int>(us1.size()) - 1;
        bg_polygon poly;
        poly.outer().reserve(size);
        for (int i = 0; i < size; ++i) {
            auto pos = us1[i].find(',');
            auto sx = us1[i].substr(0, pos);
            auto sy = us1[i].substr(pos + 1);
            if (i == 0) {
                sx = us1[i].substr(2, pos - 2);
            }
            poly.outer().emplace_back(bg_point({stod_s(sx), stod_s(sy)}));
        }
        if (!boost::geometry::is_valid(poly)) {
            boost::geometry::correct(poly);
        }
        if (sdis > 1e-1 && poly.outer().size() > 4) {
            bg_polygon spoly;
            boost::geometry::simplify(poly, spoly, sdis);
            poly.outer().swap(spoly.outer());
        }
        conf.pushSegs(poly);
        confs.emplace_back(std::move(conf));
        ++id;
    }
    return true; 
}

bool Config::readGeoData(const std::filesystem::path& file1, const std::filesystem::path& file2,
    std::vector<GeoConf>& confs, NameIds& nameIds, ShadowSet& cset, WWRs& wwrs, 
    double sdis1, double sdis2, double te, bool neglect) const {
    auto id{0};
    if (!readGeoFile(file1, confs, id, nameIds, cset, wwrs, sdis1, te, true)) return false;
    if (!neglect) {
        if (!readGeoFile(file2, confs, id, nameIds, cset, wwrs, sdis2, te, false)) return false;
    }
    return true;    
}

bool Config::confs2Phds(std::vector<GeoConf>& confs, GPhds& phds, OneMany& b2pids) const {
    std::vector<Insts> insts(confs.size());
    auto integral = [&](int i, int ipos) {
        const auto& ins = insts.at(i);
        auto it = ins.find(ipos);
        return it != ins.end() && it->second.integral;
    };
    auto size{static_cast<int>(confs.size())};
    for (int i = 0; i < size - 1; ++i) {
        for (int ipos = 0; ipos < static_cast<int>(confs[i].segs.size()); ++ipos) {
            if (integral(i, ipos)) continue;
            auto find{false};
            for (int j = i + 1; j < size; ++j) {
                for (int jpos = 0; jpos < static_cast<int>(confs[j].segs.size()); ++jpos) {
                    if (integral(j, jpos)) continue;
                    e_seg seg;
                    auto cl{e_collinear(confs[i].segs[ipos], confs[j].segs[jpos], seg)};
                    if (cl >= 0) {
                        auto minH{std::min(confs[i].h, confs[j].h)};
                        find = cl == 1;
                        auto fun = [&](SegmentInsts& ins) {
                            if (find) {
                                ins.integral = true;
                                ins.minH = minH;
                            } else {
                                InstUnit siu;
                                siu.minH = minH;
                                siu.seg = seg;
                                ins.units.emplace_back(std::move(siu));
                            }
                        };
                        fun(insts[i][ipos]);
                        fun(insts[j][jpos]);
                        if (find) break;
                    }
                }
                if (find) break;
            }
        }
    }

    auto id{0};
    auto phdId{0}; 
    for (const auto& conf : confs) {
        cxd::Vertices vts;
        for (auto& s : conf.segs) {
            vts.emplace_back(s[1]);
        }
        cxd::ConcavePolygon cpolys(vts);
        cpolys.convexDecomp();
        std::vector<cxd::ConcavePolygon> subPolys;
        cpolys.lowestLevelPolys(subPolys);
        for (const auto& sp : subPolys) {
            if (sp.valid()) {
                Polyhedron phd;
                phd.floor = conf.floor;
                if (g_GP.stretch(phd, vts, sp, insts[id], conf.h, conf.te) && phd.init(phdId, id)) {
                    b2pids[id].emplace_back(phdId);
                    phds.emplace_back(std::move(phd));
                    ++phdId;
                } else {
                    std::cout << "init err"  << phdId << " " << id << std::endl;
                }
            } 
            // else {
            //     std::cout << "subPolygons err " << conf.name << " " << phdId << " " << id << std::endl;
            //     sp.show();
            // }
        }
        ++id;
    }
    return true;
}

bool Config::readMaterials(const std::filesystem::path& file, 
    const NameIds& nameIds, BMaterials& bms) const {
    if (!std::filesystem::exists(file)) return false;
    std::ifstream f(file);
    std::string l;
    std::getline(f, l);
    while (std::getline(f, l)) {
        if (sfunc::trimr(l)) continue;
        auto&& us = sfunc::split(l, ",");
        if (us.size() < 6) continue;
        auto it = nameIds.find(std::string(us[0]));
        if (it == nameIds.end()) continue;
        auto id{it->second};
        if (id < 0) {
            id = -1 - id;
        }
        auto& ms = bms[id];
        ms.G_win = 1.0 - stod_s(us[2]);
        ms.r_roof = stod_s(us[4]);
        ms.r_wall = stod_s(us[6]);
    }
    return !bms.empty();
}

int Config::meshSrf(SrfMesh& mesh, const Surface& srf, const ep3& normi, 
    const ep3& center, double pnorm, double pinner) const {
    if (mesh.dir != srf._dir) return -1;
    if (!e_p_equal(normi, srf._normi, 3, pnorm)) return -1;
    auto find{0};
    ep3 p = srf._rmat * center;
    mesh.m.c_rot = bg_point(p(0), p(1));
    if (boost::geometry::within(mesh.m.c_rot, srf._poly)) {
        find = 1;
    } else {
        for (int i = 0; i < 4; ++i) {
            bg_point pi(mesh.m.c_rot.x() + ((i & 1) ? pinner : -pinner), mesh.m.c_rot.y() + ((i & 2) ? pinner : -pinner));
            if (boost::geometry::within(pi, srf._poly)) {
                mesh.m.c_rot = pi;
                find = 1;
                break;
            }
        }
    }
    if (find == 1 && srf.id != c_roof_Id) {
        bg_point pc{center(0), center(1)};
        if (boost::geometry::within(pc, srf.gseg)) return 1;
        auto d{boost::geometry::distance(pc, srf.gseg)};
        if (d < c_5e_2) return 1;
        find = 0;
    }
    return find;
}

void Config::meshType(SrfMesh& smesh, std::string_view typeName) const {
    auto type{0};
    for (int i = 0; i < static_cast<int>(c_srfType.size()); ++i) {
        if (c_srfType[i] == typeName) {
            type = i;
            break;
        }
    }
    smesh.m.wall = type != 1;
    smesh.setDir();
}

bool Config::readMeshs_(std::vector<SrfMesh>& bms, const std::vector<std::string_view>& us, 
    const GPhds& phds, const std::vector<int>& pids, double pnorm, double pinner) const {
    SrfMesh smesh;
    for (int i = 0; i < static_cast<int>(c_orientations.size()); ++i) {
        if (us[2] == c_orientations[i]) {
            smesh.dir = i;
            break;
        }
    }
    ep3 center(stod_s(us[5]), stod_s(us[6]), stod_s(us[7]));
    ep3 normi(stod_s(us[8]), stod_s(us[8]), stod_s(us[10]));
    smesh.m.area = stod_s(us[11]) * 1e-3;
    auto find{false};
    auto minDis{-1.0};
    auto sid{0};
    for (auto pid : pids) {
        auto& pit = phds.at(pid);
        for (auto& srf : pit.srfs) {
            auto in{meshSrf(smesh, srf, normi, center, pnorm, pinner)};
            if (in < 0) continue;
            if (in == 1) {
                find = true;
                smesh.phd = pid;
                smesh.srf = &srf;
                break;
            }
            auto dis = boost::geometry::distance(smesh.m.c_rot, srf._poly);
            if (minDis < 0.0 || minDis > dis) {
                minDis = dis;
                smesh.phd = pid;
                sid = srf.id;
            }
        }
        if (find) break;
    }
    if (!find) {
        auto& srf = phds.at(smesh.phd).srfs.at(sid);
        smesh.srf = &srf;
        boost::geometry::centroid(srf._poly, smesh.m.c_rot);
    }
    meshType(smesh, us[12]);
    bms.emplace_back(std::move(smesh));
    return true;
}

bool Config::_readMeshs(std::vector<SrfMesh>& bms, const std::vector<std::string_view>& us, 
    const GPhds& phds) const {
    auto phdId{stoi_s(us[13])};
    auto sid{stoi_s(us[14])};
    if (phdId < 0 || phdId >= static_cast<int>(phds.size())) return false;
    auto& pit = phds[phdId];
    if (sid < 0 || sid >= static_cast<int>(pit.srfs.size())) return false;
    auto& srf = pit.srfs[sid];
    SrfMesh smesh;
    smesh.phd = phdId;
    smesh.srf = &srf;
    ep3 center(stod_s(us[5]), stod_s(us[6]), stod_s(us[7]));
    ep3 p = srf._rmat * center;
    smesh.m.c_rot = bg_point(p(0), p(1));
    smesh.m.area = stod_s(us[11]) * 1e-3;
    meshType(smesh, us[12]);
    auto&& pus = sfunc::split(us[15], " ");
    if (pus.size() % 2 == 0) {
        for (int i = 0; i < static_cast<int>(pus.size()); i += 2) {
            smesh.m.poly.outer().emplace_back(bg_point(stod_s(pus[i]), stod_s(pus[i + 1])));
        }
    }
    bms.emplace_back(std::move(smesh));
    return true;
}

bool Config::readMeshs(SrfMeshs& meshs, const std::filesystem::path& dir, const GPhds& phds, 
    const OneMany& b2pids, const NameIds& nameIds, double pnorm, double pinner) const {
    std::ranges::for_each(std::filesystem::directory_iterator{dir}, 
        [&](const auto& dir_entry) {
            if (!dir_entry.is_regular_file()) return;
            auto file{dir_entry.path()};
            auto name{file.filename().string()};
            auto pos{name.find('_')};
            if (pos == std::string::npos) return;
            auto nstr{name.substr(pos + 1)};
            if (nstr != "geometry.csv") return;
            auto bstr{name.substr(0, pos)};
            auto it = nameIds.find(bstr);
            if (it == nameIds.end() || it->second < 0) return;
            auto buildId{it->second};
            auto bpit = b2pids.find(buildId);
            if (bpit == b2pids.end()) return;
            auto& bms = meshs[buildId];
            std::ifstream f(file);
            std::string l{};
            std::getline(f, l);
            while (std::getline(f, l)) {
                if (sfunc::trimr(l)) continue;
                auto&& us = sfunc::split(l, ",");
                auto size{static_cast<int>(us.size())};
                if (size == 13) {
                    readMeshs_(bms, us, phds, bpit->second, pnorm, pinner);
                } else if (size == 16) {
                    _readMeshs(bms, us, phds);
                }
            }
            if (bms.empty()) {
                std::cout << "readCuts err cuts empty " << bstr << " " << pinner << std::endl;
            }
        }
    );
    return true;   
}

void Config::writeMeshs(const SrfMeshs& smeshs, const std::filesystem::path& fdir, const BFlags& flags) const {
    static constexpr std::string_view title{"BUILDING,SURFACE,orientation,intersection,terrain_elevation,Xcoor,Ycoor,Zcoor,Xdir,Ydir,Zdir,AREA_m2,TYPE,phd,srf,poly\n"};
    static constexpr std::string_view row{"{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"};
    static constexpr std::string_view meshName{"{}_geometry.csv"}; 
    for (auto&[id, meshs] : smeshs) {
        auto it = flags.find(id);
        if (it == flags.end()) continue;
        auto file{fdir / std::format(meshName, it->second.name)};
        std::ofstream ofs(file);
        ofs << title;
        auto te{it->second.te};
        auto sid{0};
        const Surface* srf{};
        e_rmat rmat_inv{Eigen::Matrix3d::Identity()};
        for (auto& m : meshs) {
            if (srf != m.srf) {
                srf = m.srf;
                rmat_inv = srf->rmat_inv();
            }
            auto ori = c_orientations[m.srf->_dir];
            int tid{!m.m.wall};
            if (m.srf->id == c_roof_Id) {
                tid = 2;
            }
            auto type = c_srfType[tid];
            ep3 c = rmat_inv * ep3(m.m.c_rot.x(), m.m.c_rot.y(), m.srf->_rotz); //m.m.c
            auto& norm = m.srf->_normi;
            auto sname{std::format(c_srf_name, sid)};
            ofs << std::format(row, it->second.name, sname, ori, 0, te, c(0), c(1), c(2), 
                norm(0), norm(1), norm(2), m.m.area * 1e3, type, m.phd, m.srf->id, m.m.poly2str());;
            ++sid;
        }
        ofs << std::endl;
    }
}
