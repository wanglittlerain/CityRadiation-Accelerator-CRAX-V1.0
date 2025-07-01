#pragma once
#ifndef _CONCAVE_H_
#define _CONCAVE_H_
//https://github.com/mjjq/ConvexDecomposition
//Optimization Schemes https://doc.cgal.org/latest/Partition_2/index.html
#include "bg_.h"
#include "eigen_.h"
namespace cxd {
inline double signedArea(const ep2& v1, const ep2& v2) {
    return (v2.x() - v1.x()) * (v2.y() + v1.y());
}

inline double handedness(const ep2& v1, const ep2& v2, const ep2& v3) {
    return e_cross(ep2{v2 - v1}, ep2{v3 - v2});
}

inline bool intersects(ep2& interPos, const e_seg& s1, const e_seg& s2) {
    static constexpr auto TOLERANCE{1e-2};
    auto cd1d2{e_cross(s1[2], s2[2])};
    if (std::abs(cd1d2) < 1e-6)
       return false;

    auto t1{e_cross(s2[0] - s1[0], s2[2]) / cd1d2};
    if (t1 < -TOLERANCE || t1 > (1.0 + TOLERANCE))
        return false;

    interPos = s1[0] + s1[2] * t1;
    ep2 p1{interPos - s2[0]};
    auto t2{p1.dot(s2[2])};
    return !(t2 < -TOLERANCE || t2 / s2[2].dot(s2[2]) >= 1.0 - TOLERANCE);
}

struct SliceVertex {
    ep2 pos{};
    int index{};
    double dis{};
};

using Vertices = std::vector<ep2>;
class ConcavePolygon {
    static constexpr int maxTime{16};
    Vertices _verts;
    int _size{};
    std::vector<ConcavePolygon> _sons;

    size_t mod(int x, int m) const {
        return (x + m) % m;
    }

    void flip() {
        for (int i = 1; i < (_size + 1) / 2; ++i) {
            std::swap(_verts[i], _verts[_size - i]);
        }
    }

    bool rightHanded() const {
        if (_size < 3) return false;
        double sarea{};
        for (int i = 0; i < _size; ++i) {
            sarea += signedArea(_verts[i], _verts[mod(i + 1, _size)]);
        }
        return sarea < 0.0;
    }

    bool vertexInCone(e_seg const& ls1, e_seg const& ls2, const ep2& origin, const ep2& vert) const {
        ep2 relative{vert - origin};
        return e_cross(relative, ls1[2]) < 0.0 && e_cross(relative, ls2[2]) > 0.0;
    }

    std::vector<int> findVerticesInCone(e_seg const& ls1, e_seg const& ls2, const ep2& origin) const {
        std::vector<int> result;
        for (int i = 0; i < _size; ++i) {
            if (vertexInCone(ls1, ls2, origin, _verts[i])) {
                result.emplace_back(i);
            }
        }
        return result;
    }

    bool checkVisibility(const ep2& original, const ep2& vert) const {
        auto vts = verticesAlongLineSegment(e_seg(original, vert));
        return vts.size() <= 3;
    }

    int bestVertexToConnect(std::vector<int> const& indices, const ep2& origin) const {
        if (indices.size() == 1) {
            if (checkVisibility(origin, _verts[indices[0]]))
                return indices[0];
        } else if (indices.size() > 1) {
            auto firstIdx{-1};
            for (auto& index : indices) {
                auto& prev{_verts[mod(index - 1, _size)]};
                auto& curr{_verts[index]};
                auto& next{_verts[mod(index + 1, _size)]};
                if (handedness(prev, curr, next) < 0.0 && 
                    checkVisibility(origin, curr)) {
                    if (firstIdx < 0) {
                        firstIdx = index;
                    }
                    if (vertexInCone(e_seg(prev, curr), e_seg(next, curr), curr, origin)) {
                        return index;
                    }
                }            
            }
            if (firstIdx >= 0) {
                return firstIdx;
            }
            auto dis{1e+8};
            auto closest{indices[0]};
            for (auto& index : indices) {
                ep2 pos{_verts[index] - origin};
                auto currDistance{pos.dot(pos)};
                if (currDistance < dis) {
                    dis = currDistance;
                    closest = index;
                }
            }
            return closest;
        }
        return -1;
    }

    int firstReflexVertex() const {
        for (int i = 0; i < _size; ++i) {
            if (handedness(_verts[mod(i - 1, _size)], _verts[i], _verts[mod(i + 1, _size)]) < 0.0)
                return i;
        }
        return -1;
    }

    using VertexMap = std::map<int, ep2>;
    VertexMap cullByDistance(VertexMap const& input, const ep2& origin, int keep) const {
        if (keep >= static_cast<int>(input.size())) return input;
        std::vector<SliceVertex> sliceVertices;
        for (auto& [id, invert] : input) {
            SliceVertex vert;
            vert.pos = invert;
            vert.index = id;
            ep2 p{invert - origin};
            vert.dis = p.dot(p);
            sliceVertices.emplace_back(std::move(vert));
        }
        for (size_t i = 1; i < sliceVertices.size(); ++i)
            for (size_t j = i; j > 0 && sliceVertices[j].dis < sliceVertices[j - 1].dis; --j)
                std::swap(sliceVertices[j], sliceVertices[j - 1]);

        VertexMap result;
        auto idx{0};
        for (const auto& svert : sliceVertices) {
            if (idx >= keep) break;
            result.insert({svert.index, svert.pos});
            ++idx;
        }
        return result;
    }

    VertexMap verticesAlongLineSegment(e_seg const& segment) const {
        VertexMap result;
        for (int i = 0; i < _size; ++i) {
            e_seg temp(_verts[i], _verts[mod(i + 1, _size)]);
            ep2 interPos;
            if (intersects(interPos, segment, temp)) {
                result.insert({i, interPos});
            }
        }
        return result;
    }

public:
    ConcavePolygon(Vertices const& verts) : _verts{verts} {
        _size = static_cast<int>(_verts.size());
        if (_size > 2 && !rightHanded()) {
            flip();   
        }
    }
    ConcavePolygon() = default;

    void slice(int vertex1, int vertex2) {
        if (vertex1 > vertex2)
            std::swap(vertex1, vertex2);
        
        if (vertex2 - vertex1 <= 1) return;
        Vertices returnVerts;
        Vertices newVerts;
        for (int i = 0; i < _size; ++i) {
            if (i == vertex1 || i == vertex2) {
                returnVerts.emplace_back(_verts[i]);
                newVerts.emplace_back(_verts[i]);
            } else if(i > vertex1 && i < vertex2) {
                returnVerts.emplace_back(_verts[i]);
            } else
                newVerts.emplace_back(_verts[i]);
        }
        _sons.emplace_back(std::move(returnVerts));
        _sons.emplace_back(std::move(newVerts));
    }

    void slice(const e_seg& segment) {
        if (_sons.size() > 0) {
            _sons[0].slice(segment);
            _sons[1].slice(segment);
            return;
        }
        static constexpr auto TOLERANCE{1e-3};
        auto slicedVertices{verticesAlongLineSegment(segment)};
        slicedVertices = cullByDistance(slicedVertices, segment[0], 2);
        if (slicedVertices.size() < 2) return;
        Vertices leftVerts;
        Vertices rightVerts;
        for (int i = 0; i < _size; ++i) {
            ep2 relative = _verts[i] - segment[0];
            auto begin{slicedVertices.begin()};
            auto it{slicedVertices.find(i)};
            double perpDistance{std::abs(e_cross(relative, segment[2]))};
            if (perpDistance > TOLERANCE || (perpDistance <= TOLERANCE && (it == slicedVertices.end()))) {
                if ((i > begin->first) && (i <= (++begin)->first)) {
                    leftVerts.emplace_back(_verts[i]);
                } else {
                    rightVerts.emplace_back(_verts[i]);
                }
            }
            if (it != slicedVertices.end()) {
                rightVerts.emplace_back(it->second);
                leftVerts.emplace_back(it->second);
            }
        }
        _sons.emplace_back(std::move(leftVerts));
        _sons.emplace_back(std::move(rightVerts));
    }

    void convexDecomp(int time = 0) {
        if (_size > 3 && _sons.empty() && time < maxTime) {
            auto reflexIndex{firstReflexVertex()};
            if (reflexIndex == -1) return;
            auto& prev{_verts[mod(reflexIndex - 1, _size)]};
            auto& curr{_verts[reflexIndex]};
            auto& next{_verts[mod(reflexIndex + 1, _size)]};
            e_seg ls1(prev, curr);
            e_seg ls2(next, curr);
            auto&& vertsInCone{findVerticesInCone(ls1, ls2, curr)};
            auto bestVert{-1};
            if (!vertsInCone.empty()) {
                bestVert = bestVertexToConnect(vertsInCone, curr);
                if (bestVert != -1) {
                    slice(e_seg(curr, _verts[bestVert]));
                }
            }
            if (vertsInCone.empty() || bestVert == -1) {
                slice(e_seg(curr, (ls1[2] + ls2[2]) * 1e+10));
            }
            ++time;
            for (auto& son : _sons) {
                son.convexDecomp(time);
            }
        } 
    }

    const Vertices& outer() const {
        return _verts;
    }

    void lowestLevelPolys(std::vector<ConcavePolygon>& returnArr) {
        if (_sons.size() > 0) {
            _sons[0].lowestLevelPolys(returnArr);
            _sons[1].lowestLevelPolys(returnArr);
        } else
            returnArr.push_back(*this);
    }

    void show() const {
        for (const auto& it : _verts) {
            std::cout << it.x() << " " << it.y() << "\n";
        }
        std::cout << "------\n";
    }

    bool valid() const {
        bg_polygon poly;
        for (const auto& it : _verts) {
            poly.outer().emplace_back(bg_point(it.x(), it.y()));
        }
        if (!boost::geometry::is_valid(poly)) {
            boost::geometry::correct(poly);
        }
        if (boost::geometry::is_valid(poly)) {
            auto area{boost::geometry::area(poly)};
            if (area > c_1e_2) return true;
        }
        return false;
    }
};
}
#endif // CONCAVE_POLY_H
