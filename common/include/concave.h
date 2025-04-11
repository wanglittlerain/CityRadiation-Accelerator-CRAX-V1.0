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
    ep2 edge1{v2 - v1};
    ep2 edge2{v3 - v2};
    return e_cross(edge1, edge2);
}

struct LineSegment {
    e_seg seg;

    LineSegment() = default;
    LineSegment(const ep2& start, const ep2& final) {
        seg[0] = start;
        seg[1] = final;
    }

    ep2 direction() const {
        return seg[1] - seg[0];
    }
};
inline bool intersects(ep2& interPos, const LineSegment& s1, const LineSegment& s2) {
    static constexpr auto TOLERANCE{1e-2};
    ep2 d1{s1.direction()};
    ep2 d2{s2.direction()};
    auto cd1d2{e_cross(d1, d2)};
    if (std::abs(cd1d2) < 1e-6)
       return false;

    auto t1{e_cross(s2.seg[0] - s1.seg[0], d2) / cd1d2};
    if (t1 < -TOLERANCE || t1 > (1.0 + TOLERANCE))
        return false;

    interPos = s1.seg[0] + d1 * t1;
    ep2 p1{interPos - s2.seg[0]};
    ep2 p2{s2.seg[1] - s2.seg[0]};
    auto t2{p1.dot(p2)};
    if (t2 < -TOLERANCE || t2 / p2.dot(p2) >= 1.0 - TOLERANCE)
        return false;

    return true;
}

struct SliceVertex {
    ep2 pos{};
    int index{};
    double dis{};
};

using Vertices = std::vector<ep2>;
class ConcavePolygon {
    using PolygonArray = std::vector<ConcavePolygon>;
    static constexpr int maxTime{16};
    Vertices vertices;
    PolygonArray subPolygons;

    size_t mod(int x, int m) {
        int r = x % m;
        return static_cast<size_t>(r < 0 ? r + m : r);
    }

    void flip(Vertices& _verts) {
        auto iMax{static_cast<int>(_verts.size() / 2 + _verts.size() % 2)};
        for (int i = 1; i < iMax; ++i) {
            std::swap(_verts[i], _verts[_verts.size() - i]);
        }
    }

    bool rightHanded(Vertices& _verts) {
        if (_verts.size() < 3ull)
            return false;

        double sarea{};
        auto vsize{static_cast<int>(_verts.size())};
        for (int i = 0; i < vsize; ++i) {
            sarea += signedArea(_verts[i], _verts[mod(i + 1, vsize)]);
        }
        if (sarea < 0.0)
            return true;

        return false;
    }

    bool vertexInCone(LineSegment const& ls1, LineSegment const& ls2, 
        const ep2& origin, const ep2& vert) {
        ep2 relative{vert - origin};
        auto ls1Product{e_cross(relative, ls1.direction())};
        auto ls2Product{e_cross(relative, ls2.direction())};
        if (ls1Product < 0.0 && ls2Product > 0.0)
            return true;
        return false;
    }

    std::vector<int> findVerticesInCone(LineSegment const& ls1, LineSegment const& ls2, 
        const ep2& origin, const Vertices& inputVerts) {
        std::vector<int> result;
        for (int i = 0; i < static_cast<int>(inputVerts.size()); ++i) {
            if (vertexInCone(ls1, ls2, origin, inputVerts[i])) {
                result.emplace_back(i);
            }
        }
        return result;
    }

    bool checkVisibility(const ep2& original, const ep2& vert, const Vertices& polygonVertices) {
        LineSegment ls(original, vert);
        auto intersectingVerts = verticesAlongLineSegment(ls, polygonVertices);
        if (intersectingVerts.size() > 3)
            return false;
        return true;
    }

    int bestVertexToConnect(std::vector<int> const& indices, 
        const Vertices& polygonVertices, const ep2& origin) {
        if (indices.size() == 1) {
            if (checkVisibility(origin, polygonVertices[indices[0]], polygonVertices))
                return indices[0];
        } else if (indices.size() > 1) {
            auto size{static_cast<int>(polygonVertices.size())};
            auto firstIdx{-1};
            for (const auto& index : indices) {
                const auto& prev{polygonVertices[mod(index - 1, size)]};
                const auto& curr{polygonVertices[index]};
                const auto& next{polygonVertices[mod(index + 1, size)]};
                if (handedness(prev, curr, next) < 0.0 && 
                    checkVisibility(origin, curr, polygonVertices)) {
                    if (firstIdx < 0) {
                        firstIdx = index;
                    }
                    LineSegment ls1(prev, curr);
                    LineSegment ls2(next, curr);
                    if (vertexInCone(ls1, ls2, curr, origin)) {
                        return index;
                    }
                }            
            }
            if (firstIdx >= 0) {
                return firstIdx;
            }
            auto minDistance{1e+8};
            auto closest{indices[0]};
            for (const auto& index : indices) {
                ep2 pos{polygonVertices[index] - origin};
                auto currDistance{pos.dot(pos)};
                if (currDistance < minDistance) {
                    minDistance = currDistance;
                    closest = index;
                }
            }
            return closest;
        }
        return -1;
    }

    void convexDecomp(Vertices const& _vertices, int time) {
        auto reflexIndex{firstReflexVertex(_vertices)};
        if (reflexIndex == -1)
            return;
        auto size{static_cast<int>(_vertices.size())};
        const auto& prev{_vertices[mod(reflexIndex - 1, size)]};
        const auto& curr{_vertices[reflexIndex]};
        const auto& next{_vertices[mod(reflexIndex + 1, size)]};
        LineSegment ls1(prev, curr);
        LineSegment ls2(next, curr);
        auto&& vertsInCone{findVerticesInCone(ls1, ls2, curr, _vertices)};
        auto bestVert{-1};
        if (!vertsInCone.empty()) {
            bestVert = bestVertexToConnect(vertsInCone, _vertices, curr);
            if (bestVert != -1) {
                LineSegment newLine(curr, _vertices[bestVert]);
                slice(newLine);
            }
        }
        if (vertsInCone.empty() || bestVert == -1) {
            LineSegment newLine(curr, (ls1.direction() + ls2.direction()) * 1e+10);
            slice(newLine);
        }
        ++time;
        for (auto& subPoly : subPolygons) {
            subPoly.convexDecomp(time);
        }
    }

    int firstReflexVertex(Vertices const& _vertices) {
        auto size{static_cast<int>(_vertices.size())};
        for (int i = 0; i < size; ++i) {
            auto hd{handedness(_vertices[mod(i - 1, size)], 
                _vertices[i], _vertices[mod(i + 1, size)])};
            if (hd < 0.0)
                return i;
        }
        return -1;
    }

    void flip() {
        flip(vertices);
    }

    using VertexMap = std::map<int, ep2>;
    VertexMap cullByDistance(VertexMap const& input,
        const ep2& origin, int const& maxVertsToKeep) {
        if (maxVertsToKeep >= static_cast<int>(input.size()))
            return input;

        std::vector<SliceVertex> sliceVertices;
        for (const auto& [id, invert] : input) {
            SliceVertex vert;
            vert.pos = invert;
            vert.index = id;
            ep2 p{invert - origin};
            vert.dis = p.dot(p);
            sliceVertices.emplace_back(vert);
        }
        for (size_t i = 1; i < sliceVertices.size(); ++i)
            for (size_t j = i; j > 0 && sliceVertices[j].dis < sliceVertices[j - 1].dis; --j)
                std::swap(sliceVertices[j], sliceVertices[j - 1]);

        VertexMap result;
        auto idx{0};
        for (const auto& svert : sliceVertices) {
            if (idx >= maxVertsToKeep) {
                break;
            }
            result.insert({svert.index, svert.pos});
            ++idx;
        }
        return result;
    }

    VertexMap verticesAlongLineSegment(LineSegment const& segment,const Vertices& _vertices) {
        VertexMap result;
        LineSegment temp;
        auto size{static_cast<int>(_vertices.size())};
        for (int i = 0; i < size; ++i) {
            temp.seg[0] = _vertices[i];
            temp.seg[1] = _vertices[mod(i + 1, size)];
            ep2 interPos;
            if (intersects(interPos, segment, temp)) {
                result.insert({i, interPos});
            }
        }
        return result;
    }

public:
    ConcavePolygon(Vertices const& _vertices) : vertices{_vertices} {
        if (vertices.size() > 2) {
            if (!rightHanded())
                flip();
        }
    }
    ConcavePolygon() = default;

    bool rightHanded() {
        return rightHanded(vertices);
    }

    void slice(int vertex1, int vertex2) {
        if (vertex1 > vertex2)
            std::swap(vertex1, vertex2);
        
        if (vertex2 - vertex1 <= 1) return;

        Vertices returnVerts;
        Vertices newVerts;
        for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
            if (i == vertex1 || i == vertex2) {
                returnVerts.emplace_back(vertices[i]);
                newVerts.emplace_back(vertices[i]);
            } else if(i > vertex1 && i < vertex2) {
                returnVerts.emplace_back(vertices[i]);
            } else
                newVerts.emplace_back(vertices[i]);
        }
        subPolygons.emplace_back(std::move(returnVerts));
        subPolygons.emplace_back(std::move(newVerts));
    }

    void slice(const LineSegment& segment) {
        if (subPolygons.size() > 0) {
            subPolygons[0].slice(segment);
            subPolygons[1].slice(segment);
            return;
        }
        static constexpr auto TOLERANCE{1e-3};
        auto slicedVertices{verticesAlongLineSegment(segment, vertices)};
        slicedVertices = cullByDistance(slicedVertices, segment.seg[0], 2);
        if (slicedVertices.size() < 2)
            return;

        Vertices leftVerts;
        Vertices rightVerts;
        for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
            ep2 relative = vertices[i] - segment.seg[0];
            auto begin{slicedVertices.begin()};
            auto it{slicedVertices.find(i)};
            double perpDistance{std::abs(e_cross(relative, segment.direction()))};
            if (perpDistance > TOLERANCE || (perpDistance <= TOLERANCE && (it == slicedVertices.end()))) {
                if ((i > begin->first) && (i <= (++begin)->first)) {
                    leftVerts.emplace_back(vertices[i]);
                } else {
                    rightVerts.emplace_back(vertices[i]);
                }
            }
            if (it != slicedVertices.end()) {
                rightVerts.emplace_back(it->second);
                leftVerts.emplace_back(it->second);
            }
        }
        subPolygons.emplace_back(std::move(leftVerts));
        subPolygons.emplace_back(std::move(rightVerts));
    }

    void convexDecomp(int time = 0) {
        if (vertices.size() > 3 && subPolygons.empty() && time < maxTime) {
            convexDecomp(vertices, time);
        } 
    }

    const Vertices & outer() const {
        return vertices;
    }

    void lowestLevelPolys(std::vector<ConcavePolygon>& returnArr) {
        if (subPolygons.size() > 0) {
            subPolygons[0].lowestLevelPolys(returnArr);
            subPolygons[1].lowestLevelPolys(returnArr);
        } else
            returnArr.push_back(*this);
    }

    void reset() {
        if (subPolygons.size() > 0) {
            subPolygons[0].reset();
            subPolygons[1].reset();
            subPolygons.clear();
        }
    }

    void show() const {
        for (const auto& it : vertices) {
            std::cout << it.x() << " " << it.y() << "\n";
        }
        std::cout << "------" << std::endl;
    }

    bool valid() const {
        bg_polygon poly;
        for (const auto& it : vertices) {
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