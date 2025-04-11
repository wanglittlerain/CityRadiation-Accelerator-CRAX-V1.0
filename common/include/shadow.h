#pragma once
#ifndef _SHADOW_
#define _SHADOW_
#include "base.h"
class  Shadow {
public:
    Shadow() = default;
    void calc_(Solar& s, GPhds& phds, const ShadowSet& cset) const;

    void diffuse(GPhds& phds, SrfMeshs& smeshs, const ShadowSet& cset, bool alone) const;
private:
    void cfun(double sh, const ep3& sv, GPhds& phds, const ShadowSet& cset, bool save = true) const;

    void clear(GPhds& phds) const;
    void init(const ep3& sv, GPhds& phds, bool save) const;

    void neighbours(double sh, const ep3& sv, GPhds& phds, const ShadowSet& cset) const;

    bool ndis(const Polyhedron& phd1, const Polyhedron& phd2, double sh, double& dis) const;

    bool check(const ebbox& box1, const ebbox& box2, const double tanSh, 
        const ep3& sv, const double dis, bool noSolar) const;
    
    void add(std::vector<const Surface*>& res, const Polyhedron& other) const;

    void calculate(const ep3& sv, GPhds& phds, int type, bool mode, bool save) const;

    bool _projection(bg_polygon& poly, const Surface& srf, const ep3& sv, 
        double ndot, const epts& points, double fh) const;

    void forward(Surface& srf, Polyhedron& phd, const ep3& sv, int type, bool save) const;

    void srfcalc(Surface& srf, int type, bool save) const;

    bool difference(Polys& res, const bg_polygon& totalPoly, const bg_polygon& poly) const;
    double substract(const bg_polygon& totalPoly, const Polys& intersectPolys) const;

    bool projection_(const Surface* nsrf, const ep3& sv, const ep3& p) const;
    void backward(Surface& srf, Polyhedron& phd, const ep3& sv, int type, bool save) const;
};
static const Shadow sShdg;
#endif