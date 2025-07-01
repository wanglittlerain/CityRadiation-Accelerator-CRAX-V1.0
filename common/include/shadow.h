#pragma once
#ifndef _SHADOW_
#define _SHADOW_
#include "base.h"
class  Shadow {
public:
    Shadow() = default;
    void calculate(GPhds& phds, Solar& s, const ShadowSet& cset) const;
    void diffuse(GPhds& phds, SrfMeshs& smeshs, const ShadowSet& cset, int b = 0, int e = 0) const;
private:
    void cfun(GPhds& phds, double sh, const ep3& sv, const ShadowSet& cset, bool save, bool clear = false) const;
    void init(GPhds& phds, const double sh, const ep3& sv, bool save, bool clear) const;
    double ndis(const Polyhedron& phd1, const Polyhedron& phd2, double tansh) const;
    bool check(const ebbox& box1, const ebbox& box2, const ep3& sv) const;
    void add(NBHDs& res, const Surfaces& nsrfs) const;
    void checkAdd(Polyhedron& phd1, const Polyhedron& phd2, const ep3& sv, double dis, bool dfuse) const;
    void neighborhoods(double tansh, const ep3& sv, GPhds& phds, const ShadowSet& cset) const;

    bool _projection(bg_polygon& poly, const Surface& srf, const ep3& sv, 
        const epts& points, double fh, double te) const;
    void forward(Surface& srf, double te, const NBHDs& nbhds, const ep3& sv, int type, bool save) const;
    void srfcalc(Surface& srf, int type, bool save) const;

    bool difference(Polys& res, const bg_polygon& totalPoly, const bg_polygon& poly) const;
    double substract(const bg_polygon& totalPoly, const Polys& intersectPolys) const;
};
static const Shadow g_shadow;
#endif