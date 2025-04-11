#pragma once
#ifndef _RAD_FUNC_
#define _RAD_FUNC_
#include <array>
#include <cmath>
#include "const_.h"
static constexpr std::array EpsilonLimit{1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2};
static constexpr std::array<std::array<double, 6>, 8> FR = {{
    {-0.0083117, 0.5877285, -0.0620636, -0.0596012, 0.0721249, -0.0220216},
    {0.1299457, 0.6825954, -0.1513752, -0.0189325, 0.0659650, -0.0288748},
    {0.3296958, 0.4868735, -0.2210958, 0.0554140, -0.0639588, -0.0260542},
    {0.5682053, 0.1874525, -0.295129, 0.1088631, -0.1519229, -0.0139754},
    {0.873028, -0.3920403, -0.3616149, 0.2255647, -0.4620442, 0.0012448},
    {1.1326077, -1.2367284, -0.4118494, 0.2877813, -0.8230357, 0.0558651},
    {1.0601591, -1.5999137, -0.3589221, 0.2642124, -1.127234, 0.1310694},
    {0.677747, -0.3272588, -0.2504286, 0.1561313, -1.3765031, 0.2506212}
}};

struct StepCoef {
    double dhr{}; 
    double dnr{};
    double c1{};
    double c2{};
    double c3{};
    double c4{};
};

inline void dr_coef(double h, double albedo, StepCoef& c) {
    auto z{0.5 * c_pi - h};
    auto cosz{std::cos(z)};
    c.c1 = 0.0;
    c.c2 = 0.0;
    c.c3 = 0.0;
    if (c.dhr > c_1e_8) {
        auto epsilon{1.0 + c.dnr / ((1.0 + 1.041 * std::pow(z, 3)) * c.dhr)};
        auto idx{0};
        for (; idx < 7; ++idx) {
            if (epsilon < EpsilonLimit[idx]) break;
        }
        auto& arr = FR[idx];
        static constexpr auto m{1 / 1353.0};
        auto derta{c.dhr * m / (std::sin(h) + 0.15 * std::pow(3.885 + h, -1.253))};
        auto f1{std::max(0.0, std::min(1.0, arr[0] + arr[1] * derta + arr[2] * z))};
        auto f2{std::max(0.0, std::min(1.0, arr[3] + arr[4] * derta + arr[5] * z))};
        c.c1 = c.dhr * f2;
        c.c2 = c.dhr * (1 - f1) * 0.5;
        auto coef_b{(std::max)(0.0871557, cosz)};
        c.c3 = c.dhr * f1 / coef_b;
    }
    c.c4 = (c.dnr * cosz + c.dhr) * albedo;
}

inline double diffuse_(double h, const StepCoef& c, double cos_i, 
    double gdot, double sinsita, double rdome, double rhorizon, double sr) {
    if (h <= 0.0 || c.dhr < c_1e_8) return 0.0;
    auto dr{c.c1 * sinsita * rhorizon + c.c2 * (1 + gdot) * rdome + c.c3 * cos_i * (1 - sr)};
    return dr;
}

inline double ground_reflect_(double h, double c4, double gdot, double rdomeg) {
    if (h <= 0.0) return 0.0;
    auto grr{(0.5 * (1 - gdot)) * c4 * rdomeg};
    return grr;
}

inline double reflection_factor(double rdome, double rdomeg, double gdot) {
    return 1 - rdome * (1 + gdot) / 2.0 - rdomeg * (1 - gdot) / 2.0;
}
#endif