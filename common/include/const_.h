#pragma once
#ifndef _CONSTANTS_
#define _CONSTANTS_
static constexpr auto c_1e_8{1e-8};
static constexpr auto c_1e_2{1e-2};
static constexpr auto c_5e_2{5e-2};
static constexpr auto c_pi{3.14159265359};
static constexpr auto c_roof_Id{0};
static constexpr auto c_virtual_srf{-1};
static constexpr auto c_declination{2.0 * c_pi / 365.0};
static constexpr auto c_phi{6};    // Number of altitude angle steps for sky integration
static constexpr auto c_theta{24}; // Number of azimuth angle steps for sky integration
static constexpr size_t c_diffuse_step{c_phi * c_theta};

enum enShadowType : int {
    c_shadow_init,
    c_shadow_none,
    c_shadow_full,
    c_shadow_other_set
};
#endif