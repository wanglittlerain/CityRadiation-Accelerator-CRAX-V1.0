#pragma once
#ifndef _SPLIT_STRING_
#define _SPLIT_STRING_
#include <vector>
#include <string>
#include <string_view>
#include <ranges>
#include <fstream>
namespace sfunc {
inline auto split(std::string_view str, std::string_view delim) {
    auto v = str | std::views::split(delim) | std::views::transform([](auto&& unit) {
        return std::string_view(unit.begin(), unit.end());
    });
    return std::vector<std::string_view>{v.begin(), v.end()};
}

inline bool trimr(std::string& str) {
    if (!str.empty() && str.back() == '\r') {
        str.pop_back();
    }
    return str.empty();
}

inline void skip(std::ifstream& file, int n) {
    for (int i = 0; i < n; ++i) {
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}
}
#endif