#pragma once
#ifndef _ENV_H_
#define _ENV_H_
#include <format>
#include <fstream>
#include "nlohmann/json.hpp"
#include "const_.h"
#include "bg_.h"
#include "eigen_.h"
#include <chrono>
#include <array>
#include <iostream>
using json = nlohmann::json;

class Elapsed {
public:
    Elapsed(const std::string_view msg = "") {
        _all = _run = std::chrono::system_clock::now();
        _msg = msg;
    }

    ~Elapsed() {
        if (!_msg.empty()) {
            view(_msg, false);
        }
    }

    inline void view(const std::string_view str, bool run = true) {
        auto now{std::chrono::system_clock::now()};
        std::chrono::microseconds d{};
        if (run) {
            d = std::chrono::duration_cast<std::chrono::microseconds>(now - _run);
            _run = now;
        } else {
            d = std::chrono::duration_cast<std::chrono::microseconds>(now - _all);
        }
        std::cout << std::format("{}{}s\n", str, d.count() / 1000000);
    }
private:
    std::chrono::system_clock::time_point _all;
    std::chrono::system_clock::time_point _run;
    std::string _msg{};
};
#endif