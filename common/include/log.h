#pragma once
#ifndef _LOG_H_
#define _LOG_H_
#include <array>
#include <string>
#include <string_view>
#include <list>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <format>
#include <thread>
#include <atomic>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>
#include <chrono>
#ifdef _MSC_VER
	#ifndef NOMINMAX
		#define NOMINMAX
	#endif
	#include <windows.h>
#endif
using uint = unsigned int;
struct system_time {
	uint sec{};
	uint min{};
	uint hour{};
	uint mday{};
	uint mon{};
	uint year{};
	uint wday{};
};

static auto sys_time() {
	system_time ret{};
#ifdef _MSC_VER
	SYSTEMTIME t{};
	GetLocalTime(&t);
	ret.year = t.wYear;
	ret.mon	= t.wMonth;
	ret.mday = t.wDay;
	ret.hour = t.wHour;
	ret.min	= t.wMinute;
	ret.sec	= t.wSecond;
	ret.wday = t.wDayOfWeek;
#else
	auto now = std::chrono::system_clock::now();
	time_t tt = std::chrono::system_clock::to_time_t(now);
	auto ct = localtime(&tt);
	ret.year = ct->tm_year + 1900;
	ret.mon	= ct->tm_mon + 1;
	ret.mday	= ct->tm_mday;
	ret.hour = ct->tm_hour;
	ret.min	= ct->tm_min;
	ret.sec	= ct->tm_sec;
	ret.wday = ct->tm_wday;
#endif
	return ret;
}

template<class T>
struct mt_cv {
	T data;
	std::mutex mtx;
	std::condition_variable cv;
	bool ready{};

	mt_cv() = default;
	mt_cv(mt_cv const& ref) : data(ref.data) {}
	mt_cv(T const& ref) : data(ref) {}
	mt_cv& operator=(mt_cv const& ref) {
		this.data = ref.data;
		return *this;
	}
};

class log__ {
	std::array<std::string, 4> level_ = {"DBG", "INFO", "WARN", "ERR"};
	std::array<bool, 4> open_{true, true, true, true};
	mt_cv<std::list<std::string>> msgs_;
	std::atomic<bool> write_file_{};
	std::filesystem::path file_dir_{};
	std::fstream file_;
	uint day_{};
	std::mutex file_mtx_;
	std::atomic<bool> run_{true};
	std::thread th_;

	void write(std::string_view msg) {
		if (write_file_) {
			std::lock_guard lk(file_mtx_);
			if (!file_.is_open()) {
				auto&& tm = sys_time();
				std::string name{std::format("log{}-{}-{}.txt", tm.year, tm.mon, tm.mday)};
				auto f = file_dir_ / name;
				file_.open(f, std::ios_base::out | std::ios_base::app);
			}
			file_ << msg << "\n"; //or std::endl to flush
		} else {
			std::cout << msg << std::endl;
		}
	}

public:
	log__() {
		th_  = std::thread([&]() {
			while (run_) {
				std::unique_lock ul(msgs_.mtx);
				while (msgs_.data.empty() && run_) {
					msgs_.cv.wait(ul);
				}
				if (msgs_.data.empty()) {
					ul.unlock();
					continue;
				}
				std::string msg{};
				msg.swap(msgs_.data.front());
				msgs_.data.pop_front();
				ul.unlock();
				write(msg);
			}
			for (const auto& msg : msgs_.data) {
				write(msg);
			}
		});
	}

	~log__() {
		run_ = false;
		{
			std::lock_guard mlg(msgs_.mtx);
			msgs_.cv.notify_one();
		}
		th_.join();
	}

	void write_file(bool b, std::string_view dir = "") {
		write_file_ = b;
		if (b) {
			file_dir_ = dir;
		}
	}

    template <class... Args>
	void log_msg(uint8_t level, std::string_view msgFmt, Args&&... args) {
        auto&& tm = sys_time();
		{
			std::lock_guard lg(file_mtx_);
			if (tm.mday != day_ && write_file_) {
				day_ = tm.mday;
				if (file_.is_open()) file_.close();
			}
		}
		std::string msg{std::format("[{}:{}:{}][{}]", tm.hour, tm.min, tm.sec, level_[level])};
		auto fmtArgs{std::make_format_args(args...)};
		std::string buff{std::vformat(msgFmt, fmtArgs)};
		msg.append(buff);
		std::lock_guard mlg(msgs_.mtx);
		msgs_.data.emplace_back(msg);
		msgs_.cv.notify_one();
    }
};
static log__ glog;
#endif