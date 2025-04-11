#pragma once
#ifndef THREAD_POOL_H
#define THREAD_POOL_H
#include <tuple>
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool {
public:
	explicit ThreadPool(size_t threads) {
		auto work = [&]() {
			while (true) {
				std::function<void()> task;
				{
					std::unique_lock lock(_mtx);
					_cv.wait(lock, [&] {return !_run || !_tasks.empty();});
					if (!_run && _tasks.empty()) {return;}
					task = std::move(_tasks.front());
					_tasks.pop();
				}
				task();
			}
		};
		_workers.reserve(threads);
		for (size_t i = 0; i < threads; ++i) {
			_workers.emplace_back(work);
		}
	}

	template <class F, class... Args>
	auto add(F&& f, Args&&... args) {
		using return_type = std::invoke_result_t<std::decay_t<F>, std::decay_t<Args>...>;
		auto task = std::make_shared<std::packaged_task<return_type()>>(
			[f = std::forward<F>(f), args = std::make_tuple(std::forward<Args>(args)...)] {
				return std::apply(f, std::move(args));
			}
		);
		auto res = task->get_future();
		{
			std::unique_lock lock(_mtx);
			if (!_run) {throw std::runtime_error("add on stopped ThreadPool");}
			_tasks.emplace([task]() {(*task)();});
		}
		_cv.notify_one();
		return res;
	}

	~ThreadPool() {
		_run = false;
		_cv.notify_all();
		for (auto& worker : _workers) {
			if (worker.joinable()) {
				worker.join();
			}
		}
	}

	ThreadPool(const ThreadPool&) = delete;
	ThreadPool(ThreadPool&&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;
	ThreadPool& operator=(ThreadPool&&) = delete;

private:
	std::vector<std::thread> _workers;
	std::queue<std::function<void()>> _tasks;
	std::mutex _mtx;
	std::condition_variable _cv;
	std::atomic<bool> _run{true};
};

template <class F, class... Args>
inline auto reallyAsync(F&& f, Args&&... params) {
    return std::async(std::launch::async, std::forward<F>(f), std::forward<Args>(params)...);
}

inline void concurrencyF(const auto& fun, size_t size, size_t hc) {
    ThreadPool tp(hc);
    if (size < hc) {
        for (int th = 0; th < static_cast<int>(size); ++th) {
            tp.add(fun, th, th, th + 1);
        }
        return;
    }
    auto len{static_cast<int>(size / hc)};
    auto ts{static_cast<int>(hc)};
    for (int th = 0; th < ts; ++th) {
        auto begin{len * th};
        auto end{begin + len};
        if (th == ts - 1) {
            end = static_cast<int>(size);
        }
        tp.add(fun, th, begin, end);
    }
}
#endif