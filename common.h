#ifndef COMMON_H
#define COMMON_H
#include <iostream>
#include <random>
#include "json/json.hpp"
#include <chrono>

#define INIT_CONCENTRATION 0.5

static std::default_random_engine grnRand = std::default_random_engine(
    std::chrono::system_clock::now().time_since_epoch().count());
enum class ProteinType { input = 0, regul = 1, output = 2 };
template <typename T> T mix(const T& a, const T& b, const double v) {
	double r = v > 1.0 ? 1.0 : (v < 0.0 ? 0.0 : v);
	return (a * (1.0 - r)) + (r * b);
}

template <typename E>
inline static constexpr typename std::underlying_type<E>::type to_underlying(E e) {
	return static_cast<typename std::underlying_type<E>::type>(e);
}

template <typename P, typename C>
inline void eachP(const std::initializer_list<ProteinType>& l, C& container,
                  const std::function<void(P&)>& f) {
	for (const auto& ptype : l) {
		const auto t = to_underlying(ptype);
		for (auto& p : container[t]) {
			f(p.second);
		}
	}
}

template <typename P, typename C>
inline void eachP(const std::initializer_list<ProteinType>& l, C& container,
                  const std::function<void(P&, size_t pt)>& f) {
	for (const auto& ptype : l) {
		const auto t = to_underlying(ptype);
		for (auto& p : container[t]) {
			f(p.second, t);
		}
	}
}

template <typename I, typename T, unsigned int maxn> struct stackUmap {
	std::array<I, maxn> indices;
	std::array<T, maxn> values;
	size_t size;
};

#endif
