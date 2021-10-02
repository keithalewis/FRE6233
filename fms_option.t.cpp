// fms_variate_normal.t.cpp - Test fms::variate::normal
// Only test in debug mode
#ifdef _DEBUG
#include <cassert>
#include <algorithm>
#include <numeric>
#include <valarray>
#include "fms_option.h"
#include "fms_derivative.h"

using namespace fms;
using namespace fms::option;

inline std::valarray<double> sequence(double b, double e, double dx)
{
	size_t n = static_cast<unsigned>(1 + (e - b) / dx);
	std::valarray<double> a(n);

	for (size_t i = 0; i < n; ++i) {
		a[i] = b + i * dx;
	}

	return a;
}

int fms_option_test()
{
	for (auto f : sequence(50, 100, 10)) {
		for (auto sigma : sequence(0.01, 1, 0.01)) {
			for (auto k : sequence(50, 100, 10)) {
				for (auto t : sequence(0.1, 2, 0.1)) {
					for (double h : {.01, .001, .0001, .00001}) {
						{ // test delta
							auto v = [sigma, k, t](double f) { return option::value(f, sigma, k, t); };
							double dv = option::delta(f, sigma, k, t);
							derivative_test<double, double>(v, f, h, dv, 1., 1.);
						}
						{ // test gamma
							auto v = [sigma, k, t](double f) { return option::delta(f, sigma, k, t); };
							double dv = option::gamma(f, sigma, k, t);
							derivative_test<double, double>(v, f, h, dv, 1., 1.);
						}
						{ // test vega
							auto v = [f, k, t](double sigma) { return option::value(f, sigma, k, t); };
							double dv = option::vega(f, sigma, k, t);
							derivative_test<double, double>(v, sigma, h, dv, 1., 1.);
						}
					}
				}
			}
		}
	}

	return 0;
}
int fms_option_test_ = fms_option_test();

#endif // _DEBUG