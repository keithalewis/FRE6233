// fms_variate_normal.t.cpp - Test fms::variate::normal
// Only test in debug mode
#ifdef _DEBUG
#include <cassert>
#include <algorithm>
#include "fms_variate_normal.h"
#include "fms_derivative.h"

using namespace fms;
using namespace fms::variate;

int fms_variant_normal_H_test()
{
	double xs[] = { -1, 0, 1, 2 };
	{
		// H_0(x) = 1
		for (double x : xs) {
			assert(1 == normal::H(0, x));
		}

		// H_1(x) = x
		for (double x : xs) {
			assert(x == normal::H(1, x));
		}

		// H_2(x) = x^2 - 1
		for (double x : xs) {
			assert(x * x - 1 == normal::H(2, x));
		}
	}

	return 0;
}
// cause the test to be run when the dll is loaded
int fms_variant_normal_H_test_ = fms_variant_normal_H_test();

// test N^{(n)} at x
template<class X = double, class Y = double>
inline bool normal_derivative_test(int n, X x, X h)
{
	Y df = normal::N(x, n + 1);
	Y dddf = normal::N(x, n + 3);
	auto f = [n](double x) { return normal::N(x, n);  };

	return derivative_test<X,Y>(f, x, h, df, dddf);
}

int fms_variant_normal_N_test()
{
	{
		// sanity checks
		assert(0.5 == normal::N(0));
		assert(1 / M_SQRT2PI == normal::N(0, 1));
		assert(0 == normal::N(0, 2));
	}
	{
		double xs[] = { -1, 0, 1, 2 };
		double hs[] = { 0.1, 0.01, 0.001, 0.0001 };
		for (int n : { 0, 1, 2 }) {
			for (double x : xs) {
				for (double h : hs) {
					assert(normal_derivative_test(n, x, h));
				}
			}
		}
	}

	return 0;
}
int fms_variant_normal_N_test_ = fms_variant_normal_N_test();
#if 0

int fms_variant_normal_cdf_test()
{
	for (double x : sequence(-2, 2, .1)) {
		for (double s : sequence(-1, 1, .1)) {
			for (int nx : {0, 1, 2, 3}) {
				for (int ns : {0, 1, 2, 3}) {
					auto fx = [s, nx, ns](double x) { return normal::cdf(x, s, nx, ns); };
					auto fs = [x, nx, ns](double s) { return normal::cdf(x, s, nx, ns); };
					for (double h : { .01, .001, .0001, .00001}) {
						{
							auto df = normal::cdf(x, s, nx + 1, ns);
							double dddf = normal::cdf(x, s, nx + 3, ns);
							assert((derivative_test<double,double>(fx, x, h, df, dddf, 10.)));
						}
						{
							auto df = normal::cdf(x, s, nx, ns + 1);
							double dddf = normal::cdf(x, s, nx, ns + 3);
							assert((derivative_test<double, double>(fs, s, h, df, dddf, 20.)));
						}
					}
				}
			}
		}
	}

	return 0;
}
//int fms_variant_normal_cdf_test_ = fms_variant_normal_cdf_test();

#endif // _DEBUG
#endif // 0