// fms_variate_normal.t.cpp - Test fms::variate::normal
// Only test in debug mode
#ifdef _DEBUG
#include <cassert>
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
// force the test to be run when the dll is loaded
int fms_variant_normal_H_test_ = fms_variant_normal_H_test();

int fms_variant_normal_N_test()
{
	{
		// sanity checks
		assert(0.5 == normal::N(0));
		assert(1 / M_SQRT2PI == normal::N(0, 1));
		assert(0 == normal::N(0, 2));
	}
	{
		// N'(x)
		double xs[] = { -1, 0, 1, 2 };
		double hs[] = { 0.1, 0.01, 0.001, 0.0001 };
		for (double x : xs) {
			for (double h : hs) {
				double df = (normal::N(x + h) - normal::N(x - h)) / (2 * h);
				assert(fabs(normal::N(x, 1) - df) < h * h);
			}
		}
		// N''(x)
		for (double x : xs) {
			for (double h : hs) {
				double df = (normal::N(x + h, 1) - normal::N(x - h, 1)) / (2 * h);
				assert(fabs(normal::N(x, 2) - df) < h * h);
			}
		}
		// N^(k)(x)
		unsigned ks[] = { 2, 3, 4 };
		for (unsigned k : ks) {
			for (double x : xs) {
				for (double h : hs) {
					double df = (normal::N(x + h, k) - normal::N(x - h, k)) / (2 * h);
					assert(fabs(normal::N(x, k + 1) - df) < h * h);
				}
			}
		}
	}

	return 0;
}
int fms_variant_normal_N_test_ = fms_variant_normal_N_test();

int fms_variant_normal_cdf_test()
{
	for (double x : sequence(-2, 2, .1)) {
		for (double s : sequence(-1, 1, .1)) {
			for (int nx : {0, 1, 2, 3}) {
				for (int ns : {0, 1, 2, 3}) {
					for (double t : sequence(0.1, 2, .1)) {
						for (double h : { .01, .001, .0001, .00001}) {
							{
								auto f = [t, s, nx, ns](double x) { return normal::cdf(t, x, s, nx, ns); };
								auto df = normal::cdf(t, x, s, nx + 1, ns);
								double dddf = normal::cdf(t, x, s, nx + 3, ns);
								assert((derivative_test<double,double>(f, x, h, df, dddf, 10.)));
							}
							{
								auto f = [t, x, nx, ns](double s) { return normal::cdf(t, x, s, nx, ns); };
								auto df = normal::cdf(t, x, s, nx, ns + 1);
								double dddf = normal::cdf(t, x, s, nx, ns + 3);
								assert((derivative_test<double, double>(f, s, h, df, dddf, 10.)));
							}
						}
					}
				}
			}
		}
	}

	return 0;
}
int fms_variant_normal_cdf_test_ = fms_variant_normal_cdf_test();

#endif // _DEBUG