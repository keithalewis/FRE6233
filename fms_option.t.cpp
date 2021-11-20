// fms_variate_normal.t.cpp - Test fms::variate::normal
#ifdef _DEBUG
// Only test in debug mode
#include <cassert>
#include <algorithm>
#include <random>
#include "fms_option.h"
#include "fms_variate_normal.h"
#include "fms_derivative.h"

using namespace fms;
using namespace fms::option;

// s_n = (x_1 + ... + x_n)/n
// n s_n - (n-1) s_{n-1} = x_n
// s_n = s_{n-1} + (x_n - s_{n-1})/n
double monte_carlo_option_value(double f, double s, double k, int n = 10000)
{
	double v = 0;
	std::default_random_engine dre;
	std::normal_distribution<double> N;

	if (k < 0) {
		k = -k;
		for (int i = 1; i <= n; ++i) {
			double F = f * exp(s * N(dre) - s * s / 2);
			v += (std::max(k - F, 0.) - v) / i;
		}
	}
	else {
		for (int i = 1; i <= n; ++i) {
			double F = f * exp(s * N(dre) - s * s / 2);
			v += (std::max(F - k, 0.) - v) / i;
		}
	}

	return v;
}

// common to all tests
variate::normal N;
double fs[] = { 80, 90, 100, 110, 120 };
double ks[] = { 80, 90, 100, 110, 120 };
double ss[] = { .01, .02, .1, .2 };
int is[] = { 10000 };

int option_value_test()
{
	// double sd = 2; // two standard deviations
	// for fs
	// for ks
	// for ss
	// for is
	/* !!!add for loops above and fix up below
	double f = 100, s = 0.1, k = 100;
	int n = 10000;
	double v = option::black::value(N, f, s, k);
	double vn = monte_carlo_option_value(f, s, k, n);
	double sd = 2;
	assert(fabs(v - vn) <= v * sd / sqrt(n));
	*/

	return 0;
}
int option_value_test_ = option_value_test();

int option_delta_test_ = 0;
int option_gamma_test_ = 0;
int option_vega_test_ = 0;
int option_implied_test_ = 0;

#endif // _DEBUG