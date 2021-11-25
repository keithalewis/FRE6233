// fms_variate_normal.t.cpp - Test fms::variate::normal
#ifdef _DEBUG
// Only test in debug mode
#include <cassert>
#include <algorithm>
#include <random>
#include "fms_derivative.h"
#include "fms_monte_carlo.h"
#include "fms_option.h"
#include "fms_variate_normal.h"

using namespace fms;
using namespace fms::option;

// s_n = (x_1 + ... + x_n)/n
// n s_n - (n-1) s_{n-1} = x_n
// s_n = s_{n-1} + (x_n - s_{n-1})/n
double monte_carlo_option_value(double f, double s, double k, size_t n = 10000)
{
	std::function<double(double)> payoff;
	if (k < 0) {
		payoff = [k](double x) { return std::max(-k - x, 0.); };
	}
	else {
		payoff = [k](double x) { return std::max(x - k, 0.); };
	}

	std::default_random_engine dre;
	std::normal_distribution<double> X;
	auto p = [f, s, &payoff, &X, &dre]() {
		double F = f * exp(s * X(dre) - s * s / 2);

		return payoff(F);
	};

	return monte_carlo::average(n, p);
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
	// for ks !!! test both k and -k
	// for ss
	// for is
	/* !!!add for loops above and fix up below
	double var = option::black::variance(f, s, k);
	double f = 100, s = 0.1, k = 100;
	int n = 10000;
	double v = option::black::value(N, f, s, k);
	double vn = monte_carlo_option_value(f, s, k, n);
	double sd = 2;
	assert(fabs(v - vn) <= var * sd / sqrt(n));
	*/

	return 0;
}
int option_value_test_ = option_value_test();

int option_delta_test_ = 0;
int option_gamma_test_ = 0;
int option_vega_test_ = 0;
int option_implied_test_ = 0;

#endif // _DEBUG