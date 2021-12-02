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

double monte_carlo_option_variance(double f, double s, double k, size_t n = 10000)
{
	std::function<double(double)> payoff;
	if (k < 0) {
		payoff = [k](double x) { return std::pow(std::max(-k - x, 0.),2); };
	}
	else {
		payoff = [k](double x) { return std::pow(std::max(x - k, 0.),2); };
	}

	std::default_random_engine dre;
	std::normal_distribution<double> X;
	auto p = [f, s, &payoff, &X, &dre]() {
		double F = f * exp(s * X(dre) - s * s / 2);

		return payoff(F);
	};

	return monte_carlo::average(n, p) - std::pow(monte_carlo_option_value(f, s, k, n), 2);
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
	double f = 100, s = 0.1, k = 100;
	double stdev = sqrt(option::black::variance(N, f, s, k));
	int n = 10000;
	double v = option::black::value(N, f, s, k);
	double vn = monte_carlo_option_value(f, s, k, n);
	double sd = 2;
	assert(fabs(v - vn) <= stdev * sd / sqrt(n));
	*/

	return 0;
}

int option_variance_test()
{
	double sd = 2; // two standard deviations
	for (int f_num = 0; f_num < 5; f_num++) {
		double f = fs[f_num];
		for (int k_num = 0; k_num < 10; k_num++) {
			double k;
			if (k_num < 5) k = ks[k_num];
			else k = -ks[k_num - 5];
			for (int s_num = 0; s_num < 4; s_num++) {
				double s = ss[s_num];
				double stdev = sqrt(option::black:：moment4(N, f, s, k)- std::pow(option::black::variance(N, f, s, k),2));
				int n = is[0];
				double v = option::black::variance(N, f, s, k);
				double vn = monte_carlo_option_variance(f, s, k, n);
				double sd = 2;
				assert(fabs(v - vn) <= stdev * sd / sqrt(n));

			}
		}
	}

	return 0;
}

int option_value_test_ = option_value_test();

int option_delta_test_ = 0;
int option_gamma_test_ = 0;
int option_vega_test_ = 0;
int option_implied_test_ = 0;
int option_variance_test_ = option_variance_test();

#endif // _DEBUG
