// fms_binomial.t.cpp - get binomial
#include <cassert>
#include "fms_binomial.h"

using namespace fms;

int binomial_test(int n, double r, double S, double sigma, double k, double t, double tol = 1e-13)
{
	{
		double c = binomial::value(0, 0, n, r, S, sigma, 'C', k, t);
		double p = binomial::value(0, 0, n, r, S, sigma, 'P', k, t);
		double err = (c - p) - (S - k*exp(-r*t));
		assert(fabs(err) < tol);

		double ap = binomial::value(0, 0, n, r, S, sigma, 'P', k, t, true);
		assert(ap >= p);
		double ac = binomial::value(0, 0, n, r, S, sigma, 'C', k, t, true);
		assert(fabs(ac - c) < tol);
	}

	return 0;
}

int binomial_tests()
{
	binomial_test(10, 0, 100, .2, 100, .25);
	binomial_test(10, 0, 100, .1, 90, .25);
	binomial_test(10, 0, 110, .1, 100, .25);

	return 0;
}
int binomial_tests_ = binomial_tests();