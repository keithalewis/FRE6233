// xll_template.h - common includes
#pragma once
#include <cmath> // erf
#include <limits>
#include "xll/xll/xll.h"

#ifndef CATEGORY
#define CATEGORY "FRE6233"
#endif

namespace FRE6233 {

	// Return NaN to indicate error.
	constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
	// sqrt(2 pi)
	constexpr double M_SQRT2PI = 2.50662827463100050240;
	//#ifndef M_SQRT2
	// sqrt(2)
	constexpr double M_SQRT2 = 1.41421356237309504880;
	//#endif

	namespace normal {

		// P(X <= x) and derivatives
		inline double cdf(double x, int n = 0) 
		{
			if (n == 0) {
				return  (1 + ::erf(x / M_SQRT2)) / 2;
			}
			double phi = ::exp(-x * x / 2) / M_SQRT2PI;
			if (n == 1) {
				return phi;
			}
			//!!! Hermite polynomials

			return NaN;
		}
		// P_s(X <= x) = P(X <= x - s)
		inline double cdf(double x, double s, int nx = 0, int ns = 0)
		{
			return cdf(x - s, nx + ns) * (ns % 2 ? -1 : 1);
		}

		// kappa(s) = log E[e^{sX}] = s^2/2 and derivativs
		inline double cumulant(double s, int n = 0)
		{
			if (n == 0) {
				return s * s / 2;
			}
			if (n == 1) {
				return s;
			}
			if (n == 2) {
				return 1;
			}

			return 0;
		}
	}

	namespace option {
		inline double moneyness(double f, double s, double k)
		{
			if (f <= 0 || s <= 0 || k <= 0) {
				return NaN;
			}

			return (log(k / f) + normal::cumulant(s)) / s;
		}

		// put option value
		inline double value(double f, double s, double k)
		{
			double m = moneyness(f, s, k);

			return k * normal::cdf(m) - f * normal::cdf(m, s);
		}

		// put option delta
		inline double delta(double f, double s, double k)
		{
			double m = moneyness(f, s, k);

			return -normal::cdf(m, s);
		}

	}

}