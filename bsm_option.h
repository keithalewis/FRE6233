// bsm_option.h - Black-Scholes/Merton option value and greeks
#pragma once
#include <cmath>
#include <limits>
#include "xll/xll/xll.h"

#ifndef CATEGORY
#define CATEGORY "FRE6233"
#endif

namespace bsm {

	// Return NaN to indicate error.
	constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

	// Machine epsilon
	constexpr double epsilon = std::numeric_limits<double>::epsilon();

	// sqrt(2 pi)
	constexpr double M_SQRT2PI = 2.50662827463100050240;
#ifndef M_SQRT2
	// sqrt(2)
	constexpr double M_SQRT2 = 1.41421356237309504880;
#endif

	// standard normal random variate
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

		//  moneyness
		inline double moneyness(double f, double s, double k)
		{
			if (f <= 0 || s <= 0 || k <= 0) {
				return NaN;
			}

			return (log(k / f) + normal::cumulant(s)) / s;
		}

		// put (k < 0) or call (k > 0) option value
		inline double value(double f, double s, double k)
		{
			if (k < 0) { // put
				double m = moneyness(f, s, -k);

				return (-k) * normal::cdf(m) - f * normal::cdf(m, s);
			}
			else { // call
				// c = p + f - k
				return option::value(f, s, -k) + f - k;
			}
		}

		// put (k < 0) or call (k > 0) option delta, dv/dk
		inline double delta(double f, double s, double k)
		{
			if (k < 0) { // put
				double m = moneyness(f, s, -k);

				return -normal::cdf(m, s);
			}
			else { // call
				// dc/df = dp/df + 1
				return option::delta(f, s, -k) + 1;
			}
		}

		// put (k < 0) or call (k > 0) option vega, dv/ds
		inline double vega(double f, double s, double k)
		{
			k = std::fabs(k); // same for put or call
			double m = moneyness(f, s, k);

			return -f * normal::cdf(m, s, 0, 1);
		}

		// implied volatility using initial guess, max number of iterations, and tolerance
		inline double implied(double f, double v, double k,
			double s = 0, unsigned n = 0, double eps = 0)
		{
			if (s == 0) {
				s = 0.1; // initial vol guess
			}
			if (n == 0) {
				n = 100; // maximum number of iterations
			}
			if (eps == 0) {
				eps = sqrt(epsilon); // absolute tolerance
			}

			double v_ = value(f, s, k);
			double dv_ = vega(f, s, k);
			double s_ = s - (v_ - v) / dv_; // Newton-Raphson
			if (s_ < 0) {
				s_ = s / 2;
			}
			while (fabs(s_ - s) > eps) {
				v_ = value(f, s_, k);
				dv_ = vega(f, s_, k);
				s = s_ - (v_ - v) / dv_;
				if (s < 0) {
					s = s_ / 2;
				}
				std::swap(s_, s);
				if (n == 0) {
					return NaN;
				}
				--n;
			}

			return s_;
		}
	}

}