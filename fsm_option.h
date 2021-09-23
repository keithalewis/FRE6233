// fsm_option.h - Option value and greeks
#pragma once
#include <cmath>
#include <limits>
#include <tuple>
#include "fms_variate_normal.h"

using namespace fms::variate;

namespace fms {

	// Return NaN to indicate error.
	constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

	// Machine epsilon
	constexpr double epsilon = std::numeric_limits<double>::epsilon();

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
				double x = moneyness(f, s, -k);

				return (-k) * normal::cdf(x) - f * normal::cdf(x, s);
			}
			else if (k > 0) { // call
				// c = p + f - k
				return value(f, s, -k) + f - k;
			}

			// k = -/+ 0
			return signbit(k) ? 0 : f;
		}

		// put (k < 0) or call (k > 0) option delta, dv/df
		inline double delta(double f, double s, double k)
		{
			if (k < 0) { // put
				double x = moneyness(f, s, -k);

				return -normal::cdf(x, s);
			}
			else if (k > 0) { // call
				// dc/df = dp/df + 1
				return delta(f, s, -k) + 1;
			}

			return signbit(k) ? 0 : 1;
		}

		// put (k < 0) or call (k > 0) option gamma, d^2v/df^2
		inline double gamma(double f, double s, double k)
		{
			double x = moneyness(f, s, std::fabs(k));

			return normal::cdf(x, s, 1) / (f * s);
		}

		// put (k < 0) or call (k > 0) option vega, dv/ds
		inline double vega(double f, double s, double k)
		{
			double x = moneyness(f, s, fabs(k));

			return -normal::cdf(x, s, 0, 1) * f;
		}

		// put (k < 0) or call (k > 0) option theta, dv/dt
		inline double theta(double f, double sigma, double k, double t)
		{
			return vega(f, sigma * sqrt(t), k) * sigma / (2 * sqrt(t));
		}

		// implied volatility using initial guess, max number of iterations, and tolerance
		inline double implied(double f, double v, double k,
			double s = 0, unsigned n = 0, double tol = 0)
		{
			// max(k - f,0) >= k - f
			// max(k - f,0) <= k
			if (k < 0) {
				if (v <= std::max(-k - f, 0.) || v >= k) {
					return NaN;
				}
			}
			// max(f - k,0) >= f - k
			// max(f - k,0) <= f
			else if (k > 0) {
				if (v <= std::max(f - k, 0.) || v >= f) {
					return NaN;
				}
			}

			if (s == 0) {
				s = 0.1; // initial vol guess
			}
			if (n == 0) {
				n = 100; // maximum number of iterations
			}
			if (tol == 0) {
				tol = sqrt(epsilon); // absolute tolerance
			}

			double v_ = value(f, s, k);
			double dv_ = vega(f, s, k); // dv/ds
			double s_ = s - (v_ - v) / dv_; // Newton-Raphson
			if (s_ < 0) {
				s_ = s / 2;
			}
			while (fabs(s_ - s) > tol) {
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

		struct contract {
			double k; // strike
			double t; // expiration
		};
		// different types with the same data
		struct put : contract {};
		struct call : contract {};
		struct digital_put : contract {};
		struct digital_call : contract {};

		// Black-Scholes/Mertion option value and greeks
		namespace bsm {

			// Convert B-S/M parameters to Black forward parameters.
			inline auto Dfsk(double r, double S, double sigma, const contract& o)
			{
				double D = exp(-r * o.t);
				double f = S / D;
				double s = sigma * sqrt(o.t);

				return std::tuple(D, f, s, o.k);
			}

			// call using moneyness(r, S, sigma, contract({k, t}))
			// or moneyness(r, S, sigma, (contract){.k = k, .t = t})
			inline double moneyness(double r, double S, double sigma, const contract& o)
			{
				auto [D, f, s, k] = Dfsk(r, S, sigma, o);

				return option::moneyness(f, s, fabs(o.k));
			}

			// call using value(r, S, sigma, put({k, t}))
			inline double value(double r, double S, double sigma, put o)
			{
				auto [D, f, s, k] = Dfsk(r, S, sigma, o);

				return D * option::value(f, s, -o.k);
			}

			inline double value(double r, double S, double sigma, call o)
			{
				auto [D, f, s, k] = Dfsk(r, S, sigma, o);

				return D * option::value(f, s, o.k);
			}

			// delta, ...

		} // namespace bsm
	}

}