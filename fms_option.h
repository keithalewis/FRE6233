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

	namespace option {

		//  moneyness
		inline double moneyness(double f, double sigma, double k, double t)
		{
			if (f <= 0 || sigma <= 0 || k <= 0 || t <= 0) {
				return NaN;
			}

			return (log(k / f) + normal::cumulant(t, sigma)) / sigma;
		}

		// put (k < 0) or call (k > 0) option value
		inline double value(double f, double sigma, double k, double t)
		{
			if (k < 0) { // put
				double x = moneyness(f, sigma, -k, t);

				return (-k) * normal::cdf(t, x, 0) - f * normal::cdf(t, x, sigma);
			}
			else if (k > 0) { // call
				// c = p + f - k
				return value(f, sigma, -k, t) + f - k;
			}

			// k = -/+ 0
			return signbit(k) ? 0 : f;
		}

		// put (k < 0) or call (k > 0) option delta, dv/df
		inline double delta(double f, double sigma, double k, double t)
		{
			if (k < 0) { // put
				double x = moneyness(f, sigma, -k, t);

				return -normal::cdf(t, x, sigma, 0, 0);
			}
			else if (k > 0) { // call
				// dc/df = dp/df + 1
				return delta(f, sigma, -k, t) + 1;
			}

			return signbit(k) ? 0 : 1;
		}

		// put (k < 0) or call (k > 0) option gamma, d^2v/df^2
		inline double gamma(double f, double sigma, double k, double t)
		{
			double x = moneyness(f, sigma, std::fabs(k), t);

			return normal::cdf(t, x, sigma, 1, 0) / (f * sigma);
		}

		// put (k < 0) or call (k > 0) option vega, dv/ds
		inline double vega(double f, double sigma, double k, double t)
		{
			double x = moneyness(f, sigma, fabs(k), t);

			return -f * normal::cdf(x, sigma, 0, 1);
		}

		// put (k < 0) or call (k > 0) option theta, -dv/dt
		inline double theta(double f, double sigma, double k, double t, double dt = 1./250)
		{
			double v = value(f, sigma, k, t);
			double v_ = value(f, sigma, k, t - dt);

			return (v_ - v) / dt;
		}

		// implied volatility using initial guess, max number of iterations, and tolerance
		inline double implied(double f, double v, double k, double t,
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
				tol = sqrt(std::numeric_limits<double>::epsilon()); // absolute tolerance
			}

			double v_ = value(f, s, k, t);
			double dv_ = vega(f, s, k, t); // dv/ds
			double s_ = s - (v_ - v) / dv_; // Newton-Raphson
			if (s_ < 0) {
				s_ = s / 2;
			}
			while (fabs(s_ - s) > tol) {
				v_ = value(f, s_, k, t);
				dv_ = vega(f, s_, k, t);
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

			return option::moneyness(f, s, fabs(o.k), o.t);
		}

		// call using value(r, S, sigma, put({k, t}))
		inline double value(double r, double S, double sigma, put o)
		{
			auto [D, f, s, k] = Dfsk(r, S, sigma, o);

			return D * option::value(f, s, -o.k, o.t);
		}
		inline double value(double r, double S, double sigma, call o)
		{
			auto [D, f, s, k] = Dfsk(r, S, sigma, o);

			return D * option::value(f, s, o.k, o.t);
		}

		// delta, ...

	} // namespace bsm

}