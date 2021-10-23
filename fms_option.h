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
		enum contract {
			PUT = 'P',
			DIGITAL_PUT = 'Q',
			CALL = 'C',
			DIGITAL_CALL = 'D',
		};

		//  moneyness
		inline double moneyness(double f, double s, double k)
		{
			if (f <= 0 || s <= 0 || k <= 0) {
				return NaN;
			}

			return (log(k / f) + normal::cumulant(s)) / s;
		}

		namespace black {
			// put (k < 0) or call (k > 0) option value
			inline double value(double f, double s, double k)
			{
				if (k < 0) { // put
					double x = moneyness(f, s, -k);

					return (-k) * normal::cdf(x, 0) - f * normal::cdf(x, s);
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

					return -normal::cdf(x, s, 0, 0);
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

				return normal::cdf(x, s, 1, 0) / (f * s);
			}

			// n-th derivative with respect to f
			inline double value(double f, double s, double k, unsigned n)
			{
				if (n == 0) {
					return value(f, s, k);
				}
				if (n == 1) {
					return delta(f, s, k);
				}

				double x = moneyness(f, s, std::fabs(k));

				return normal::cdf(x, s, n - 1, 0) / pow(f * s, n - 1);
			}

			// put (k < 0) or call (k > 0) option vega, dv/ds
			inline double vega(double f, double s, double k)
			{
				double x = moneyness(f, s, fabs(k));

				return -f * normal::cdf(x, s, 0, 1);
			}

			// put (k < 0) or call (k > 0) option theta, -dv/dt
			inline double theta(double f, double sigma, double k, double t, double dt = 1. / 250)
			{
				double s = sigma * sqrt(t);
				double v = value(f, s, k);
				s = sigma * sqrt(t - dt);
				double v_ = value(f, s, k);

				return (v_ - v) / dt;
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
					tol = sqrt(std::numeric_limits<double>::epsilon()); // absolute tolerance
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
		}

		namespace digital {

			inline double value(double f, double s, double k)
			{
				double x = moneyness(f, s, fabs(k));
				double v = normal::cdf(x, 0);

				if (k < 0) {
					return v;
				}
				else if (k > 0) {
					return 1 - v;
				}

				return signbit(k) ? 0 : 1;
			}
			inline double delta(double f, double s, double k)
			{
				double x = moneyness(f, s, fabs(k));
				double v = -normal::cdf(x, 0, 1) / (f * s);

				if (k < 0) {
					return v;
				}
				else if (k > 0) {
					return -v;
				}

				return 0;
			}
			inline double gamma(double f, double s, double k)
			{
				double x = moneyness(f, s, fabs(k));
				double v = (normal::cdf(x, 0, 1) - s * normal::cdf(x, 0, 2)) / (f * f* s);

				if (k < 0) {
					return v;
				}
				else if (k > 0) {
					return -v;
				}

				return 0;
			}
			inline double vega(double f, double s, double k)
			{
				double x = moneyness(f, s, fabs(k));
				double v = normal::cdf(x, 0, 1) * (s * normal::cumulant(s, 1) - x) / s;

				if (k < 0) {
					return v;
				}
				else if (k > 0) {
					return -v;
				}

				return 0;
			}
		}

		namespace bsm {

			// Convert B-S/M parameters to Black forward parameters.
			inline auto Dfs(double r, double S, double sigma, double t)
			{
				double D = exp(-r * t);
				double f = S / D;
				double s = sigma * sqrt(t);

				return std::tuple(D, f, s);
			}

			inline double moneyness(double r, double S, double sigma, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				return option::moneyness(f, s, fabs(k));
			}

			inline double value(double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return D * black::value(f, s, -k);
				case option::contract::CALL:
					return D * black::value(f, s, k);
				case option::contract::DIGITAL_PUT:
					return D * digital::value(f, s, -k);
				case option::contract::DIGITAL_CALL:
					return D * digital::value(f, s, k);
				}

				return NaN;
			}

			// delta
			inline double delta(double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return black::delta(f, s, -k);
				case option::contract::CALL:
					return black::delta(f, s, k);
				case option::contract::DIGITAL_PUT:
					return digital::delta(f, s, -k);
				case option::contract::DIGITAL_CALL:
					return digital::delta(f, s, k);
				}

				return NaN;
			}
			// gamma
			inline double gamma(double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return black::gamma(f, s, -k) / D;
				case option::contract::CALL:
					return black::gamma(f, s, k) / D;
				case option::contract::DIGITAL_PUT:
					return digital::gamma(f, s, -k) / D;
				case option::contract::DIGITAL_CALL:
					return digital::gamma(f, s, k) / D;
				}

				return NaN;
			}
			// vega
			inline double vega(double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return black::vega(f, s, -k) * D * sqrt(t);
				case option::contract::CALL:
					return black::vega(f, s, k) * D * sqrt(t);
				case option::contract::DIGITAL_PUT:
					return digital::vega(f, s, -k) * D * sqrt(t);
				case option::contract::DIGITAL_CALL:
					return digital::vega(f, s, k) * D * sqrt(t);
				}

				return NaN;
			}

		} // namespace bsm
	}
}