// fms_binomial.h - binomial model
#pragma once
#include <cmath>
#include <functional>

// indicate error
#define ensure(e) if (!(e)) { return std::numeric_limits<double>::quiet_NaN(); }

namespace fms::binomial {

	// e^{-rt(n - j)/n}E[nu(S_n) | W_j = i], V_j = j - 2 W_j
	// S_n = S e^{rt + sigma sqrt(t)/sqrt(n) V_n}/cosh(sigma sqrt(t)/sqrt(n))^n
	inline double value(int i, int j, int n, double r, double S, double sigma, double t, const std::function<double(double)>& nu, bool american = false)
	{
		double s = sigma * sqrt(t) / sqrt(n);
		double S_j = S * exp(r * t * j / n + s * (j - 2. * i)) / pow(cosh(s), j);

		if (j == n) {
			return nu(S_j);
		}

		double v = value(i, j + 1, n, r, S * exp(r * t / n), sigma, t - 1./n, nu, american) / 2
			+ value(i + 1, j + 1, n, r, S * exp(r * t / n), sigma, t - 1./n, nu, american) / 2;

		if (american) {
			//!!! optimal exercise code goes here !!!
		}

		return v;
	}

	// American put (k < 0) or call (p > 0) value at time j given W_j = i
	inline double value(int i, int j, int n, double r, double S, double sigma, double k, double t, bool american = false)
	{
		std::function<double(double)> nu;

		if (k < 0) { // put
			nu = [k](double F) { return std::max(-k - F, 0.); };
		}
		else if (k > 0) { // call
			nu = [k](double F) { return std::max(F - k, 0.); };
		}
		else {
			return signbit(k) ? 0 : f;
		}

		return value(i, j, n, f, s, nu, american);
	}

	// American B-S/M put or call value at time j given W_j = i, where S_j = e^{-rtj/n} F_j
	inline double value(int i, int j, int n, double r, double S, double sigma, int c, double k, double t, bool american = false)
	{
		double Dn = exp(-r * t * (n - j) / n);
		double Dj = exp(-r * t * j / n);
		double s = sigma * sqrt(t);

		if (c == 'P') {
			return Dn * value(i, j, n, S / Dj, s, -k, american);
		}
		else if (c == 'C') {
			return Dn * value(i, j, n, S / Dj, s, k, american);
		}

		return std::numeric_limits<double>::quiet_NaN();
	}

} // namespace fms

#undef ensure