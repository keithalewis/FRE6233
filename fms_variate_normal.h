// fms_variate_normal.h - Normally distributed random variate
#pragma once
#include <cmath>

namespace fms::variate {

	// sqrt(2 pi)
	constexpr double M_SQRT2PI = 2.50662827463100050240;
#ifndef M_SQRT2
	// sqrt(2)
	constexpr double M_SQRT2 = 1.41421356237309504880;
#endif

	// standard normal random variate
	struct normal {

		// Hermite polynomials
		// H_{n+1}(x) = x H_n(x) - n H_{n-1}(x), H_0(x) = 1, H_1(x) = x
		static constexpr double H(unsigned n, double x) noexcept
		{
			if (n == 0) {
				return 1;
			}
			if (n == 1) {
				return x;
			}

			return x * H(n - 1, x) - (n - 1) * H(n - 2, x);
		}

		// P(X <= x) and derivatives
		static double N(double x, unsigned n = 0)
		{
			if (n == 0) {
				return  (1 + erf(x / M_SQRT2)) / 2;
			}

			double phi = exp(-x * x / 2) / M_SQRT2PI;

			if (n == 1) {
				return phi;
			}

			// (d/dx)^n phi(x) = (-1)^n phi(x) H_n(x)
			return phi * H(n - 1, x) * ((n & 1) ? 1 : -1);
		}

		// P_t^sigma(X_t <= x) = N((x - sigma t)/sqrt(t))
		// D_x^n N = N^{(n)}(...) / (sqrt(t))^n
		// D_sigma^n N = N^{(n)}(...) * (sqrt(t))^n (-1)^n
		static double cdf(double t, double x, double sigma, unsigned nx = 0, unsigned nsigma = 0)
		{
			double srt = sqrt(t);

			return N(x/srt - sigma * srt, nx + nsigma) 
				* pow(srt, nsigma - nx) * (nsigma & 1 ? -1 : 1);
		}

		// kappa(t, sigma) = log E[e^{sigma X_t}] = sigma^2 t/2 and derivativs
		static double cumulant(double t, double sigma, unsigned n = 0)
		{
			if (n == 0) {
				return sigma * sigma * t/ 2;
			}
			if (n == 1) {
				return sigma * t;
			}
			if (n == 2) {
				return t;
			}

			return 0;
		}
	};

} // namespace fms::variate
