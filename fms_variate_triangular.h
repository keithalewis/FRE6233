// fms_variate_triangular.h - Triangular variate
// A triangular distribution is defined by three values:
// l, m, h. Its density is a linear function connecting
// (l, 0), (m, 2/(h - l)), and (h, 0).
//
// f(x) = 0, x < l;
// f(x) = 2(x - l)/(m - l)(h - l), l <= x <= m;
// f(x) = 2(h - x)/(h - m)(h - l), m <= x <= h;
// f(x) = 0, x > h;
#pragma once
#include <limits>
#include "fms_variate.h"

namespace fms::variate {

	// triangular distribution
	struct triangular : public fms::variate::base {
		double l, m, h;
		triangular(double l, double m, double h)
			: l(l), m(m), h(h)
		{ }
		// P^s(X <= x) = E[e^{s X - kappa(s)} 1(X <= x)] and derivatives
		double _cdf(double x, double s, unsigned nx = 0, unsigned ns = 0) const override
		{
			x = x; s = s; nx = nx, ns = ns;

			//!!! implement for nx = 0, ns = 0
			//!!! implement for nx = 1, ns = 0
			//!!! implement for nx = 1, ns = 1

			return std::numeric_limits<double>::quiet_NaN();
		}

		// kappa(s) = log E[e^{s X}] and derivatives
		// E[e^{s X}] =
		//     int_l^m e^{sx} a(x - l) dx
		//   + int_m^h e^{sx} b(h - x) dx
		//
		// int e^{sx} dx = e^{sx}/s
		// int x e^{sx} dx = e^{sx}(x/s - 1/s^2)
		//
		// E[e^{s X}] =
		//   [a e^{sx}(x/s - 1/s^2) - al e^{sx}/s]_l^m
		// + [bh e^{sx}/s - b e^{sx}(x/s - 1/s^2)]_m^h
		//
		double mgf(double s) const
		{
			double a = 2 / ((m - l) * (h - l));
			double b = 2 / ((h - m) * (h - l));

			auto I = [s](double x) {
				return exp(s * x) / s;
			};
			auto Ix = [s](double x) {
				return exp(s * x) * (x / s - 1 / (s * s));
			};

			double Esx = (a * Ix(m) - a * l * I(m)) - (a * Ix(l) - a * l * I(l));
			Esx += (b * h * I(h) - b * Ix(h)) - (b * h * I(m) - b * Ix(m));

			return Esx;
		}
		double _cumulant(double s, unsigned n = 0) const override
		{
			if (n != 0) {
				return std::numeric_limits<double>::quiet_NaN();
			}

			return log(mgf(s));
		}
	};

} // namespace fms::variate
