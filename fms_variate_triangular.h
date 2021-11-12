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
		// P^s(X <= x) and derivatives
		double _cdf(double x, double s, unsigned nx = 0, unsigned ns = 0) const override
		{
			x = x; s = s; nx = nx, ns = ns;

			//!!! implement for nx = 0, ns = 0
			//!!! implement for nx = 1, ns = 0
			//!!! implement for nx = 1, ns = 1

			return std::numeric_limits<double>::quiet_NaN();
		}

		// kappa(s) = log E[e^{s X}] and derivativs
		double _cumulant(double s, unsigned n = 0) const override
		{
			s = s; n = n;

			// !!! implement for n = 0

			return std::numeric_limits<double>::quiet_NaN();
		}
	};

} // namespace fms::variate
