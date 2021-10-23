// fms_pwflat.h - Piecewise flat left-continuous curve
#pragma once
#ifdef _DEBUG
#include <cassert>
#endif
#include <algorithm>
#include <limits>

namespace fms::pwflat {

	template<class X>
	using NaN = std::numeric_limits<X>::quiet_NaN();

	// f(u) = f[i], t[i-1] < u <= t[i], 0 <= i < n
	// f(u) = NaN, u < 0
	// f(u) = _f, u > t[n-1]
	template<class T, class F>
	inline F value(T u, unsigned n, const T* t, const F* f, F _f = NaN<F>)
	{
		if (u < 0) {
			return NaN<F>;
		}

		auto i = std::lower_bound(t, t + n, u);
		if (i == t + n) {
			return _f;
		}

		return f[i - t];
	}
#ifdef _DEBUG
	inline int value_test()
	{
		double t[] = { 1,2,3 };
		double f[] = { .1,.2,.3 };
		assert(_isnan(value(-1., 3, t, f)));

	}
#endif // _DEBUG

} // namespace fms
