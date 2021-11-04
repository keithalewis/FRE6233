// fms_pwflat.h - Piecewise flat left-continuous curve
#pragma once
#ifdef _DEBUG
#include <cassert>
#endif
#include <algorithm>
#include <limits>

namespace fms::pwflat {

	// f(u) = f[i], t[i-1] < u <= t[i], 0 <= i < n
	// f(u) = NaN, u < 0
	// f(u) = _f, u > t[n-1]
	template<class T, class F>
	inline F value(T u, unsigned n, const T* t, const F* f, 
		F _f = std::numeric_limits<F>::quiet_NaN())
	{
		if (u < 0) {
			return std::numeric_limits<F>::quiet_NaN();
		}

		auto i = std::lower_bound(t, t + n, u);

		return i == t + n ? _f : f[i - t];
	}

#ifdef _DEBUG
	inline int value_test()
	{
		double t[] = { 1,2,3 };
		double f[] = { .1,.2,.3 };
		assert(_isnan(value(-1., 3, t, f)));
		assert(f[0] == value(0., 3, t, f));
		double t0 = 0;
		for (int i : {0, 1, 2}) {
			assert(f[i] == value(t[i], 3, t, f));
			double ti = (t[i] + t0) / 2;
			assert(f[i] == value(ti, 3, t, f));
			t0 = t[i];
		}
		assert(4 == value(t[2] + 1, 3, t, f, 4.));

		return 0;
	}
#endif // _DEBUG

	// int_u^v f(t) dt
	template<class T, class F>
	inline F integral(T u, T v, unsigned n, const T* t, const F* f,
		F _f = std::numeric_limits<F>::quiet_NaN())
	{
		if (u < 0 || u > v) {
			return std::numeric_limits<F>::quiet_NaN();
		}

		auto i = std::lower_bound(t, t + n, u);
		if (i == t + n) {
			return _f;
		}

		return f[i - t];
	}

#ifdef _DEBUG
	inline int integral_test()
	{
		return 0;
	}
#endif // _DEBUG

} // namespace fms
