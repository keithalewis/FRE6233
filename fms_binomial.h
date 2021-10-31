// fms_binomial.h - binomial model
// P(X_n = +/-1) = 1/2
// W_n = X_1 + ... + X_n, W_0 = 0
// E[W_n] = 0, Var(W_n) = n
// P(W_n = n - 2k) = C(n,k)/2^n
// C(n,k) = n!/(n-k)!k!
#pragma once
#include <functional>

namespace fms::binomial {

	// E_k[f(W_n)], f:int -> double
	inline double E(int k, int n, const std::function<double(int)>& f)
	{
		if (k == n) {
			return f(n - 2 * k);
		}

		return (E(k + 1, n, f) + E(k - 1, n, f)) / 2;
	}

	// F = f e^{s W_n/sqrt(n)}/cosh^n(s/sqrt(n))
	inline double value(double f, double s, double k)
	{

	}

} // namespace fms
