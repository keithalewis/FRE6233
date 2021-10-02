// fms_derivative.h - Test derivatives
#pragma once
#include <functional>
#include <limits>
#include <valarray>
namespace fms {

	inline constexpr double epsilon = std::numeric_limits<double>::epsilon();

	// sequence from b to e with increment dx
	inline std::valarray<double> sequence(double b, double e, double dx)
	{
		size_t n = static_cast<unsigned>(1 + (e - b) / dx);
		std::valarray<double> a(n);

		for (size_t i = 0; i < n; ++i) {
			a[i] = b + i * dx;
		}

		return a;
	}

	// symmetric difference quotient
	template<class X, class Y>
	inline auto difference_quotient(const std::function<Y(X)>& f, X dx)
	{
		return [dx, &f](X x) { 
			return (f(x + dx) - f(x - dx)) / (2 * dx); 
		};
	}

	// (f(x + h) - f(x - h))/2h = f'(x) + f'''(x) h^2/3! + ...
	template<class X, class Y>
	inline bool derivative_test(const std::function<Y(X)>& f, X x, X h, X df, X dddf, X tol = 1)
	{
		auto Df = difference_quotient(f, h);

		double lhs = Df(x) - df;
		double rhs = dddf * h * h / 6;

		bool b = std::fabs(lhs) < tol * std::fabs(rhs);
		if (!b) {
			tol = lhs / rhs;
			if (std::fabs(lhs) < sqrt(epsilon) && std::fabs(rhs) < sqrt(epsilon)) {
				b = true;
			}
		}

		return b;
	}

}
