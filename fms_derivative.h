// fms_derivative.h - Test derivatives
#pragma once
#include <functional>

namespace fms {

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

		return std::fabs(Df(x) - df) < tol * std::fabs(dddf * h * h / 6);
	}

}
