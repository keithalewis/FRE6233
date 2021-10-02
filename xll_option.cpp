// xll_option.cpp - Black-Scholes/Merton option value and greeks.
#include "fms_option.h"
#include "xll/xll/xll.h"

#ifndef CATEGORY
#define CATEGORY "FRE6233"
#endif

using namespace xll;
using namespace fms;

// Create XML documentation and index.html in `$(TargetPath)` folder.
// Use `xsltproc file.xml -o file.html` to create HTML documentation.
#ifdef _DEBUG
Auto<Open> xao_template_docs([]() {

	return Documentation(CATEGORY, "Documentation for " CATEGORY ".");

});
#endif // _DEBUG

AddIn xai_option_moneyness(
	Function(XLL_DOUBLE, "xll_option_moneyness", "OPTION.MONEYNESS")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		})
	.FunctionHelp("Return the option moneyness.")
	.Category(CATEGORY)
	.Documentation(R"(
Moneyness is \((\log(k/f) + \sigma^2 t/2)/\sigma\).
)")
);
double WINAPI xll_option_moneyness(double f, double s, double k, double t)
{
#pragma XLLEXPORT
	return option::moneyness(f, s, fabs(k), t);
}

AddIn xai_option_value(
	Function(XLL_DOUBLE, "xll_option_value", "OPTION.VALUE")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		})
	.FunctionHelp("Return the option call (k > 0) or put (k < 0) value.")
	.Category(CATEGORY)
	.Documentation(R"(
Option value is \(E[\max\{F - k, 0\}]\) for a call
and \(E[\max\{k - F, 0\}]\) for a put.
)")
);
double WINAPI xll_option_value(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return option::value(f, sigma, k, t);
}

AddIn xai_option_delta(
	Function(XLL_DOUBLE, "xll_option_delta", "OPTION.DELTA")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		})
	.FunctionHelp("Return the option call (k > 0) or put (k < 0) delta.")
	.Category(CATEGORY)
	.Documentation(R"(
Option delta is the derivative of option value with respect to forward.
)")
);
double WINAPI xll_option_delta(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return option::delta(f, sigma, k, t);
}

AddIn xai_option_gamma(
	Function(XLL_DOUBLE, "xll_option_gamma", "OPTION.GAMMA")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		})
		.FunctionHelp("Return the option call (k > 0) or put (k < 0) gamma.")
	.Category(CATEGORY)
	.Documentation(R"(
Option gamma is the second derivative of option value with respect to forward.
)")
);
double WINAPI xll_option_gamma(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return option::gamma(f, sigma, k, t);
}

AddIn xai_option_vega(
	Function(XLL_DOUBLE, "xll_option_vega", "OPTION.VEGA")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		})
		.FunctionHelp("Return the option call (k > 0) or put (k < 0) vega.")
	.Category(CATEGORY)
	.Documentation(R"(
Option vega is the derivative of option value with respect to vol.
)")
);
double WINAPI xll_option_vega(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	return option::vega(f, sigma, k, t);
}

AddIn xai_option_theta(
	Function(XLL_DOUBLE, "xll_option_theta", "OPTION.THETA")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration"),
		Arg(XLL_DOUBLE, "_dt", "is an optional time increment. Default is 1/250."),
		})
		.FunctionHelp("Return the option call (k > 0) or put (k < 0) theta.")
	.Category(CATEGORY)
	.Documentation(R"(
Option theta is the derivative of option value with respect to time.
)")
);
double WINAPI xll_option_theta(double f, double sigma, double k, double t, double dt)
{
#pragma XLLEXPORT
	if (dt == 0) {
		dt = 1. / 250;
	}

	return option::theta(f, sigma, k, t, dt);
}

AddIn xai_option_implied(
	Function(XLL_DOUBLE, "xll_option_implied", "OPTION.IMPLIED")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "v", "is the option value."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		Arg(XLL_DOUBLE, "_sigma", "is an optional initial guess. Default is 0.1."),
		Arg(XLL_WORD, "_n", "is an optional maximum number of iterations. Default is 100."),
		Arg(XLL_DOUBLE, "_tol", "is an optional absolute tolerance. Default is square root of machine epsilon."),
		})
	.FunctionHelp("Return the option call (k > 0) or put (k < 0) implied vol.")
	.Category(CATEGORY)
	.Documentation(R"(
Option implied volatility is the inverse of value.
)")
);
double WINAPI xll_option_implied(double f, double v, double k, double t, double sigma, unsigned n, double tol)
{
#pragma XLLEXPORT
	return option::implied(f, v, k, t, sigma, n, tol);
}
