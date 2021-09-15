// xll_template.cpp - Sample xll project.
#include "xll_template.h"

using namespace xll;
using namespace FRE6233;

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
		Arg(XLL_DOUBLE, "s", "is the vol."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		})
	.FunctionHelp("Return the option moneyness.")
	.Category(CATEGORY)
	.Documentation(R"(
Moneyness is \((\log(k/f) + s^2/2)/s\).
)")
);
double WINAPI xll_option_moneyness(double f, double s, double k)
{
#pragma XLLEXPORT
	return option::moneyness(f, s, fabs(k));
}

AddIn xai_option_value(
	Function(XLL_DOUBLE, "xll_option_value", "OPTION.VALUE")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "s", "is the vol."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		})
	.FunctionHelp("Return the option call (k > 0) or put (k < 0) value.")
	.Category(CATEGORY)
	.Documentation(R"(
Option value is \(E[\max\{F - k, 0\}]\) for a call
and \(E[\max\{k - F, 0\}]\) for a put.
)")
);
double WINAPI xll_option_value(double f, double s, double k)
{
#pragma XLLEXPORT
	// c - p = f - k
	return k < 0 ? option::value(f, s, -k) : option::value(f, s, k) + f - k;
}

AddIn xai_option_delta(
	Function(XLL_DOUBLE, "xll_option_delta", "OPTION.DELTA")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "s", "is the vol."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		})
	.FunctionHelp("Return the option call (k > 0) or put (k < 0) delta.")
	.Category(CATEGORY)
	.Documentation(R"(
Option delta is the derivative of option value with respect to forward.
)")
);
double WINAPI xll_option_delta(double f, double s, double k)
{
#pragma XLLEXPORT
	// dc - dp = 1
	return k < 0 ? option::delta(f, s, -k) : 1 + option::delta(f, s, k);
}
