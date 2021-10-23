// xll_pwflat.cpp - Piecewise constant curves
#include "fms_pwflat.h"
#include "xll/xll/xll.h"

using namespace fms;
using namespace xll;
/*
AddIn xai_pwflat_value(
	Function(XLL_DOUBLE, "xll_pwflat_value", "PWF.VALUE")
	.Arguments({
		Arg(XLL_DOUBLE, "u", "is the value at which to calculate the function."),
		Arg(XLL_FP, "t", "is an array of times."),
		Arg(XLL_FP, "f", "is an array of values."),
		Arg(XLL_DOUBLE, "_f", "is an optional extrapolated value.")
		})
	.Category("")
	.FunctionHelp("Return the value of a piecewise constant function.")
);
double WINAPI xll_pwflat_value(double u, const _FPX* pt, const _FPX* pf, double _f)
{
#pragma XLLEXPORT
	return pwflat::value(u, size(*pt), pt->array, pf->array, _f);
}
*/