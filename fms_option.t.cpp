// fms_variate_normal.t.cpp - Test fms::variate::normal
// Only test in debug mode
#ifdef _DEBUG
#include <cassert>
#include <algorithm>
#include "fms_option.h"
#include "fms_derivative.h"

using namespace fms;
using namespace fms::option;

int fms_option_vega_test()
{
	for (auto f : sequence(50, 100, 10)) {
		for (auto s : sequence(0.01, 1, 0.01)) {
			for (auto k : sequence(50, 100, 10)) {
				auto vs = [f, k](double s) { return option::value(f, s, k); };
				for (double h : {.001, .0001, .00001}) {
					{ // test vega
						double dv = option::vega(f, s, k);
						assert((derivative_test<double, double>(vs, s, h, dv, 1., 100000.)));
					}
				}
			}
		}
	}

	return 0;
}
//int fms_option_test_ = fms_option_test();

int fms_option_value_test(double f = 100, double s = 0.1, double k = 100, unsigned n = 0)
{
	auto v = [s,k,n](double f) { return value(f, s, k, n); };
	{
		for (double h : {0.01, 0.001, 0.0001}) {
			double dv = value(f, s, k, n + 1);
			double dddv = value(f, s, k, n + 3);
			assert((derivative_test<double, double>(v, f, h, dv, dddv, 100.)));
		}
	}

	return 0;
}
int fms_option_value_test_ = fms_option_value_test();

#endif // _DEBUG