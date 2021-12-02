// fms_variate_normal.t.cpp - Test fms::variate::normal
#ifdef _DEBUG
// Only test in debug mode
#include <cassert>
#include <algorithm>
#include <random>
#include "fms_derivative.h"
#include "fms_monte_carlo.h"
#include "fms_option.h"
#include "fms_variate_normal.h"

using namespace fms;
using namespace fms::option;

double monte_carlo_option_value(double f, double s, double k, size_t n = 10000)
{
    std::function<double(double)> payoff;
    if (k < 0) {
        payoff = [k](double x) { return std::max(-k - x, 0.); };
    }
    else {
        payoff = [k](double x) { return std::max(x - k, 0.); };
    }

    std::default_random_engine dre;
    std::normal_distribution<double> X;
    auto p = [f, s, &payoff, &X, &dre]() {
        double F = f * exp(s * X(dre) - s * s / 2);

        return payoff(F);
    };

    return monte_carlo::average(n, p);
}

namespace fms::monte_carlo {
    template<class X, class S = std::invoke_result<X>::type>
    inline S std(size_t n, X& x) //Welford's online algo, compute the var in one pass
    {
        S s = 0;
        S avglast = 0;
        S avg = x();
        for (unsigned int m = 2; m <= n; ++m) {
            S temp = x();
            avglast = avg;
            avg += (temp - s) / m;
            s += ((temp - avglast) * (temp - avg) - s) / m;
        }

        return std::sqrt(s);
    }
}

double monte_carlo_option_implied(double f, double s, size_t n = 10000)
{

    std::default_random_engine dre;
    std::normal_distribution<double> X;
    auto p = [f, s, &X, &dre]() {
        double F = f * exp(s * X(dre) - s * s / 2);

        return F;
    };

    return monte_carlo::std(n, p);
}

// common to all tests
variate::normal N;
double fs[] = { 80, 90, 100, 110, 120 };
double ks[] = { 80, 90, 100, 110, 120 };
double ss[] = { .01, .02, .1, .2 };
int is[] = { 10000 };
int option_value_test()
{
	// double sd = 2; // two standard deviations
	// for fs
	// for ks !!! test both k and -k
	// for ss
	// for is
	// !!!add for loops above and fix up below
	double sd = 2;

	for (double f : fs) {
		for (double k : ks) {
			for (int k_sign : {1, -1}) {
				k = k * k_sign;
				for (double s : ss) {
					for (double n : is) {
						double stdev = sqrt(option::black::variance(N, f, s, k));
						double v = option::black::value(N, f, s, k);
						double vn = monte_carlo_option_value(f, s, k, n);
						assert(fabs(v - vn) <= stdev * sd / sqrt(n));
					}
				}
			}
		}
	}

	return 0;
}
int option_implied_test()
{
    double sd = 2; // two standard deviations
    for (int i = 0; i < sizeof(fs) / sizeof(fs[0]); i++)
    {
        for (int j = 0; j < sizeof(ks) / sizeof(ks[0]); j++) //!!! test both k and -k
        {
            for (int m = 0; m < sizeof(ss) / sizeof(ss[0]); m++)
            {
                for (int z = 0; z < sizeof(is) / sizeof(is[0]); z++)
                {
                    double f = fs[i], s = ss[m], k = ks[j];
                    double stdev = sqrt(option::black::variance(N, f, s, k));
                    int n = is[z];
                    double v0 = option::black::value(N, f, s, k);
                    double v = option::black::implied(N, f, v0, k, s);
                    double vn = monte_carlo_option_implied(f, s, n);
                    //double sd = 2;
                    assert(fabs(v - vn) <= stdev * sd / sqrt(n));

                    stdev = sqrt(option::black::variance(N, f, s, -k));
                    v0 = option::black::value(N, f, s, -k);
                    v = option::black::implied(N, f, v0, -k, s);
                    vn = monte_carlo_option_implied(f, s, n);
                    assert(fabs(v - vn) <= stdev * sd / sqrt(n));
                }
            }
        }
    }

    return 0;
}
int option_value_test_ = option_value_test();

int option_delta_test_ = 0;
int option_gamma_test_ = 0;
int option_vega_test_ = 0;
int option_implied_test_ = option_implied_test();
int option_variance_test_ = 0;

#endif // _DEBUG