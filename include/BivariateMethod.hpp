#ifndef BIVARIATE_METHOD_HPP_
#define BIVARIATE_METHOD_HPP_

#include <functional>

using Function = std::function<double (const double)>;

class BivariateNormal {
    private:
        Function Pdf_;
        Function Cdf_;
    public:
        BivariateNormal(Function Pdf, Function Cdf)
        : Pdf_(Pdf), Cdf_(Cdf) {};

        double operator () (const double a, const double b, const double rho, const double ALower, const long NX) {
            double res = 0.0;
            double hx = (a - ALower) / NX;
            double fac = 1.0 / std::sqrt(1 - rho * rho);
            double s = ALower + (0.5 * hx);
            double x = s;
            double w;
            for (long n = 0; n < NX; n += 1) {
                w = fac * (b - rho * (x));
                res += hx * Pdf_(x) * Cdf_(w);
                x += hx;
            }
            return res;
        }
};

class BivariateNormalRichardson {
    private:
        BivariateNormal Normal;
    public:
        BivariateNormalRichardson(Function Pdf, Function Cdf)
        : Normal(Pdf, Cdf) {};

        double operator () (const double a, const double b, const double rho, const double ALower, const long NX) {
            // Repeated Richardson extrapolation
            double uh = Normal(a, b, rho, ALower, NX);
            double uh2 = Normal(a, b, rho, ALower, 2 * NX);
            auto vh = (4 * uh2 - uh) / 3.0;
            double uh4 = Normal(a, b, rho, ALower, 4 * NX);
            auto vh2 = (4 * uh4 - uh2) / 3.0;
            return (16 * vh2 - vh) / 15;
        }
};

class BivariateNormalGaussLegendre {
    private:
        const std::vector<double> x = { -0.0765265211334973,
                                        0.0765265211334973, -0.2277858511416451,
                                        0.2277858511416451, -0.3737060887154195,
                                        0.3737060887154195, -0.5108670019508271,
                                        0.5108670019508271, -0.6360536807265150,
                                        0.6360536807265150, -0.7463319064601508,
                                        0.7463319064601508, -0.8391169718222188,
                                        0.8391169718222188, -0.9122344282513259,
                                        0.9122344282513259, -0.9639719272779138,
                                        0.9639719272779138, -0.9931285991850949,
                                        0.9931285991850949 };
        const std::vector<double> w = { 0.1527533871307258,
                                        0.1527533871307258, 0.1491729864726037,
                                        0.1491729864726037, 0.1420961093183820,
                                        0.1420961093183820, 0.1316886384491766,
                                        0.1316886384491766, 0.1181945319615184,
                                        0.1181945319615184, 0.1019301198172404,
                                        0.1019301198172404, 0.0832767415767048,
                                        0.0832767415767048, 0.0626720483341091,
                                        0.0626720483341091, 0.0406014298003869,
                                        0.0406014298003869, 0.0176140071391521,
                                        0.0176140071391521 };
        Function Pdf_;
        Function Cdf_;
    public:
        BivariateNormalGaussLegendre(Function Pdf, Function Cdf)
        : Pdf_(Pdf), Cdf_(Cdf) {};

        double operator () (const double a, const double b, const double rho, const double ALower) {
            // Using Gauss Legendre 20-point method
            long N = 20;
            double res = 0.0;
            double alpha = (a - ALower) / 2.0;
            double beta = (a + ALower) / 2.0;
            const double fac = 1.0 / std::sqrt(1 - rho * rho);
            double xl, wl;
            // Using loop unrolling for (possible) improved performance
            for (long i = 0; i < N; i += 4) {
                xl = alpha * x[i] + beta;
                wl = fac * (b - rho * (xl));
                res += w[i] * Pdf(xl) * Cdf(wl);
                xl = alpha * x[i + 1] + beta;
                wl = fac * (b - rho * (xl));
                res += w[i + 1] * Pdf(xl) * Cdf(wl);
                xl = alpha * x[i + 2] + beta;
                wl = fac * (b - rho * (xl));
                res += w[i + 2] * Pdf(xl) * Cdf(wl);
                xl = alpha * x[i + 3] + beta;
                wl = fac * (b - rho * (xl));
                res += w[i + 3] * Pdf(xl) * Cdf(wl);
            }
            return alpha * res;
        }
};

#endif // BIVARIATE_METHOD_HPP_