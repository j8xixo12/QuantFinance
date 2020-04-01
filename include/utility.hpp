#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <vector>
#include <cmath>
#include <random>
#include <boost/math/constants/constants.hpp>
#include "Sde.hpp"
#include "Fdm.hpp"

std::vector<double> CreateDiscreteFunction(const std::vector<double>& xarr,
                                const std::function<double (double)>& fun) {
    std::size_t nx = xarr.size();

    std::vector<double> Vec (nx, 0.0);

    for (std::size_t i = 0; i < nx; ++i) {
        Vec[i] = fun(xarr[i]);
    }
    return Vec;
}

std::vector<double> CreateMesh(std::size_t n, double a, double b) { // Create a mesh of size n+1 on closed interval [a,b]
    std::vector<double> x(n + 1); 
    x[0] = a; 
    x[x.size()-1] = b;
    double h = (b - a) / static_cast<double>(n); 
    for (std::size_t j = 1; j < x.size() - 1; ++j) {
        x[j] = x[j - 1] + h; 
    }
    return x;
}

std::vector<double> CreateRefinedMesh(const std::vector<double> &mesh) {
    std::size_t n = 2 * mesh.size();
    std::vector<double> x(n + 1);
    x[0] = mesh[0];
    x[n] = mesh[mesh.size() - 1];
    double h = (mesh[mesh.size() - 1] - mesh[0]) / static_cast<double>(n);
    for (std::size_t i = 1; i < x.size() -1; ++i) {
        x[i] = x[i - 1] + h;
    }

    return x;
}

// Creating a path of an SDE approximated by FDM
template <class Derived, typename T = double>
std::vector<T> Path(const Sde<T>& sde,
                  Fdm<Derived, T>& fdm,
                  long NT) {
    std::vector<T> result(NT + 1);
    result[0] = sde.ic;
    double dt = sde.B / static_cast<T>(NT); 
    double tn = dt;
    for (std::size_t n = 1; n < result.size(); ++n) {
        result[n] = fdm.advance(result[n - 1], tn, dt);
        tn += dt; 
    }
    return result;
}


double N(double x) { 
    // aka CdfN(x)
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

template <typename T>  T generateRN(T a, T b) { // Generate a uniform random number in interval [a,b]
    std::default_random_engine eng; 
    std::random_device rd; 
    eng.seed(rd());
    std::uniform_real_distribution<T> uniformVar(a, b);
    return uniformVar(eng);
}

// Normal variates etc.
double Pdf(double x) {
    const double A = 1.0 / std::sqrt(2.0 * boost::math::constants::pi<double>());
    return A * std::exp(-x * x * 0.5);
}
// C++11 supports the error function
double Cdf(double x) { // The approximation to the cumulative normal distribution
    return (0.5 * (1.0 - std::erf(-x / std::sqrt(2.0))));
}

double BivariateNormalGaussLegendre(double a, double b, double rho,
    double ALower) { // 26.3.3 Abramowitz and Stegun
    // 20-point Gauss Legendre rule
    // Nodes
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
    // Weights
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
    // Using Gauss Legendre 20-point method
    long N = 20;
    double res = 0.0;
    double alpha = (a - ALower) / 2.0;
    double beta = (a + ALower) / 2.0;
    const double fac = 1.0 / std::sqrt(1 - rho*rho);
    double xl, wl;
    // Using loop unrolling for (possible) improved performance
    for (long i = 0; i < N; i += 4) {
        xl = alpha*x[i] + beta;
        wl = fac*(b - rho * (xl));
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
    return alpha*res;
}

double BivariateNormal(double a, double b, double rho,
                        double ALower, long NX) { // 26.3.3 Abramowitz and Stegun
    double res = 0.0;
    double hx = (a - ALower) / NX;
    const double fac = 1.0 / std::sqrt(1 - rho*rho);
    double s = ALower + (0.5 * hx);
    // double hx2 = hx / 2;
    double x = s;
    double w;
    for (long n = 0; n < NX; n += 1) {
        w = fac * (b - rho * (x));
        res += hx * Pdf(x) * Cdf(w);
        x += hx;
    }
    return res;
}

double BivariateNormalExtrapolate(double a, double b, double rho,
                                    double ALower, long NX) { // 26.3.3 Abramowitz and Stegun
    // Repeated Richardson extrapolation
    double uh = BivariateNormal(a, b, rho, ALower, NX);
    double uh2 = BivariateNormal(a, b, rho, ALower, 2 * NX);
    auto vh = (4 * uh2 - uh) / 3.0;
    double uh4 = BivariateNormal(a, b, rho, ALower, 4 * NX);
    auto vh2 = (4 * uh4 - uh2) / 3.0;
    return (16 * vh2 - vh) / 15;
}

#endif // UTILITY_HPP_