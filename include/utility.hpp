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
double Pdf(const double x) {
    const double A = 1.0 / std::sqrt(2.0 * boost::math::constants::pi<double>());
    return A * std::exp(-x * x * 0.5);
}
// C++11 supports the error function
double Cdf(const double x) { // The approximation to the cumulative normal distribution
    return (0.5 * (1.0 - std::erf(-x / std::sqrt(2.0))));
}

#endif // UTILITY_HPP_