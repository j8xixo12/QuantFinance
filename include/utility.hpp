#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include "Sde.hpp"
#include "Fdm.hpp"
#include "IBvp.hpp"

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

double N(double x) { 
    // aka CdfN(x)
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double n(double x) {
    // aka pdf(x)
    return std::exp((-x * x) / 2.) / std::sqrt(2.0 * boost::math::constants::pi<double>());
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

std::vector<double> OptionPrice(IBvpSolver& fdm) { // Compute option price using a FD scheme
    return fdm.result();
}

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

std::tuple<std::size_t, std::size_t >
      FindMeshValues(const std::vector<double>& xarr,
                     const std::vector<double>& yarr,
                     double x, double y) { 
    // Compute indices by searching in an array xarr for a 'threshold'
    // value x.
    // Find position of first element in vector that satisfies // the predicate d >= x.
    // Logarithmic complexity for random-access iterators
    auto posA = std::upper_bound(xarr.begin(), xarr.end(), x); 
    auto posB = std::upper_bound(yarr.begin(), yarr.end(), y);
    auto maxA = std::distance(xarr.begin(), posA); 
    auto maxB = std::distance(yarr.begin(), posB);
    return std::make_tuple(maxA, maxB);
}

template<typename T>
T mean(std::vector<T> &x) {
    auto n = x.size();
    T sum = 0.0;
    for (const auto &i : x) {
        sum += i;
    }
    return sum / static_cast<double>(n);
}

template<typename T>
std::vector<T> ACF(std::vector<T> &x, const int &lag) {
    T Mean = mean(x);
    std::vector<T> ret(lag, 0.0);

    auto N = x.size();
    T n = 0.0;
    T variance = 0.0;
    T x_i = 0.0;

    for (auto i = 0; i < N; ++i) {
        variance += (x[i] - Mean) * (x[i] - Mean);
    }
    variance = variance / N;
    int temp = 0;
    for (auto i = 0; i <= ret.size(); ++i) {
        n = 0.0;
        for (auto j = 0; j < N - temp; ++j) {
            x_i = x[j] - Mean;
            n += (x_i * (x[(j + i) % N] - Mean));
        }
        ret[i] = n / variance / (N - temp);
        temp++;
    }
    return ret;
}

template<typename T>
std::vector<T> MA(std::vector<T> &x, const int &span) {
    std::vector<T> ret(x.size(), 0.0);
    auto sum = 0.0;
    for (std::size_t i = 0; i < span; ++i) {
        sum += x[i];
        ret[i] = sum / span;
    }
    sum = 0.0;
    for (std::size_t i = span - 1; i < x.size(); ++i) {
        for (std::size_t j = 0; j < span; ++j) {
            sum += x[i - j];
        }
        ret[i] = sum / span;
        sum = 0.0;
    }
    return ret;
}
#endif // UTILITY_HPP_