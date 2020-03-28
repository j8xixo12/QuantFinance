#include "CubicSpline.hpp"
#include "utility.hpp"
#include <random>
#include <iostream>


double Pdf(double x) { // Probability function for Gauss fn.
    double A = 1.0 / sqrt(2.0 * 3.14159265358979323846);
    return A * exp(-x*x*0.5); 
}
double dPdfdx(double x) { // 1st derivative of probability function for Gauss fn.
    return -x * Pdf(x); 
}
double d2Pdfdx2(double x) { // 2nd derivative of probability function for Gauss fn.
    return -Pdf(x) + x*x*Pdf(x); 
}
double Cdf(double x) { // The approximation to the cumulative normal distribution
// C++11
    return 0.5* (1.0 + std::erf(x / std::sqrt(2.0))); 
}

int main(int argc, char* argv[]) {
    // Infinite interval, truncated to [a,b]
    double a = -6.0;
    double b = 6.0;
    std::size_t n = 200000;
    // double h = (b - a) / static_cast<double>(n);
    std::vector<double> xarr = CreateMesh(n, a, b);

    auto fun = Pdf;
    auto yarr = CreateDiscreteFunction(xarr, fun);
    // Generate a random number in [a,b] 
    std::default_random_engine eng; 
    std::random_device rd; 
    eng.seed(rd());
    std::uniform_real_distribution<double> dist(a, b);
    double xvar = dist(eng);

    CubicSplineInterpolator csi(xarr, yarr, SecondDeriv);
    try {
        double result = csi.Solve(xvar);
        std::cout << "Interpolated value at " << xvar << " " << std::setprecision(16) << result << ", is " << fun (xvar) << std::endl;
        std::cout << "Integral: " << csi.Integral() << std::endl; 
    } catch (std::exception& e) { // Catch not in range values
        std::cout << e.what() << std::endl; 
    }

    // Numerical Differentiation
    auto derivs = csi.ExtendedSolve(xvar);
    std::cout << "Derivatives, approx: " << std::get<0>(derivs) << ", "
    << std::get<1>(derivs) << ", "<< std::get<2>(derivs) << std::endl;
    std::cout << "Derivatives, exact: " << Pdf(xvar) << ", "
    << dPdfdx(xvar) << ", " << d2Pdfdx2(xvar) << std::endl;
    // 1st derivative
    xvar = dist(eng);
    std::cout << "1st derivative: " << xvar << ", " << dPdfdx(xvar)
    << ", " << csi.Derivative(xvar) << std::endl;

    return 0;
}