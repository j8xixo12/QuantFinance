#include <iostream>
#include <iomanip>
#include "utility.hpp"

int main(int argc, char* argv[]) {
    double a = -1.946981; 
    double b = 3.1876528; 
    double rho = 0.054292586;
    double AL = -8.0;
    int N = 20; // Number of subdivisions
    std::cout << "*A&S 26.3.3 Gauss Legendre: " << std::setprecision(16)
    << BivariateNormalGaussLegendre(a, b, rho, AL) << std::endl;
    std::cout << "*A&S 26.3.3: " << std::setprecision(16)
    << BivariateNormal(a, b, rho, AL, N) << std::endl;
    std::cout << "*A&S 26.3.3 Extrapolated: " << std::setprecision(16)
    << BivariateNormalExtrapolate(a, b, rho, AL, N) << std::endl;
    return 0;
}