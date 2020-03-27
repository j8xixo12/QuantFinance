#include "LUTridiagonalSolver.hpp"
#include "DoubleSweep.hpp"
#include <iostream>
#include <vector>

/* Solve BVP u" = 1 in (0, 1) with u(0) = u(1) = 0 Solution u(x) = x(1-x)
Solve by FDM */

int main(int argc, char* argv[]) {
    using value_type = double;
    using Vector = std::vector<value_type>;

    std::size_t J = 20;
    std::cout << "Number of subdivisions J ";
    std::cin >> J;
    double h = 1.0 / static_cast<double>(J);
    // Boundary conditions 
    double BCL = 0.0; 
    double BCR = 0.0;

    // Double Sweep
    Vector a(J + 1, 1.0);
    Vector b(J + 1, -2.0);
    Vector c(J + 1, 1.0);
    Vector r(J + 1, -2.0 * h * h);// Right-hand side
    // Thomas algorithm 
    Vector a2(J - 1, 1.0);
    Vector b2(J - 1, -2.0);
    Vector c2(J - 1, 1.0);
    Vector r2(J - 1, -2.0 * h * h); // Right-hand side

    // Take the boundary conditions into consideration 
    r2[0] -= BCL;
    r2[r2.size()-1] -= BCR;
    LUTridiagonalSolver<double> mySolver2(a2, b2, c2, r2); 
    std::cout << "Matrix has a solution? " << mySolver2.diagonalDominant() << '\n'; 
    
    DoubleSweep<value_type> mySolver(a, b, c, r, BCL, BCR);

    Vector result = mySolver();
    Vector result2 = mySolver2();

    auto exact = [](double x) { return x*(1.0 - x); }; 
    double val = h; // Double Sweep
    // Compare output from Double Sweep and Thomas with each other 
    for (std::size_t j = 1; j < result.size()-1; ++j) { 
        // The values should be zero
            std::cout<< j << ", " << result[j]-result2[j-1] <<", "<< exact(val)<< '\n';
            val += h;
    }
    return 0;
}