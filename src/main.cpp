#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
#include "LUTridiagonalSolver.hpp"
#include "DoubleSweep.hpp"

double InitialCondition(double x) {
    if (x >= 0.0 && x <= 0.5) {
        return 2.0 * x; 
    }
    return 2.0 * (1.0 - x); 
}
template<typename T>
using Vector = std::vector<T>;

int main(int argc, char *argv[]) {
    // BTCS scheme for the heat equation
    long J = 20;
    long N = 100;
    long choice = 1;
    double theta = (choice == 2) ? 0.5 : 1.0;
    double T = 1.0;

    double A = 0.0; // LHS
    double B = 1.0; // RHS
    double h = (B - A) / static_cast<double>(J); 
    double k = T / static_cast<double>(N);
    // Boundary conditions 
    double lambda = k / (h * h);

    // Constructors with size, start index and value (Diagonals of matrix) 
    // J intervals, thus J-1 UNKNOWN internal points
    // Dirichlet boundary conditions
    Vector<double> a(J-1,-lambda*theta);
    Vector<double> b(J-1,(1.0 + 2.0*lambda*theta));
    Vector<double> c(J-1,-lambda*theta);
    Vector<double> r(J-1, 0); // Right-hand side NOT CONSTANT ANYMORE
                                // Boundary conditions into consideration
    // Create mesh in space
    Vector<double> xarr(J + 1); 
    xarr[0] = A;
    for (std::size_t j = 1; j < xarr.size(); ++j) {
        xarr[j] = xarr[j - 1] + h; 
    }

    Vector<double> vecOld(xarr.size()); // At time n 
    Vector<double> vecNew(xarr.size()); // At time n+1
    // Initial condition
    for (std::size_t j = 0; j < vecOld.size(); j++) {
        vecOld[j] = InitialCondition(xarr[j]);
    }
    // We start at 1st time point
    double current = k;

    auto r_fun = [&vecOld, lambda, theta] (Vector<double> &r) -> void {
        for (std::size_t j = 1; j < r.size() - 1; ++j) {
            r[j] = (lambda * (1.0 - theta) * vecOld[j + 1])
            + (1.0 - (2.0 * lambda * (1.0 - theta))) * vecOld[j]
            + (lambda * (1.0 - theta) * vecOld[j - 1]);
        }
    };
    
    LUTridiagonalSolver<double> mySolver(a, b, c, r, r_fun);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    while (current <= T) {
     // Update at new time level n+1
        // Compute inhomogeneous term
        vecNew = mySolver.solve();
        vecOld = vecNew;
        current += k;
    }
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::cout << "Elapsed time:  " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us." << std::endl;

}