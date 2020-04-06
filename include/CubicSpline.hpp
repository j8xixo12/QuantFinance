#ifndef CUBICSPLINE_HPP_
#define CUBICSPLINE_HPP_

#include <vector>
#include <functional>
#include <string>
#include <iomanip>
#include "LUTridiagonalSolver.hpp"
#include "utility.hpp"

enum CubicSplineBC {SecondDeriv, FirstDeriv};

class CubicSplineInterpolator {
    private:
        std::vector<double> x; // Abscissa x-values x_0,...,x_N
        std::vector<double> y; // Function values
        std::vector<double> h; // Array of mesh sizes in x direction
        CubicSplineBC type; // Type of BC
        std::vector<double> M; // Moments of spline
        std::vector<double> A, B, C, r_; // Input arrays for LU
        
        // For first order derivatives
        double a, b;

        void CalculateVectors(std::vector<double>& r) {
            // A, B, C and r
            std::size_t N = x.size() - 1;
            if (type == SecondDeriv) {
                C[0] = 0.0; 
                r[0] = 0.0;
                A[A.size() - 1] = 0.0;
                r[r.size() - 1] = 0.0; 
            } else {
                C[0] = 1.0;
                r[0] = 6.0 * ((y[1] - y[0]) / h[1] - a) / h[1];
                A[A.size() - 1] = 1.0;
                r[r.size() - 1] = 6.0 * (b - ((y[N] - y[N-1]) / h[N])) / h[N]; 
            }
            double tmp;
            for (std::size_t j = 1; j < x.size() - 1; ++j) {
                double fac = 1.0 / (h[j] + h[j + 1]); 
                C[j] = h[j + 1] * fac;
                A[j] = h[j] * fac;
                tmp = ((y[j + 1] - y[j]) / h[j + 1]) - ((y[j] - y[j - 1]) / h[j]);
                r[j] = (6.0 * tmp) * fac; 
            }
        }

        std::size_t findAbscissa(const std::vector<double>& x, double xvar) const { // Will give index of LHS value <= xvar.
            if (xvar < x[0] || xvar > x[x.size() - 1]) {
                std::string s = "\nValue " + std::to_string(xvar)
                + " not in range "
                + "(" + std::to_string(x[0]) + ","
                + std::to_string(x[x.size() - 1]) + ")";
                
                throw std::out_of_range(s);
            }
            // This algo has log complexity
            auto posA = std::lower_bound(std::begin(x), std::end(x), xvar);
            std::size_t index = std::distance(std::begin(x), posA);
            return index;
        }

    public:
        CubicSplineInterpolator(const std::vector<double>& xarr,
                                const std::vector<double>& yarr,
                                CubicSplineBC BCType,
                                double alpha = 0., double beta = 0.) {
            // Arrays must have the same size
            x = xarr;
            y = yarr;
            type = BCType;
            a = alpha; // LHS
            b = beta; // RHS 
            std::size_t N = xarr.size();
            // Calculate array of offset
            h = std::vector<double>(N, 0.0);
            for (std::size_t j = 1; j < h.size(); ++j) {
                h[j] = x[j] - x[j-1]; 
            }

            // All arrays have start index 1
            // Compared to the equations in the book, M(j) ––> M(j+1) 
            M = std::vector<double>(N, 0.0); // Solution
            // LU Coefficients
            A = std::vector<double>(N, 0.0);
            B = std::vector<double>(N, 2.0);// Diagonal vector, constant == 2 
            C = std::vector<double>(N, 0.0);
            r_ = std::vector<double>(N, 0.0);

            LUTridiagonalSolver<double> mySolver(A, B, C, r_, std::bind(&CubicSplineInterpolator::CalculateVectors, this, std::placeholders::_1));
            M = mySolver.solve();
        }

        CubicSplineInterpolator(const std::vector<double>& xarr,
                                const std::function<double (double)>& fun,
                                CubicSplineBC BCType,
                                double alpha = 0.0, double beta = 0.0) {
            // Continuous function value case
            std::size_t N = xarr.size();
            // Arrays must have the same size
            x = xarr;
            y = CreateDiscreteFunction(xarr, fun);
            type = BCType;
            a = alpha;   // LHS
            b = beta;    // RHS
            // Calculate array of offset
            h = std::vector<double>(N, 0.0);
            for (std::size_t j = 1; j < h.size(); ++j) {
                h[j] = x[j] - x[j - 1]; 
            }
            // All arrays have start index 1
            // Compared to the equations in the book, M(j) ––> M(j+1)
            M = std::vector<double>(N, 0.0); // Solution
            // LU Coefficients
            A = std::vector<double>(N, 0.0);
            B = std::vector<double>(N, 2.0);// Diagonal vector, constant == 2
            C = std::vector<double>(N, 0.0); 
            r_ = std::vector<double>(N, 0.0);

            LUTridiagonalSolver<double> mySolver(A, B, C, r_, std::bind(&CubicSplineInterpolator::CalculateVectors, this, std::placeholders::_1));
            mySolver.solve();
        }

        double Solve(double xvar) const { // Find the interpolated valued at a value x)
            std::size_t j = findAbscissa(x, xvar);// Index of LHS value <= x // Now use the formula
            double tmp = xvar - x[j];
            double tmpA = x[j + 1] - xvar;
            double tmp3 = tmp * tmp * tmp;
            double tmp4 = tmpA * tmpA * tmpA;
            double A = (y[j + 1] - y[j]) / h[j + 1] - (h[j + 1] * (M[j + 1] - M[j])) / 6.0; 
            double B = y[j] - (M[j] * h[j + 1] * h[j + 1]) / 6.0;
            double result = (M[j] * tmp4)/(6.0 * h[j + 1])
            + (M[j + 1] * tmp3) / (6.0 * h[j + 1])
            + (A * tmp) + B;
            return result;
        }

        double Integral() const { // ANW page 45
            double r1 = 0.0; 
            double r2 = 0.0;
            for (std::size_t j = 1; j < x.size(); ++j) {
                double hj = h[j];
                r1 += (y[j - 1] + y[j]) * hj;
                r2 -= (M[j - 1] + M[j]) * hj * hj * hj;
            }
            r1 *= 0.5;
            r2 /= 24.0;
            return r1 + r2;
        }

        std::tuple<double, double, double> ExtendedSolve(double xvar) const {
            // Solve for S and derivatives S', S"
            auto j = findAbscissa(x, xvar);
            double Mj = M[j]; 
            double Mjp1 = M[j + 1]; 
            double hjp1 = 1.0 / h[j + 1];

            double xj = x[j]; double xjp1 = x[j + 1]; 
            double tmp = xvar - x[j];
            double tmpA = x[j + 1] - xvar;
            double tmp3 = tmp * tmp * tmp;
            double tmp4 = tmpA * tmpA * tmpA;
            // S"
            double s2 = hjp1*(Mj*(xjp1 - xvar) + Mjp1*(xvar - xj));
            double A = (y[j + 1] - y[j]) / h[j + 1] - (h[j + 1] * (M[j + 1] - M[j])) / 6.0;
            double B = y[j] - (M[j] * h[j + 1] * h[j + 1]) / 6.0;

            // S'
            double s1 = -0.5 * hjp1 * Mj * tmpA * tmpA + 0.5 * hjp1 * Mjp1 * tmp * tmp + A;
            // S
            double s0 = (M[j] * tmp4) / (6.0 * h[j + 1])
            + (M[j + 1] * tmp3) / (6.0 * h[j + 1]) + (A * tmp)
            + B;
                return std::make_tuple(s0, s1, s2);
        }

        double Derivative(double xvar) const { // Solve for derivative S'
            auto j = findAbscissa(x, xvar);
            double Mj = M[j]; 
            double Mjp1 = M[j + 1];
            double hjp1 = 1.0 / h[j + 1];
            // double xj = x[j]; 
            // double xjp1 = x[j + 1]; 
            double tmp = xvar - x[j];
            double tmpA = x[j + 1] - xvar;
            // double tmp3 = tmp * tmp * tmp;
            // double tmp4 = tmpA * tmpA * tmpA;
            double A = (y[j + 1] - y[j]) / h[j + 1] - (h[j + 1] * (M[j + 1] - M[j])) / 6.0;
            // S'
            double s1 = -0.5 * hjp1 * Mj * tmpA * tmpA + 0.5 * hjp1 * Mjp1 * tmp * tmp + A;
            return s1; 
        }
};
#endif // CUBICSPLINE_HPP_