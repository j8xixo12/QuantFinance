#ifndef TWOFACTOR_PDE_HPP_
#define TWOFACTOR_PDE_HPP_

#include <functional>
#include <memory>
#include <iostream>
#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include "Range.hpp"
#include "utility.hpp"


template <typename T> struct TwoFactorPdeDomain {
    // 1. Domain 
    Range<T> rx; 
    Range<T> ry; 
    Range<T> rt;
    // 2. (Dirichlet) Boundary conditions, 
    std::function<T (T x,T t)> LowerBC; 
    std::function<T (T x,T t)> UpperBC; 
    std::function<T (T y,T t)> LeftBC; 
    std::function<T (T y,T t)> RightBC;
    // 3. Initial condition
    std::function<T (T x,T y)> IC;
};

template <typename T> class TwoFactorPde { 
    // Model a convection-diffusion-reaction PDE as a
    // composition of universal function wrappers.
    public:
        // U_t = a11U_xx + a22U_yy + b1U_x + b2U_y + cU_xy + dU + F;
        // f = f(x,y,t) in general
        std::function<T (T,T,T)> a11;
        std::function<T (T,T,T)> a22;
        std::function<T (T,T,T)> c;
        std::function<T (T,T,T)> b1;
        std::function<T (T,T,T)> b2;
        std::function<T (T,T,T)> d;
        std::function<T (T,T,T)> F;
        TwoFactorPde() = default;
};

template <typename T> class TwoFactorAsianPde : public std::enable_shared_from_this<TwoFactorAsianPde<T>> { 
    // Model a convection-diffusion-reaction PDE as a composition
    // of universal function wrappers. Asian-style PDE.
    public:
        // U_t = a11U_xx + b1U_x + b2U_y + dU + F;
        // f = f(x,y,t) in general
        std::function<T (T,T,T)> a11;
        std::function<T (T,T,T)> b1;
        std::function<T (T,T,T)> b2;
        std::function<T (T,T,T)> d;
        std::function<T (T,T,T)> F;
        TwoFactorAsianPde() = default;
        TwoFactorAsianPde(std::function<T (T,T,T)> a11in, std::function<T (T,T,T)> b1in,
                        std::function<T (T,T,T)> b2in, std::function<T (T,T,T)> din,
                        std::function<T (T,T,T)> Fin) : a11(a11in), b1(b1in),
                                                            b2(b2in), d(din), F(Fin) {}
        TwoFactorAsianPde(const TwoFactorAsianPde& pde) : a11(pde.a11), b1(pde.b1),
                                                            b2(pde.b2), d(pde.d), F(pde.F) {}
};

enum Type {
    Call,
    Put
};

class AnchorPde {
    using Function = std::function<double(double S1, double S2)>;
    private:
        // Scale factors (hotspots)
        double scale1;
        double scale2;
        // Parameters of PDE
        double r;
        double lambda;
        double K;
        // Domain information
        double T;
        double xMax;
        double yMax;
        Function payoff;
        Type type;

    public:
        AnchorPde(double sc1, double sc2, double rin, double lambdain, double Kin,
                    double Tin, double xMax, double yMax, Function Fin, Type typein)
                    : scale1(sc1), scale2(sc2), r(rin), lambda(lambdain),
                      K(Kin), T(Tin), xMax(xMax), yMax(yMax), payoff(Fin), type(typein) {};
                      
        // Similar to an Asian PDE.
        double SIGMOID2(double t) { 
            // t = S/A
            // Example of Table 2
            constexpr double a = 0.6; // (vol for low S/I)
            constexpr double b = 0.15; // (vol for high S/I)
            constexpr double c = 10.0;
            constexpr double d = -0.3; // (when S=I vol is in the middle) const double ecd = std::exp(c*d);
            // Sigmoid2 function
            double val = std::exp(-c * (std::log(t) - d));
            double val2 = a + (b - a) / (1.0 + val);
            return val2;
        }

        double a11(double x, double y, double t) {
            double tx = 1.0 - x; 
            double ty = 1.0 - y;
            double S = scale1 * x / tx;
            double I = scale2 * y / ty;
            // Straight computation, i.e. function call 
            double s1 = SIGMOID2(S / I);
            return 0.5 * s1 * s1 * x * x * tx * tx;
        }

        double b1(double x, double y, double t) {
            double tx = 1.0 - x; 
            double ty = 1.0 - y;
            double S = scale1 * x / tx;
            double I = scale2 * y / ty;
            double s1 = SIGMOID2(S / I); 
            return r * x  *tx - s1 * s1 * x * x * tx;
        }

        double b2(double x, double y, double t) {
            double tx = 1.0 - x; 
            double ty = 1.0 - y;
            double S = scale1 * x / tx;
            double I = scale2 * y / ty;
            return lambda * (S - I) * ty * ty / scale2; 
        }

        double d(double x, double y, double t) {
            return -r; 
        }

        double F(double x, double y, double t) {
            return 0.; 
        }

        double IC(double x, double y) { 
            // Initial condition
            // Workaround/fudge to avoid spike 
            if (x == 1.0) x -= 0.0001;
            if (y == 1.0) y -= 0.0001;
            return payoff(scale1 * x / (1.0 - x), scale2 * y / (1.0 - y)); 
        }

        double BCLower(double x, double t) {
            return 0.;
        }

        double BCUpper(double x, double t) {
            // Tricky one for variable I
            return 0; // Force the solution to be BS classic. 
        }
        double BCLeft(double y, double t) { 
            // y is S2
            if (type == Call) {
                return 0.; // C 
            } else {
                return K * std::exp(-r * t); // P
            }
        }
        
        double BCRight(double y, double t) {
            if (type == Put) {
                return 0.; // P
            } else {
                double S = scale1 * xMax / (1 - xMax + 0.001); 
                return S - K * std::exp(-r * t); // C
            }
        }
};

namespace u = boost::numeric::ublas;

class TwoFactorAsianADESolver {
    private:
        TwoFactorPdeDomain<double> pdeDomain;
        std::shared_ptr<TwoFactorAsianPde<double>> pde;

        u::matrix<double> U;
        u::matrix<double> V;
        // Mesh-related data
        double hx, hy, delta_k, hx1, hy1, hx2, hy2;
        // Mesh-point values of coefficients
        double A, D, E, F, G; 
        // double B, C;  
        double t2, tx1, ty1;
        // Other variables
        double tprev, tnow, T;
        // std::size_t NX, NY, NT;

        public:
            u::matrix<double> MatNew;
            std::vector<double> xmesh;
            std::vector<double> ymesh;
            std::vector<double> tmesh;

            TwoFactorAsianADESolver(const TwoFactorPdeDomain<double>& domain,
                            const std::shared_ptr<TwoFactorAsianPde<double>> &pdein,
                            const std::vector<double>& xarr,
                            const std::vector<double>& yarr,
                            const std::vector<double>& tarr)
            : xmesh(xarr), ymesh(yarr), tmesh(tarr) {
                std::cout << "Two-factor ADE Asian classic version\n"; 
                pdeDomain = domain;
                pde = pdein;
                hx = pdeDomain.rx.spread() / static_cast<double>(xmesh.size() - 1); 
                hy = pdeDomain.ry.spread() / static_cast<double>(ymesh.size() - 1); 
                delta_k = pdeDomain.rt.spread() / static_cast<double>(tmesh.size() - 1);
                T = pdeDomain.rt.High();
                // Extra, handy variables 
                hx1 = 1.0 / hx;
                hx2 = 1.0 / (hx * hx);
                hy1 = 1.0 / hy;
                hy2 = 1.0 / (hy * hy);
                // Some optimising variables
                t2 = delta_k * hx2;
                tx1 = delta_k * hx1;
                ty1 = delta_k * hy1;
                // Initialise U, UOld, V, VOld, MatNew data structures
                // NumericMatrix(I rows, I columns, I rowStart, I columnStart);
                // Constructor with size & start index
                U = u::matrix<double>(xmesh.size(), ymesh.size());
                V = u::matrix<double>(xmesh.size(), ymesh.size());
                MatNew = u::matrix<double>(xmesh.size(), ymesh.size());
                initIC();
            }

            ~TwoFactorAsianADESolver() {}

            void initIC() { 
                // Utility function to initialise payoff function and BCs at t = 0
                tprev = tnow = tmesh[0];
                // Now initialise values in interior of interval using // the initial function 'IC' from the PDE
                for (std::size_t i = 0; i < xmesh.size(); ++i) {
                    for (std::size_t j = 0; j < ymesh.size(); ++j) {
                        MatNew(i,j) =  U(i, j) = V(i, j) = pdeDomain.IC(xmesh[i], ymesh[j]); 
                    }
                }
            }

            void calculateBC() {
                // Calculate the discrete BC on the FOUR edges of the boundary
                // Lower and Upper BC
                std::size_t lower = 0;
                std::size_t upper = ymesh.size()-1;
                for (std::size_t i = 0; i < xmesh.size(); ++i) {
                    MatNew(i, lower) = U(i, lower) = V(i, lower) = pdeDomain.LowerBC(xmesh[i], tnow);// N/A
                    MatNew(i, upper) = U(i, upper) = V(i, upper) = pdeDomain.UpperBC(xmesh[i], tnow);
                }
                // Left and Right BC
                std::size_t left = 0;
                std::size_t right = xmesh.size() -1;
                for (std::size_t j = 0; j < ymesh.size(); ++j) {
                    MatNew(left, j) = U(left, j) = V(left, j) = pdeDomain.LeftBC(ymesh[j], tnow);
                    MatNew(right, j) = U(right, j) = V(right, j) = pdeDomain.RightBC(ymesh[j], tnow);
                }
            }

            void calculate() {
                // Tells how to calculate sol. at n+1, Explicit ADE schemes
                double mx, my;
                int sgn;
                for (std::size_t i = 1; i <= xmesh.size()-2; ++i) {
                    mx = xmesh[i];
                    for (std::size_t j = 1; j <= ymesh.size()-2; ++j) {
                        // Create coefficients
                        // hU_t = aU_xx + bU_yy + cU_xy + dU_x + eU_y + fU + G;
                        // f = f(x,y,t)
                        my = ymesh[j];
                        A = pde->a11(mx, my, tnow) * t2;
                        D = 0.5 * pde->b1(mx, my, tnow) * tx1; 
                        // double sgnD = MySign(D);
                        E = pde->b2(mx, my, tnow) * ty1;
                        sgn = sign(E);
                        F = pde->d(mx, my, tnow) * delta_k; 
                        G = pde->F(mx, my, tnow) * delta_k;
                        // Larkin + 1st order upwind in y
                        U(i, j) = (U(i, j) * (1.0 - A) + U(i + 1, j) * (A + D)
                        + sgn * E * U(i, j + sgn) + U(i - 1, j) * (A - D) + G) / (1.0 + A - F + sgn * E);
                    } 
                }

                for (std::size_t i = xmesh.size()-2; i >= 1 ; --i) {
                    mx = xmesh[i];
                    for (std::size_t j = ymesh.size() - 2; j >= 1; --j) {
                        // Create coefficients
                        my = ymesh[j];
                        A = pde->a11(mx, my, tnow) * t2;
                        D = 0.5 * pde->b1(mx, my, tnow) * tx1; 
                        E = pde->b2(mx, my, tnow) * ty1;
                        sgn = sign(E);
                        F = pde->d(mx, my, tnow) * delta_k;
                        G = pde->F(mx, my, tnow) * delta_k;
                        // Larkin + 1st order upwind
                        V(i, j) = (V(i, j) * (1.0 - A) + V(i - 1, j) * (A - D)
                        + V(i + 1, j) * (A + D) + sgn * E * V(i, j + sgn) + G) / (1.0 + A - F + sgn*E);
                    }
                }

                for (std::size_t i = 0; i < MatNew.size1() ; ++i) {
                    for (std::size_t j = 0; j < MatNew.size2(); ++j) {
                        MatNew(i,j) = 0.5 * (U(i,j) + V(i,j)); 
                        U(i, j) = MatNew(i, j);
                        V(i, j) = MatNew(i, j);
                    } 
                }

            }

            u::matrix<double> result() { 
                // The result of the calculation
                std::cout << "result";
                for (std::size_t n = 1; n < tmesh.size(); ++n) {
                    if ((n / 100) * 100 == n) {
                        std::cout << n << ", ";
                    }
                    tnow = tmesh[n];
                    calculateBC(); // Calculate the BC at n+1
                    calculate(); // Calculate the solution at n+1
                    tprev = tnow;
                }

                return MatNew;
            }
};

std::shared_ptr<TwoFactorAsianPde<double>> CreateAsianPde(double r, double T, double S, double I,
                                                            double K, double lambda, double xMax,
                                                            double yMax, double xHotSpot, double yHotSpot,
                                                            double scale1, double scale2, AnchorPde &anchorpde) {

    auto a11 = std::bind(&AnchorPde::a11, &anchorpde, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    auto b1 = std::bind(&AnchorPde::b1, &anchorpde, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    auto b2 = std::bind(&AnchorPde::b2, &anchorpde, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    auto d = std::bind(&AnchorPde::d, &anchorpde, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    auto F = std::bind(&AnchorPde::F, &anchorpde, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

    return std::shared_ptr<TwoFactorAsianPde<double>>
      (std::make_shared<TwoFactorAsianPde<double>>(a11, b1, b2, d, F));
}

void CreateAnchorPdeDomain(TwoFactorPdeDomain<double>& pdeDomain,
                AnchorPde &pde,
                double xMax, double yMax, double T,
                const std::function<double(double, double)>& IC) {
    pdeDomain.rx = Range<double>(0.0, xMax);
    pdeDomain.ry = Range<double>(0.0, yMax);
    pdeDomain.rt = Range<double>(0.0, T);
    pdeDomain.LeftBC = std::bind(&AnchorPde::BCLeft, &pde, std::placeholders::_1, std::placeholders::_2); 
    pdeDomain.RightBC = std::bind(&AnchorPde::BCRight, &pde, std::placeholders::_1, std::placeholders::_2); 
    pdeDomain.UpperBC = std::bind(&AnchorPde::BCUpper, &pde, std::placeholders::_1, std::placeholders::_2); 
    pdeDomain.LowerBC = std::bind(&AnchorPde::BCLower, &pde, std::placeholders::_1, std::placeholders::_2);
    pdeDomain.IC = std::bind(&AnchorPde::IC, &pde, std::placeholders::_1, std::placeholders::_2);
}
#endif // TWOFACTOR_PDE_HPP_