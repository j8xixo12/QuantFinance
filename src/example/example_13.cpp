#include <iostream>
#include <functional>
#include <boost/numeric/ublas/matrix.hpp>
#include <fstream>
#include "TwoFactorPDE.hpp"

int main(int argc, char* argv[]) {
    double r = 0.049;
    double T = 0.5;
    // Values where price is calculated S/I == 1 
    double S = 100.; 
    double I = 100.;
    double K = 110.;
    double lambda = 0.5;
        // Define payoff as a lambda function
    double xMax = 1.0; // x far field 
    double yMax = 1.0; // y FF
    double xHotSpot = 0.5; 
    double yHotSpot = 0.5; 
    double scale1 = S * (1.0 - xHotSpot) / xHotSpot; 
    double scale2 = I * (1.0 - yHotSpot) / yHotSpot;
    // Values where price is calculated
    // S/I == 1
    // Define payoff as a lambda function 
    // Type cp = Type::Put;
    // if (cp == Type::Call){
    // auto payoff = [=](double S, double A)->double {return std::max(S - K, 0.0);};
    // } else {
        auto payoff = [=](double S, double I)->double {return std::max(K - S, 0.0);};
    // }

    AnchorPde anchorpde(scale1, scale2, r, lambda, K, T, xMax, yMax, payoff, Type::Put);
    // Domain
    TwoFactorPdeDomain<double> pdeDomain;
    CreateAnchorPdeDomain(pdeDomain, anchorpde, xMax, yMax, T, payoff);
    std::shared_ptr<TwoFactorAsianPde<double>> asianpde = CreateAsianPde(r, T, S, I, K, lambda, xMax, yMax,
                                                                            xHotSpot, yHotSpot, scale1, scale2, anchorpde);

    long NX = 100; 
    long NY = 100; 
    long NT = 100;

    std::vector<double> xmesh = pdeDomain.rx.mesh(NX); 
    std::vector<double> ymesh = pdeDomain.ry.mesh(NY); 
    std::vector<double> tmesh = pdeDomain.rt.mesh(NT);
    std::cout << "Now creating solver\n";
    TwoFactorAsianADESolver solver(pdeDomain, asianpde, xmesh, ymesh, tmesh);
    boost::numeric::ublas::matrix ans = solver.result();

    std::ofstream output("out.csv", std::ios::out);
    output << std::endl;
    for (std::size_t i = 0; i < ans.size1(); ++i) {
        for (std::size_t j = 0; j < ans.size2(); ++j) {
            output << xmesh[i] << " " << ymesh[j] << " " << ans(i, j) << std::endl;
        }
    }

    output.close();
    return 0;
}