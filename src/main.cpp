#include <iostream>
#include <functional>
#include "TwoFactorPDE.hpp"

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

int main(int argc, char* argv[]) {
    double r = 0.049;
    double T = 0.5;
    double S = 100; 
    double I = 100;
    double lambda = 0.5;
    double K = 100.;
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
    auto payoff = [=](double S, double A)->double {return std::max(S - K, 0.0);};
    // } else {
    //     auto payoff = [=](double S, double I)->double {return std::max(K - S, 0.0);};
    // }

    AnchorPde anchorpde(scale1, scale2, r, lambda, K, T, xMax, yMax, payoff, Type::Call);
    // Domain
    TwoFactorPdeDomain<double> pdeDomain;
    CreateAnchorPdeDomain(pdeDomain, anchorpde, xMax, yMax, T, payoff);
    std::shared_ptr<TwoFactorAsianPde<double>> asianpde = CreateAsianPde();

    long NX = 100; 
    long NY = 100; 
    long NT = 100;

    std::vector<double> xmesh = pdeDomain.rx.mesh(NX); 
    std::vector<double> ymesh = pdeDomain.ry.mesh(NY); 
    std::vector<double> tmesh = pdeDomain.rt.mesh(NT);
    std::cout << "Now creating solver\n";
    TwoFactorAsianADESolver solver(pdeDomain, asianpde, xmesh, ymesh, tmesh);
    solver.result();

    auto values = FindMeshValues(solver.xmesh, solver.ymesh, xHotSpot, yHotSpot);
    std::cout << "MaxA, MaxB: " << std::get<0>(values) << ", " << std::get<1>(values) << std::endl;
    return 0;
}