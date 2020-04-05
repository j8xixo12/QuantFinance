#include <iostream>
#include <fstream>
#include "IBvp.hpp"
#include "Option.hpp"
#include "BlackScholes.hpp"
#include "Range.hpp"
#include "utility.hpp"
#include "OptionCommand.hpp"
#include "CubicSpline.hpp"

int main(int argc, char* argv[]) {
    std::ofstream output("out.csv", std::ios::out);
    Option myOption;
    myOption.sig_ = 0.3; 
    myOption.K_ = 65.0; 
    myOption.T_ = 0.25;
    myOption.r_ = 0.08; 
    myOption.b_ = 0.08; 
    myOption.beta_ = 1.0;
    myOption.SMax_ = 325.0; 
    myOption.type = OptionType::Call; 

    BlackScholesPde myImp(myOption);
    Range<double> rangeX(0.0, myOption.SMax_); 
    Range<double> rangeT(0.0, myOption.T_);
    IBvp currentImp(myImp, rangeX, rangeT);

    long J = 500;
    long N = 500;
    std::cout << "Number of space divisions: ";
    std::cout << J << std::endl;
    std::cout << "Number of time divisions: ";
    std::cout << N << std::endl;

    CNIBVP fdmCN(currentImp, N, J);

    auto vCN = OptionPrice(fdmCN);
    auto sol = fdmCN.result(); 
    auto xarr = fdmCN.XValues();

    double h = rangeX.spread() / static_cast<double> (J);
    std::vector<double> zarr(xarr.size() - 2); 
    for (std::size_t j = 0; j < zarr.size(); ++j) {
        zarr[j] = xarr[j + 1];
    }

    // Delta array: eq. (22.14) but centred difference variant 
    std::vector<double> delta(zarr.size());

    for (std::size_t j = 0; j < zarr.size(); ++j) {
        delta[j] = (sol[j + 2] - sol[j]) / (2. * h); 
    }
    // Exact solution
    // OptionCommand(double strike, double expiration,
    // double riskFree, double costOfCarry, double volatility) 
    CallDelta cDelta (myOption.K_, myOption.T_, myOption.r_, myOption.b_, myOption.sig_);
    // PutDelta cDelta(myOption.K, myOption.T, myOption.r,
    // myOption.b, myOption.sig);
    std::vector<double> cDeltaPrices(zarr.size()); 
    for (std::size_t j = 0; j < zarr.size(); ++j) {
        cDeltaPrices[j] = cDelta.execute(zarr[j]); 
    }
    // Compute delta from cubic splines 
    CubicSplineInterpolator csi2(xarr, sol, SecondDeriv); 
    std::vector<double> splineDelta(zarr.size());
    for (std::size_t j = 0; j < zarr.size(); ++j)
    {
        splineDelta[j] = csi2.Derivative(zarr[j]); 
    }

    for (size_t j = 0; j < zarr.size(); ++j) {
        output << j << '\t' << cDeltaPrices[j] << '\t' << delta[j] << '\t' << splineDelta[j] << std::endl;
    }

    output.close();
    return 0;
}
