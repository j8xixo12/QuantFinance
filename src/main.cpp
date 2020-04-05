#include <iostream>
#include <fstream>
#include "IBvp.hpp"
#include "Option.hpp"
#include "BlackScholes.hpp"
#include "Range.hpp"
#include "utility.hpp"

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
    myOption.type = OptionType::Put; 

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
    ADE_BC_CRTP fdmADEBC(currentImp, N, J);
    auto vCN = OptionPrice(fdmCN);
    auto vAD = OptionPrice(fdmADEBC);

    for (std::size_t i = 0; i <= N; ++i) {
        output << i << '\t' << vCN[i] << '\t' << vAD[i] <<std::endl;
    }

    output.close();
    return 0;
}
