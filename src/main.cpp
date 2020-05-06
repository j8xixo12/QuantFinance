#include <iostream>
#include <random>
#include <fstream>
#include "utility.hpp"
#include "RNGenerator.hpp"
#include "Autoregressive.hpp"
#include "Sde.hpp"
#include "SdeFactory.hpp"
#include "Pricer.hpp"
#include "SUD.hpp"
#include "Option.hpp"
#include "FdmEuler.hpp"
#include "FdmHeun.hpp"

int main(int argc, char* argv[]) {

    std::ofstream output("out.dat", std::ios::out);
    std::ofstream output1("out2.dat", std::ios::out);
    RNGenerator<double, std::normal_distribution> rng(0.0, 1.0);

    std::vector<double> data(1000000, 0.0);

    for(auto i = 1; i < data.size(); ++i) {
        data[i] = 0.5 * data[i - 1] + rng();
    }

    // int NSim = 1;
    // int NT = 36500;
    // Option opt;
    // opt.r_ = 0.08;
    // opt.sig_ = 0.3;
    // opt.D_ = 0.0; 
    // double initialCondition = 60.0; 
    // opt.T_ = 1.0;
    // double K = 65.0;

    // auto sde = SdeFactory<double>::GetSde(opt, initialCondition, GBM);
    
    // // Factories for objects in context diagram
    // std::function<double(double)> payoffPut = [&K](double x) {return std::max<double>(0.0, K - x); }; 
    // std::function<double(double)> payoffCall = [&K](double x) {return std::max<double>(0.0, x - K); }; 

    // double r = opt.r_; 
    // double T = opt.T_;

    // std::function<double()> discounter = [&r, &T]() { return std::exp(-r * T); };
    // auto pricerCall = std::shared_ptr<EuropeanPricer> (std::make_shared<EuropeanPricer>(payoffCall, discounter));
    // auto pricerPut = std::shared_ptr<EuropeanPricer> (std::make_shared<EuropeanPricer>(payoffPut, discounter));
    // auto fdm = std::shared_ptr<FdmHeun<double>> (std::make_shared<FdmHeun<double>>(sde, NT));
    // auto rngPM = std::shared_ptr<PolarMarsaglia<double, std::uniform_real_distribution>>
    //                 (std::make_shared<PolarMarsaglia<double, std::uniform_real_distribution>>(0.0, 1.0));
    
    // SUD<double, Sde, EuropeanPricer, FdmHeun, PolarMarsaglia<double, std::uniform_real_distribution>>
    //     s(sde, pricerPut, fdm, rngPM, NSim, NT);

    // s.start();

    // auto data = s.result();

    AR ar(data, 5);
    ar.YuleWalker();
    ar.FindAcf();
    ar.MinimizeAIC();
    std::cout << "AIC min: " << ar.Get_AIC_min() << '\t' << "Best lag: " << ar.Get_lag() << std::endl;
    auto vec = ar.Get_phi();

    std::vector<std::vector<double>> outarry;

    for (auto i = 0; i < vec.size(); ++i) {
        output1 << i << ", " << vec[i] << std::endl;
    }

    outarry.push_back(data);

    for (auto i = 0; i < outarry[0].size(); ++i) {
        output << i << ", ";
        for(auto j = 0; j < outarry.size(); ++j) {
            output << outarry[j][i] << ", ";
        }
        output << std::endl;
    }

    output.close();
    output1.close();

    return 0;
}