#include <iostream>
#include <random>
#include <memory>
#include <fstream>
#include "RNGenerator.hpp"
#include "Sde.hpp"
#include "SdeFactory.hpp"
#include "Pricer.hpp"
#include "SUD.hpp"
#include "FdmEuler.hpp"
#include "Option.hpp"

int main(int argc, char* argv[]) {
    int NSim = 50;
    int NT = 1000;
    Option opt;
    opt.r_ = 0.08;
    opt.sig_ = 0.3;
    opt.D_ = 0.0; 
    double initialCondition = 60.0; 
    opt.T_ = 1.0;
    double K = 65.0;

    auto sde = SdeFactory<double>::GetSde(opt, initialCondition, CEV);
    
    // Factories for objects in context diagram
    std::function<double(double)> payoffPut = [&K](double x) {return std::max<double>(0.0, K - x); }; 
    std::function<double(double)> payoffCall = [&K](double x) {return std::max<double>(0.0, x - K); }; 

    double r = 0.08; 
    double T = 0.25;

    std::function<double()> discounter = [&r, &T]() { return std::exp(-r * T); };
    auto pricerCall = std::shared_ptr<EuropeanPricer> (std::make_shared<EuropeanPricer>(payoffCall, discounter));
    auto pricerPut = std::shared_ptr<EuropeanPricer> (std::make_shared<EuropeanPricer>(payoffPut, discounter));
    auto fdm = std::shared_ptr<EulerFdm<double>> (std::make_shared<EulerFdm<double>>(sde, NT));
    auto rngPM = std::shared_ptr<PolarMarsaglia<double, std::uniform_real_distribution>>
                    (std::make_shared<PolarMarsaglia<double, std::uniform_real_distribution>>(0.0, 1.0));
    
    SUD<double, Sde, EuropeanPricer, EulerFdm, PolarMarsaglia<double, std::uniform_real_distribution>>
        s(sde, pricerPut, fdm, rngPM, NSim, NT);

    s.start();

    auto res = s.result();

    std::ofstream output("out.dat", std::ios::out);

    for (std::size_t i = 0; i < res[0].size(); ++i) {
        output << fdm->x[i] << ", ";
        for (std::size_t j = 0; j < res.size(); ++j) {
            output << res[j][i] << ", ";
        }
        output << std::endl;
    }

    return 0;
}