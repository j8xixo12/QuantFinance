#include <iostream>
#include <random>
#include <memory>
#include <fstream>
#include "RNGenerator.hpp"
#include "Sde.hpp"
#include "Pricer.hpp"
#include "SUD.hpp"
#include "FdmEuler.hpp"

int main(int argc, char* argv[]) {
    int NSim = 50;
    int NT = 500;
    double driftCoefficient = 0.08;
    double diffusionCoefficient = 0.3;
    double dividendYield = 0.0; 
    double initialCondition = 60.0; 
    double expiry = 0.25;

    auto sde = std::shared_ptr<GBM<double>> (std::make_shared<GBM<double>>(driftCoefficient, 
                                                diffusionCoefficient,
                                                dividendYield,
                                                initialCondition,
                                                expiry));
    
    double K = 65.0;
    // Factories for objects in context diagram
    std::function<double(double)> payoffPut = [&K](double x) {return std::max<double>(0.0, K - x); }; 
    std::function<double(double)> payoffCall = [&K](double x) {return std::max<double>(0.0, x - K); }; 

    double r = 0.08; 
    double T = 0.25;

    std::function<double()> discounter = [&r, &T]() { return std::exp(-r * T); };
    auto pricerCall = std::shared_ptr<EuropeanPricer> (std::make_shared<EuropeanPricer>(payoffCall, discounter));
    auto pricerPut = std::shared_ptr<EuropeanPricer> (std::make_shared<EuropeanPricer>(payoffPut, discounter));
    auto fdm = std::shared_ptr<EulerFdm<GBM<double>>> (std::make_shared<EulerFdm<GBM<double>>>(sde, NT));
    auto rngPM = std::shared_ptr<PolarMarsaglia<double, std::uniform_real_distribution>>
                    (std::make_shared<PolarMarsaglia<double, std::uniform_real_distribution>>(0.0, 1.0));
    
    SUD<GBM<double>, EuropeanPricer, EulerFdm<GBM<double>>, PolarMarsaglia<double, std::uniform_real_distribution>>
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