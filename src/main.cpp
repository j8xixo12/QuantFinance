#include <iostream>
#include <random>
#include <fstream>
#include <RInside.h>
#include <Rcpp.h>
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

Rcpp::NumericVector acfc(Rcpp::NumericVector x, bool plot = true, int lagmax = 1) {
    Rcpp::Environment stats("package:stats");
    Rcpp::Function ri = stats["acf"];
    Rcpp::Function na_pass = stats["na.pass"];
    Rcpp::List result = ri(x, lagmax, "correlation", plot, na_pass);
    // Rcpp::NumericVector res = result["acf"];
    return result["acf"];
}

int main(int argc, char* argv[]) {

    int maxlag = 20;
    int NSim = 1;
    int NT = 31536;
    Option opt;
    opt.r_ = 0.08;
    opt.sig_ = 0.3;
    opt.D_ = 0.0; 
    double initialCondition = 60.0; 
    opt.T_ = 1.0;
    double K = 65.0;

    auto sde = SdeFactory<double>::GetSde(opt, initialCondition, GBM);
    
    // Factories for objects in context diagram
    std::function<double(double)> payoffPut = [&K](double x) {return std::max<double>(0.0, K - x); }; 
    std::function<double(double)> payoffCall = [&K](double x) {return std::max<double>(0.0, x - K); }; 

    double r = opt.r_; 
    double T = opt.T_;

    std::function<double()> discounter = [&r, &T]() { return std::exp(-r * T); };
    auto pricerCall = std::shared_ptr<EuropeanPricer> (std::make_shared<EuropeanPricer>(payoffCall, discounter));
    auto pricerPut = std::shared_ptr<EuropeanPricer> (std::make_shared<EuropeanPricer>(payoffPut, discounter));
    auto fdm = std::shared_ptr<FdmHeun<double>> (std::make_shared<FdmHeun<double>>(sde, NT));
    auto rngPM = std::shared_ptr<PolarMarsaglia<double, std::uniform_real_distribution>>
                    (std::make_shared<PolarMarsaglia<double, std::uniform_real_distribution>>(0.0, 1.0));
    
    SUD<double, Sde, EuropeanPricer, FdmHeun, PolarMarsaglia<double, std::uniform_real_distribution>>
        s(sde, pricerPut, fdm, rngPM, NSim, NT);

    s.start();

    auto data = s.result();

    AR ar(data[1], maxlag);
    ar.YuleWalker();
    ar.FindAcf();
    ar.MinimizeAIC();
    std::cout << "AIC min: " << ar.Get_AIC_min() << '\t' << "Best lag: " << ar.Get_lag() << std::endl;
    auto vec = ar.Get_phi();

    RInside R(argc, argv);

    R["acf"] = acfc(Rcpp::wrap(data[1]), false, maxlag);
    R["data"] = Rcpp::wrap(data[1]);
    R["phi"] = vec;

    R.parseEvalQ("x11();plot(acf, type='h')");
    R.parseEvalQ("x11();plot(data, type='l')");
    R.parseEvalQ("x11();plot(phi, type='h');Sys.sleep(30)");

    return 0;
}