#include <random>
#include <iostream>
#include <fstream>
#include "utility.hpp"
#include "SdeFactory.hpp"
#include "Option.hpp"
#include "FdmEuler.hpp"
#include "FdmHeun.hpp"

int main(int argc, char* argv[]) {
    std::ofstream output("out.csv", std::ios::out);

    Option myOption(100.0, 10.0, 0.1, 0.9, 0.03);
    double S0 = 20.0;
    // B. Get the SDE
    // Using factories
    std::shared_ptr<Sde<double>> sde = SdeFactory<double>::GetSde (myOption, S0, CEV);

    // C. Set up the FDM classes; 1 to N relationship between SDE and FDM 
    std::default_random_engine eng;
    std::normal_distribution<double> nor(0.0, 1.0); 
    std::random_device rd;
    eng.seed(rd());

    FdmEuler<double> fdm(sde, nor, eng);
    FdmHeun<double> fdm2(sde, nor, eng);
    // D. Joining up SDE and FDM in Mediator class 
    long NT = 200000;
    // auto vec = Path(*sde, fdm, NT);
    auto vec2 = Path(*sde, fdm2, NT);

    for (std::size_t i = 0; i < vec2.size(); ++i) {
        output << vec2[i] << std::endl;
    }
    output.close();

    return 0;
}