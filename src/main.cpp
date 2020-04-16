#include <iostream>
#include <memory>
#include <random>
#include <fstream>
#include "SdeFactory.hpp"
#include "Option.hpp"


int main(int argc, char* argv[]) {

    Option opt;
    opt.K_ = 65.0;
    opt.T_ = 1.0;
    opt.sig_ = 0.4;
    opt.r_  = 0.02;

    std::ofstream output("out.csv", std::ios::out);

    auto PutPayoff = [&opt] (double S)-> double {return std::max<double> (opt.K_ - S, 0.0);};


    double S_0 = 60.0;
    std::shared_ptr<Sde<double>> sde_ptr = SdeFactory<double>::GetSde(opt, S_0, SdeType::GBM);
    
    double x;
    double VOld = S_0;
    double VNew;
    long NT = 1000;
    std::cout << "Number of time steps: " << NT << std::endl;
    long NSIM = 50;
    std::cout << "Number of simulations: " << NSIM << std::endl;

    std::vector<std::vector<double>> v_data(NT, std::vector<double>(NSIM, 0.0));

    double M = static_cast<double>(NSIM);
    double dt = opt.T_ / static_cast<double> (NT); 
    double sqrdt = std::sqrt(dt);

    // Normal random number
    double dW;
    double price = 0.0; // Option price 
    double payoffT;
    double avgPayoffT = 0.0;
    double squaredPayoff = 0.0;
    double sumPriceT = 0.0;
    // Normal (0,1) rng
    std::default_random_engine dre; 
    std::normal_distribution<double> nor(0.0, 1.0);
    // Create a random number
    dW = nor(dre);
    long coun = 0; // Number of times S hits origin
    
    for (long i = 0; i < M; ++i) {
        if ((i / 100'00) * 100'00 == i) {// Give status after each 10000th iteration
            std::cout << i << std::endl;
        }

        VOld = S_0;
        x = 0.0;

        for (long index = 0; index < NT; ++index) {
            // Create a random number
            dW = nor(dre);
            // The FDM (in this case explicit Euler) 
            VNew = VOld + (dt * sde_ptr->drift(x, VOld))
                    + (sqrdt * sde_ptr->diffusion(x, VOld) * dW);
            VOld = VNew;
            x += dt;
            v_data[index][i] = VNew;
        }

        payoffT = PutPayoff(VNew); 
        sumPriceT += payoffT;
        avgPayoffT += payoffT / M;
        avgPayoffT *= avgPayoffT;
        squaredPayoff += (payoffT * payoffT);
    }

    for (long index = 0; index < NT; ++index) {
        output << index * dt << ",";
        for (long i = 0; i < M; i++) {
            output << v_data[index][i] << ","; 
        }
        output << std::endl;
    }
    
    output.close();
    // Finally, discounting the average price
    price = std::exp(-opt.r_ * opt.T_) * sumPriceT / M;
    std::cout << "Price, after discounting: " << price << ", " << std::endl;
    double SD = std::sqrt((squaredPayoff / M) - sumPriceT * sumPriceT / (M * M));
    std::cout << "Standard Deviation: " << SD << ", " << std::endl;
    double SE = SD / std::sqrt(M);
    std::cout << "Standard Error: " << SE << ", " << std::endl;

    return 0;
}