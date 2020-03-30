#ifndef OPTION_COMMAND_HPP_
#define OPTION_COMMAND_HPP_

#include <iostream>
#include "Option.hpp"
#include "utility.hpp"

class OptionCommand : protected Option {
    public:
        OptionCommand() : Option() {};
        explicit OptionCommand(double strike, double expiration,
                            double riskFree, double costOfCarry, double
                            volatility)
                                   : Option(strike, expiration, riskFree, volatility, costOfCarry) {}
        
    // Want to forbid copy constructor and assignment operator
        OptionCommand(const OptionCommand& c) = delete;
        OptionCommand& operator = (const OptionCommand& c) = delete;
        // The abstract interface
        virtual void execute(double S) = 0;
        // Implement as function object; example of Template Method Pattern
        virtual void operator () (double S)
        {
            // Call derived class' execute()
            execute(S);
        }

        virtual ~OptionCommand() {};
};

class CallPrice final : public OptionCommand {
    public:
        explicit CallPrice(double strike, double expiration, double riskFree,
                            double costOfCarry, double volatility)
            : OptionCommand(strike, expiration, riskFree, costOfCarry,
            volatility) {}

        virtual void execute(double S) override {
            double tmp = sig_ * std::sqrt(T_);
            double d1 = (std::log(S / K_) + (b_ + (sig_ * sig_) * 0.5) * T_) / tmp; 
            double d2 = d1 - tmp;
            std::cout << "Call Price: " <<
            (S * std::exp((b_ - r_) * T_) * N(d1)) - (K_ * std::
            exp(-r_ * T_) * N(d2)) << std::endl;
        } 
};

class CallDelta final : public OptionCommand {
    public:
        explicit CallDelta(double strike, double expiration, double riskFree,
                            double costOfCarry, double volatility)
            : OptionCommand(strike, expiration, riskFree, costOfCarry,
            volatility) {}

        virtual void execute(double S)  override {
            double tmp = sig_ * std::sqrt(T_);
            double d1 = (std::log(S / K_) + (b_ + (sig_ * sig_) * 0.5) * T_) / tmp;
            std::cout << "Call delta: " << std::exp((b_ - r_) * T_) * N(d1) << std::endl; 
        }
};

class CallGamma final : public OptionCommand {
    public:
        explicit CallGamma(double strike, double expiration, double riskFree,
                        double costOfCarry, double volatility)
        : OptionCommand(strike, expiration, riskFree, costOfCarry,
        volatility) {}

        virtual void execute(double S) override {
            double tmp = sig_ * std::sqrt(T_);
            double d1 = (std::log(S / K_) + (b_ + (sig_ * sig_) * 0.5) * T_) / tmp;
            std::cout << "Call Gamma: " << N(d1) * std::exp((b_ - r_) * T_) / (S * tmp) << std::endl;
        }   
};

#endif // OPTION_COMMAND_HPP_