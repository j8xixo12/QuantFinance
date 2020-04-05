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
        virtual double execute(double S) = 0;
        // Implement as function object; example of Template Method Pattern
        virtual double operator () (double S)
        {
            // Call derived class' execute()
            return execute(S);
        }

        virtual ~OptionCommand() {};
};

class CallPrice final : public OptionCommand {
    public:
        explicit CallPrice(double strike, double expiration, double riskFree,
                            double costOfCarry, double volatility)
            : OptionCommand(strike, expiration, riskFree, costOfCarry,
            volatility) {}

        virtual double execute(double S) override {
            double tmp = sig_ * std::sqrt(T_);
            double d1 = (std::log(S / K_) + (b_ + (sig_ * sig_) * 0.5) * T_) / tmp; 
            double d2 = d1 - tmp;
            
            return (S * std::exp((b_ - r_) * T_) * N(d1)) - (K_ * std::
            exp(-r_ * T_) * N(d2));
        } 
};

class CallDelta final : public OptionCommand {
    public:
        explicit CallDelta(double strike, double expiration, double riskFree,
                            double costOfCarry, double volatility)
            : OptionCommand(strike, expiration, riskFree, costOfCarry,
            volatility) {}

        virtual double execute(double S)  override {
            double tmp = sig_ * std::sqrt(T_);
            double d1 = (std::log(S / K_) + (b_ + (sig_ * sig_) * 0.5) * T_) / tmp;
            return std::exp((b_ - r_) * T_) * N(d1); 
        }
};

class CallGamma final : public OptionCommand {
    public:
        explicit CallGamma(double strike, double expiration, double riskFree,
                        double costOfCarry, double volatility)
        : OptionCommand(strike, expiration, riskFree, costOfCarry,
        volatility) {}

        virtual double execute(double S) override {
            double tmp = sig_ * std::sqrt(T_);
            double d1 = (std::log(S / K_) + (b_ + (sig_ * sig_) * 0.5) * T_) / tmp;
            return n(d1) * std::exp((b_ - r_) * T_) / (S * tmp);
        }   
};

class CallVega final : public OptionCommand {
    public:
        explicit CallVega(double strike, double expiration, double riskFree,
                        double costOfCarry, double volatility)
        : OptionCommand(strike, expiration, riskFree, costOfCarry,
        volatility) {}

        virtual double execute(double S) override {
            double tmp = sig_ * std::sqrt(T_);
            double d1 = (std::log(S / K_) + (b_ + (sig_ * sig_) * 0.5) * T_) / tmp;
            return S * std::sqrt(T_) * std::exp((b_ - r_) * T_) * n(d1);
        }   
};

class CallTheta final : public OptionCommand {
    public:
        explicit CallTheta(double strike, double expiration, double riskFree,
                        double costOfCarry, double volatility)
        : OptionCommand(strike, expiration, riskFree, costOfCarry,
        volatility) {}

        virtual double execute(double S) override {
            double d1 = (std::log(S / K_) + (b_ + (sig_ * sig_) * 0.5) * T_) / (sig_ * std::sqrt(T_));
            double d2 = d1 - sig_ * std::sqrt(T_);
            return -((S * sig_ * std::exp(std::pow((b_ - r_), T_) * n(d1))) / (2. * std::sqrt(T_)))
            - (b_ - r_) * S * std::exp((b_ - r_) * T_) * N(d1) - r_ * K_ * std::exp(-r_ * T_) * N(d2);
        }   
};

#endif // OPTION_COMMAND_HPP_