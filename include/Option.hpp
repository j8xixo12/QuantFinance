#ifndef Option_hpp_
#define Option_hpp_
#include <cmath>

enum OptionType {
	Call,
	Put
};

class Option {
	public:
		Option() : r_(0.0), sig_(0.0), K_(0.0), T_(0.0), S_(0.0), b_(0.0), q_(0.0), D_(0.0),
					beta_(0.0), SMax_(0.0) {};
		explicit Option(const double r, const double sig, const double K,
						const double T)
						: r_(r), sig_(sig), K_(K), T_(T) {};
		
		explicit Option(const double K, const double T, const double r, const double sig,
						const double b)
						: r_(r), sig_(sig), K_(K), T_(T), b_(b) {};
		Option(const Option &other) : r_(other.r_), sig_(other.sig_), K_(other.K_), T_(other.T_), b_(other.b_),
								q_(other.q_), D_(other.D_), beta_(other.beta_), SMax_(other.SMax_), type(other.type) {}

		double r_;		// Interest rate
		double sig_;	// Volatility
		double K_;		// Strike price
		double T_;		// Expiry date
		double S_;		// Current underlying price (e.g. stock, forward)
		double b_;		// Cost of carry
		double q_;
		double D_;
		double beta_;    // elasticity factor
		double SMax_;     // Far field condition
		OptionType type;
};

class PerpetualOption : public Option {
	public:
		double PutPrice() {
			double sig2 = sig_ * sig_;
			double fac = b_ / sig2 - 0.5; 
			fac *= fac;
			double y1 = 0.5 - b_ / sig2 + std::sqrt(fac + 2.0*r_ / sig2);
			double fac2 = ((y1 - 1.0)*S_) / (y1 * K_);
			return K_ * std::pow(fac2, y1) / (y1 - 1.0);
		}

		double CallPrice() {
			double sig2 = sig_ * sig_;
			double fac = b_ / sig2 - 0.5; 
			fac *= fac;
			double y2 = 0.5 - b_ / sig2 - std::sqrt(fac + 2.0*r_ / sig2);
			double fac2 = ((y2 - 1.0) * S_) / (y2 * K_);
			return K_ * std::pow(fac2, y2) / (1.0 - y2);
		}
};
#endif // Option_hpp_