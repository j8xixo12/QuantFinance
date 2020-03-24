#ifndef Option_hpp_
#define Option_hpp_

class Option
{
public:
    explicit Option(const double r, const double sig, const double K,
                    const double T_)
                    : r_(r), sig_(sig), K_(K), T_(T_) {};

	double r_;		// Interest rate
	double sig_;		// Volatility
	double K_;		// Strike price
	double T_;		// Expiry date
	// double U_;		// Current underlying price (e.g. stock, forward)
	// double b_;		// Cost of carry
};

#endif // Option_hpp_