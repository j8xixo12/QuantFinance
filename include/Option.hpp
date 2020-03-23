#ifndef Option_hpp_
#define Option_hpp_

class Option
{
public:
    explicit Option(const double r, const double sig, const double K,
                    const double T_, int type)
                    : r_(r), sig_(sig), K_(K), T_(T_), type_(type) {};

	double r_;		// Interest rate
	double sig_;		// Volatility
	double K_;		// Strike price
	double T_;		// Expiry date
	// double U_;		// Current underlying price (e.g. stock, forward)
	// double b_;		// Cost of carry

    int type_;		// 1 == Call, 2 == Put

	double payoff(double S)const
	{
		if (type_ == 1)
		{
			if (S > K_)
				return (S - K_);
		
			return 0.0;
		}
		else
		{
			if (S < K_)
				return -(S - K_);
			return 0.0;
		}
	}
};

#endif // Option_hpp_