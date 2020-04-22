#ifndef SUD_HPP_
#define SUD_HPP_

#include "Sde.hpp"
#include "RNGenerator.hpp"
#include "Fdm.hpp"
#include "Pricer.hpp"

template <typename T, template<class TT> typename Sde, typename Pricer, 
            template<class TTT> typename Fdm, typename RNGenerator>
class SUD : private Sde<T>, private Pricer, private Fdm<Sde<T>>, private RNGenerator { // System under discussion
    private:
        // Four main components
        std::shared_ptr<Sde<T>> sde;
        std::shared_ptr<Pricer> pricer;
        std::shared_ptr<Fdm<Sde<T>>> fdm;
        std::shared_ptr<RNGenerator> rng;

        int NSim;
        std::vector<std::vector<double>> res;
    public:
        SUD() {}
        // Generated path per simulation
        SUD(const std::shared_ptr<Sde<T>>& s,
            const std::shared_ptr<Pricer>& p,
            const std::shared_ptr<Fdm<Sde<T>>>& f,
            const std::shared_ptr<RNGenerator>& r,
            int numberSimulations, int NT)
            : sde(s), pricer(p), fdm(f), rng(r) {
                NSim = numberSimulations;
                res = std::vector(NSim + 1, std::vector<double>(NT + 1, 0.0));
        }

        void start() {
            double VOld, VNew;
            for (int i = 1; i <= NSim; ++i) { 
                // Calculate a path at each iteration
                if ((i / 100000) * 100000 == i) { // Give status after a given numbers of iterations
                    std::cout << i << ",";
                }
                VOld = sde->InitialCondition(); 
                res[i][0] = VOld; 
                for (std::size_t n = 1; n < res[i].size(); n++) {
                    VNew = fdm->advance(VOld, fdm->x[n - 1], fdm->k, (*rng)());
                    res[i][n] = VNew; 
                    VOld = VNew;
                }
                // Send path data to the Pricer(s) 
                // This step can be optimised 
                // pricer->ProcessPath(res[i]); 
                // pricer->PostProcess();
            }
        }

        std::vector<std::vector<double>> result() { return res; }
};

#endif // SUD_HPP_