#ifndef AUTOREGRESSIVE_HPP_
#define AUTOREGRESSIVE_HPP_

#include <vector>
#include <limits>
#include <armadillo>
#include "utility.hpp"

template <typename T>
class AR {
    public:
        explicit AR(std::vector<T> series, int Order) : time_series(series), NumOrder(Order), R(NumOrder + 1, arma::fill::zeros)
                                                        , std_epsilon(0.0), AIC_min(0.0), lag(0) {
            phi.reserve(Order);
            AIC_min = std::numeric_limits<double>::max();

        }

        void FindAcf() {
            phi = ACF(time_series, NumOrder);
        }

        void YuleWalker() {
            T sum = 0.0;
            T C0 = 0.0;
            for(auto j = NumOrder; j < time_series.size(); ++j) {
                C0 += time_series[j] * time_series[j];
            }

            for (auto i = 0; i <= NumOrder; ++i) {
                for(auto j = NumOrder; j < time_series.size(); ++j) {
                    sum += time_series[j] * time_series[j - i];
                }
                R(i) = sum / C0;
                sum = 0.0;
            }
            
            arma::vec lhs(NumOrder, arma::fill::zeros);
            arma::mat Rmat(NumOrder, NumOrder, arma::fill::eye);

            for(auto i = 0; i < NumOrder; ++i) {
                lhs(i) = R(i + 1);
            }

            // lhs.print();

            for (auto i = 0; i < NumOrder; ++i) {
                for (auto j = i + 1; j < NumOrder; ++j) {
                    Rmat(i, j) = R(j);
                    Rmat(j, i) = R(j);
                }
            }

            auto temp = inv(Rmat) * lhs;
            psi = arma::conv_to<std::vector<T>>::from(temp);
        }

        void MinimizeAIC() {
            T sum = 0.0;
            T RSS = 0.0;
            T AIC = 0.0;
            T temp = 0.0;
            std::size_t i, j, k;
            for (i = 1; i < NumOrder; ++i) {
                for (j = i; j < time_series.size(); ++j) {
                    for(k = 1; k <= i; ++k) {
                        sum += psi[k] * time_series[j - k];
                    }
                    temp = (time_series[j] - sum);
                    RSS += temp * temp;
                    sum = 0.0;
                }
                AIC = std::log(RSS / static_cast<double>(time_series.size() - i)) + (2. * i) / static_cast<double>(time_series.size() - i);
                // std::cout << "AIC: " << AIC << '\t' << "Sigma: " <<(RSS / static_cast<double>(time_series.size() - i)) << std::endl;
                if (AIC < AIC_min) {
                    AIC_min = AIC;
                    lag = i;
                }
                RSS = 0.0;
            }
        }

        std::vector<T> Get_phi() { return phi; }
        T Get_std_epsilon() { return std_epsilon; }
        T Get_AIC_min() { return AIC_min; }
        int Get_lag() { return lag; }

    private:
        std::vector<T> time_series;
        std::vector<T> phi;
        std::vector<T> psi;
        int NumOrder;
        arma::vec R;
        T std_epsilon;
        T AIC_min;
        int lag;
};

#endif // AUTOREGRESSIVE_HPP_