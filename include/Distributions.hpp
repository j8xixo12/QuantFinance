#ifndef DISTRIBUTIONS_HPP_
#define DISTRIBUTIONS_HPP_

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions.hpp>

template<typename T1, template<typename T2> class Dist>
class Distribution {
    private:
        Dist<T1> dist;
    public:
        template<typename... Args>
        Distribution(Args... args) : dist(Dist<T1>(args ...)) {};

        T1 operator () (T1 x) const { 
            // pdf, as function object
            return boost::math::pdf(dist, x);
        }

        T1 pdf (T1 x) const { 
            // pdf
            return (*this)(x);
        }

        T1 mean() const {
            return boost::math::mean(dist);
        }

        T1 mode() const {
            return boost::math::mode(dist);
        }

        T1 median() const {
            return boost::math::median(dist);
        }

        T1 standard_deviation() const {
            return boost::math::standard_deviation(dist);
        }

        T1 variance() const {
            return boost::math::variance(dist);
        }

        T1 cdf(T1 x) const { 
            // cdf
            return boost::math::cdf(dist, x);
        }

        T1 hazard(T1 x) const { 
            // Hazard function pdf/(1-cdf)
            return boost::math::hazard(dist, x);
        }

        T1 chf(T1 x) const { 
            // Cumulative hazard function integral of hazard
            // from -inf to x
            return boost::math::chf(dist, x);
        }

        // The interval [A,B] in which distribution is defined
        T1 A() const {
            return boost::math::range(dist).first; 
        }

        T1 B() const {
            return boost::math::range(dist).second; 
        }

        std::pair<T1, T1> range() const {
            return boost::math::range(dist);
        }

        std::pair<T1, T1> support() const {
            return boost::math::support(dist);
        }

        T1 quantile(T1 p) const {
            return boost::math::quantile(dist, p);
        }

        T1 kurtosis() const {
            return boost::math::kurtosis(dist);
        }

        T1 skewness() const {
            return boost::math::skewness(dist);
        }
};

#endif // DISTRIBUTIONS_HPP_