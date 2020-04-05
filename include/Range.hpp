#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <vector>

template<typename Type = double> class Range {
    private:
        Type high;
        Type low;
    public:
        Range() {
            high = 0.0;
            low = 0.0;
        } // Default constructor
        Range(const Type& lo, const Type& hi) : high(hi), low(lo) {} // Low and high value 
        Range(const Range<Type>& other) : high(other.high), low(other.low) {} // Copy constructor
        ~Range() {}

        // Modifier functions
        void Low(const Type& t1) { low = t1; } // Sets the low value of current range 
        void High(const Type& t1) { high = t1; }// Sets the high value of current range
        // Accessing functions
        Type Low() const { return low; }
        Type High() const { return high; }
        Type spread() const { return high - low; }
        // Lowest value in range
        // Highest value in the range
        // High - Low value
        // Boolean functions
        bool left(const Type& value) const { return (value < low) ? true : false; } // Value to the left? 
        bool right(const Type& value) const { return (value > high) ? true : false; }// Value to the right? 
        bool contains(const Type& value) const { return ((low <= value) && (high >= value)) ? true : false; }// Contains value?
        std::vector<Type> mesh(std::size_t N) const {
            std::vector<Type> x(N + 1); 
            x[0] = low; 
            x[x.size()-1] = high;
            Type h = (high - low) / static_cast<Type>(N); 
            for (std::size_t j = 1; j < x.size() - 1; ++j) {
                x[j] = x[j - 1] + h; 
            }
            return x;
        }
};

#endif // RANGE_HPP_