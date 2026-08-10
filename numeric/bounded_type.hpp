#ifndef FLICK_BOUNDED_TYPE_HPP
#define FLICK_BOUNDED_TYPE_HPP

#include <cmath>
#include <iostream>
#include <ratio>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace flick {

using negative_one = std::ratio<-1>;
using zero         = std::ratio<0>;
using one          = std::ratio<1>;
using two          = std::ratio<2>;
using three        = std::ratio<3>;
using four         = std::ratio<4>;
using five         = std::ratio<5>;
using fifteen      = std::ratio<15>;
using hundred      = std::ratio<100>;

using one_pi  = std::ratio<
    3141592653589793238,
    1000000000000000000>;

using two_pi  = std::ratio_multiply<std::ratio<2>, one_pi>;
using four_pi = std::ratio_multiply<std::ratio<4>, one_pi>;

using pi_half = std::ratio<
    3141592653589793238,
    2000000000000000000>;


template<class T, class LowerBound, class UpperBound>
class bounded_type {
    static_assert(std::is_arithmetic_v<T>,
                  "bounded_type requires an arithmetic type");

    static constexpr T lower_bound_value()
    {
        return static_cast<T>(LowerBound::num) /
               static_cast<T>(LowerBound::den);
    }

    static constexpr T upper_bound_value()
    {
        return static_cast<T>(UpperBound::num) /
               static_cast<T>(UpperBound::den);
    }

    static void validate(T value)
    {
        constexpr T lower = lower_bound_value();
        constexpr T upper = upper_bound_value();

        bool invalid = value < lower || value > upper;

        // Comparisons with NaN are always false, so check explicitly.
        if constexpr (std::is_floating_point_v<T>) {
            invalid = invalid || std::isnan(value);
        }

        if (invalid) {
            std::ostringstream ss;
            ss << "bounded type value " << value
               << " is not within the required interval "
               << lower << " to " << upper;

            throw std::invalid_argument(ss.str());
        }
    }

    T value_{};

public:
    bounded_type() = default;

    explicit bounded_type(T value)
        : value_{value}
    {
        validate(value_);
    }

    T operator()() const noexcept
    {
        return value_;
    }

    T value() const noexcept
    {
        return value_;
    }

    bounded_type& operator=(T new_value)
    {
        // Validate first, so the object remains unchanged if this throws.
        validate(new_value);
        value_ = new_value;

        return *this;
    }

    friend std::ostream& operator<<(
        std::ostream& os,
        const bounded_type& bt)
    {
        return os << bt.value_;
    }
};

} // namespace flick

#endif // FLICK_BOUNDED_TYPE_HPP
