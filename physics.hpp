#ifndef PHYSICS
#define PHYSICS

#include <utility>
#include <cmath>

struct {
    double g = 9.80665;
} Constants;

template <typename ArithType> std::pair<ArithType, double> makeVector(const ArithType &x, const ArithType &y, const bool &useRadians = true) {
    static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
    return {std::is_integral<ArithType>::value ? std::round(std::sqrt(x * x + y * y)) : std::sqrt(x * x + y * y), std::atan2(y, x) * (useRadians ? 1 : 180 / M_PI)};
}
template <typename ArithType> std::pair<ArithType, ArithType> breakVector(const ArithType &mag, const double &angle, const bool &useRadians = true) {
    static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
    if (std::is_integral<ArithType>::value) {return {std::round(mag * std::cos(angle * (useRadians ? 1 : 180 / M_PI))), std::round(mag * std::sin(angle * (useRadians ? 1 : 180 / M_PI)))};}
    return {mag * std::cos(angle * (useRadians ? 1 : M_PI / 180)), mag * std::sin(angle * (useRadians ? 1 : M_PI / 180))};
}

namespace kmtcs {
    template <typename ArithType> ArithType airTime(const ArithType &mag, const double &angle, const ArithType &dy, const bool &useRadians = true) {
        static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        const double vy = mag * std::sin(angle * (useRadians ? 1 : M_PI / 180));
        const double root = std::sqrt(vy * vy - 2 * Constants.g * dy);
        return std::is_integral<ArithType>::value ? std::round(std::fmax((-vy + root) / -Constants.g, (-vy - root) / -Constants.g)) : std::fmax((-vy + root) / -Constants.g, (-vy - root) / -Constants.g);
    }

    template <typename ArithType> ArithType peakTime(const ArithType &mag, const double &angle, const bool &useRadians = true) {
        static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        return std::is_integral<ArithType>::value ? std::round((mag * std::sin(angle * (useRadians ? 1 : M_PI / 180))) / Constants.g) : ((mag * std::sin(angle * (useRadians ? 1 : M_PI / 180))) / Constants.g);
    }

    template <typename ArithType> ArithType maxHeight(const ArithType &mag, const double &angle, const ArithType &y0, const bool &useRadians = true) {
        static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        const double peak = peakTime<double>(mag, angle, useRadians);
        return std::is_integral<ArithType>::value ? std::round(y0 + mag * std::sin(angle * (useRadians ? 1 : M_PI / 180)) * peak - 0.5 * Constants.g * peak * peak) : (y0 + mag * std::sin(angle * (useRadians ? 1 : M_PI / 180)) * peak - 0.5 * Constants.g * peak * peak);
    }

    template <typename ArithType> double launchAngle(const ArithType &dx, const ArithType &dy, const ArithType &mag, const bool &minimizeHeight = true, const bool &useRadians = true) {
        static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        const double root = std::sqrt(mag * mag * mag * mag - Constants.g * (Constants.g * dx * dx + 2 * dy * mag * mag));
        const double a1 = std::atan2(mag * mag + root, Constants.g * dx);
        const double a2 = std::atan2(mag * mag - root, Constants.g * dx);
        if (minimizeHeight) {
            if (maxHeight<ArithType>(mag, a1, 0, useRadians) < maxHeight<ArithType>(mag, a2, 0, useRadians)) {return a1 * (useRadians ? 1 : 180 / M_PI);}
            return a2 * (useRadians ? 1 : 180 / M_PI);
        }
        if (maxHeight<ArithType>(mag, a1, 0, useRadians) < maxHeight<ArithType>(mag, a2, 0, useRadians)) {return a2 * (useRadians ? 1 : 180 / M_PI);}
        return a1 * (useRadians ? 1 : 180 / M_PI);
    }

    template <typename ArithType> std::pair<ArithType, double> landingVector(const ArithType &mag, const double &angle, const ArithType &dy, const bool &useRadians = true) {
        static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        return makeVector<ArithType>(std::is_integral<ArithType>::value ? std::round(mag * std::cos(angle * (useRadians ? 1 : M_PI / 180))) : (mag * std::cos(angle * (useRadians ? 1 : M_PI / 180))), mag * std::sin(angle * (useRadians ? 1 : M_PI / 180)) - Constants.g * airTime<double>(mag, angle, dy, useRadians), useRadians);
    }
}

#endif /* PHYSICS */
