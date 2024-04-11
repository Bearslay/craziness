#ifndef PHYSICS
#define PHYSICS

#include "astr.hpp"

#include <string>
#include <utility>
#include <cmath>
#include <vector>

/**
 * Collection of constants that can be used in various equations (and changed too in some cases)
 */
struct {
    /**
     * Acceleration due to gravity
     */
    double g = 9.80665;
} Constants;

template <typename ArithType> std::pair<ArithType, double> makeVector(ArithType x, ArithType y, bool useRadians = true) {
    static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
    return {std::sqrt(x * x + y * y), std::atan(y / x) * (useRadians ? 1 : 180 / M_PI)};
}
template <typename ArithType> std::pair<ArithType, ArithType> breakVector(ArithType mag, double dir, bool useRadians = true) {
    static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
    return {mag * std::cos(dir), mag * std::sin(dir) * (useRadians ? 1 : 180 / M_PI)};
}

namespace kmtcs {
    double launchAngle_Rad(double dx, double dy, double mag, bool useMin = true) {
        double angle1 = atan((pow(mag, 2) + sqrt(pow(mag, 4) - Constants.g * (Constants.g * pow(dx, 2) + 2 * dy * pow(mag, 2)))) / (Constants.g * dx));
        double angle2 = atan((pow(mag, 2) - sqrt(pow(mag, 4) - Constants.g * (Constants.g * pow(dx, 2) + 2 * dy * pow(mag, 2)))) / (Constants.g * dx));
        return useMin ? fmin(angle1, angle2) : fmax(angle1, angle2);
    }
    double launchAngle_Deg(double dx, double dy, double mag, bool useMin = true) {
        return launchAngle_Rad(dx, dy, mag, useMin) * 180 / M_PI;
    }

    template <typename ArithType> ArithType launchAngle(ArithType dx, ArithType dy, ArithType mag, bool useMinAngle = true, bool useRadians = true) {
        static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        const double root = std::sqrt(mag * mag * mag * mag - Constants.g * (Constants.g * dx * dx + 2 * dy * mag * mag));
        const double a1 = std::atan((mag * mag + root) / (Constants.g * dx));
        const double a2 = std::atan((mag * mag - root) / (Constants.g * dx));
        const double output = (useMinAngle ? std::fmin(a1, a2) : std::fmax(a1, a2)) * (useRadians ? 1 : 180 / M_PI);
        return std::is_integral<ArithType>::value ? std::round(output) : output;
    }
    
    double airTime_Deg(double mag, double angle, double dy) {
        double vy = mag * sin(angle * M_PI / 180);
        double time1 = (-vy + sqrt(pow(vy, 2) - 2 * Constants.g * dy)) / -Constants.g;
        double time2 = (-vy - sqrt(pow(vy, 2) - 2 * Constants.g * dy)) / -Constants.g;
        return fmax(time1, time2);
    }

    template <typename ArithType> ArithType airTime(ArithType mag, double angle, ArithType dy, bool useRadians = true) {
        static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        const double vy = mag * std::sin(angle * (useRadians ? 1 : M_PI / 180));
        const double root = std::sqrt(vy * vy - 2 * Constants.g * dy);
        const double t1 = (-vy + root) / -Constants.g;
        const double t2 = (-vy - root) / -Constants.g;
        return std::is_integral<ArithType>::value ? std::round(std::fmax(t1, t2)) : std::fmax(t1, t2);
    }

    double peakTime_Deg(double mag, double angle) {
        return -(mag * sin(angle * M_PI / 180)) / -Constants.g;
    }

    template <typename ArithType> ArithType peakTime(ArithType mag, double angle, bool useRadians = true) {

    }

    double maxHeight_Deg(double mag, double angle, double y0) {
        double peakTime = peakTime_Deg(mag, angle);
        return y0 + (mag * sin(angle * M_PI / 180)) * peakTime - 0.5 * Constants.g * pow(peakTime, 2);
    }

    double landingAngle_Deg(double mag, double angle, double dy) {
        double vy = (mag * sin(angle * M_PI / 180)) - (Constants.g * airTime_Deg(mag, angle, dy));
        return makeVector_Deg(mag * cos(angle * M_PI / 180), vy).second;
    }

    double landingVelocity_Deg(double mag, double angle, double dy) {
        double vy = (mag * sin(angle * M_PI / 180)) - (Constants.g * airTime_Deg(mag, angle, dy));
        return makeVector_Rad(mag * cos(angle * M_PI / 180), vy).first;
    }

    template <typename ArithType> std::pair<ArithType, double> landingVector(ArithType mag, double angle, ArithType dy, bool useRadians = true) {
        static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        const double vy = (mag * std::sin(angle) - (Constants.g * airTime<ArithType>(mag, angle, dy, useRadians)));
        return makeVector<ArithType>(mag * std::cos(angle * (useRadians ? 1 : M_PI / 180)), vy);
    }
}

#endif /* PHYSICS */
