#ifndef PHYSICS
#define PHYSICS

#include "astr.hpp"

#include <string>
#include <utility>
#include <cmath>

struct {
    double g = 9.80665;
} Constants;

template <typename ArithType> class Vector_3D {
    private:
        /**
         * The x-component of the vector with respect to the origin (0, 0, 0)
         */
        ArithType XAxis = 0;
        /**
         * The y-component of the vector with respect to the origin (0, 0, 0)
         */
        ArithType YAxis = 0;
        /**
         * The z-component of the vector with respect to the origin (0, 0, 0)
         */
        ArithType ZAxis = 0;
        
        /**
         * The magnitude of the vector
         */
        ArithType Magnitude = 0;
        /**
         * The polar angle (Θ) within a spherical coordinate system stored in radians
         * This context is using the physics convention and would instead be phi using the mathematics conventions
         */
        double Theta = 0;
        /**
         * The azimuth angle (φ) within a spherical coordinate system stored in radians
         * This context is using the physics convention and would instead be theta using the mathematics conventions
         */
        double Phi = 0;

    public:
        /**
         * Default Constructor (set all values to zero)
         */
        Vector_3D() {}

        /**
         * Parametric Constructor - Input cartesian coordinates
         * 
         * @param xaxis The x-component of the vector
         * @param yaxis The y-component of the vector
         * @param zaxis The z-component of the vector
         */
        Vector_3D(ArithType xaxis, ArithType yaxis, ArithType zaxis) {
            static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");

            XAxis = xaxis;
            YAxis = yaxis;
            ZAxis = zaxis;

            double mag = sqrt(pow(XAxis, 2) + pow(YAxis, 2) + pow(ZAxis, 2));
            if (std::is_integral<ArithType>::value) {
                Magnitude = round(mag);
            } else {
                Magnitude = mag;
            }

            Theta = acos(ZAxis / mag);
            Phi = acos(XAxis / sqrt(pow(XAxis, 2) + pow(YAxis, 2)));
        }

        /**
         * Parametric Constructor - Input spherical coordinates
         */
        Vector_3D(ArithType magnitude, double theta, double phi, bool degrees) {
            static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");

            Magnitude = magnitude;
            Theta = theta;
            Phi = phi;

            if (degrees) {
                Theta *= M_PI / 180;
                Phi *= M_PI / 180;
            }

            XAxis = Magnitude * sin(Phi) * cos(Theta);
            YAxis = Magnitude * sin(Phi) * sin(Theta);
            ZAxis = Magnitude * cos(Phi);
        }

        ArithType getX() const {return XAxis;}
        ArithType getY() const {return YAxis;}
        ArithType getZ() const {return ZAxis;}
        ArithType getMag() const {return Magnitude;}
        double getTheta() const {return Theta;}
        double getPhi() const {return Phi;}

        std::string rToString(bool specifyPositive = false) const {return "(" + astr::toString(XAxis, specifyPositive) + ", " + astr::toString(YAxis, specifyPositive) + ", " + astr::toString(ZAxis, specifyPositive) + ")";}
        std::string sToString(bool specifyPositive = false, bool degrees = true) const {return "(" + astr::toString(Magnitude, specifyPositive) + ", " + astr::toString(Theta * (degrees ? 180 / M_PI : 1), specifyPositive) + ", " + astr::toString(Phi * (degrees ? 180 / M_PI : 1), specifyPositive) + ")";}
        std::string rToString_Places(unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) const {return "(" + astr::toString_Places(XAxis, beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(YAxis, specifyPositive) + ", " + astr::toString_Places(ZAxis, specifyPositive) + ")";}
        std::string sToString_Places(unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false, bool degrees = true) const {return "(" + astr::toString_Places(Magnitude, beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(Theta * (degrees ? 180 / M_PI : 1), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(Phi * (degrees ? 180 / M_PI : 1), beforeDecimal, afterDecimal, add, specifyPositive) + ")";}
        std::string rToString_Length(unsigned long length, bool leading = true, bool specifyPositive = false) const {return "(" + astr::toString_Length(XAxis, length, leading, specifyPositive) + ", " + astr::toString_Length(YAxis, length, leading, specifyPositive) + ", " + astr::toString_Length(ZAxis, length, leading, specifyPositive) + ")";}
        std::string sToString_Length(unsigned long length, bool leading = true, bool specifyPositive = false, bool degrees = true) const {return "(" + astr::toString_Length(Magnitude, length, leading, specifyPositive) + ", " + astr::toString_Length(Theta * (degrees ? 180 / M_PI : 1), length, leading, specifyPositive) + ", " + astr::toString_Length(Phi * (degrees ? 180 / M_PI : 1), length, leading, specifyPositive) + ")";}
        std::wstring rToWideString(bool specifyPositive = false) const {return astr::toWideString(rToString(specifyPositive));}
        std::wstring sToWideString(bool specifyPositive = false, bool degrees = true) const {return astr::toWideString(sToString(specifyPositive, degrees));}
        std::wstring rToWideString_Places(unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) const {return astr::toWideString(rToString_Places(beforeDecimal, afterDecimal, add, specifyPositive));}
        std::wstring sToWideString_Places(unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false, bool degrees = true) const {return astr::toWideString_Places(sToString_Places(beforeDecimal, afterDecimal, add, specifyPositive, degrees));}
        std::wstring rToWideString_Length(unsigned long length, bool leading = true, bool specifyPositive = false) const {return astr::toWideString(rToString_Length(length, leading, specifyPositive));}
        std::wstring sToWideString_Length(unsigned long length, bool leading = true, bool specifyPositive = false, bool degrees = true) const {return astr::toWideString(sToString_Length(length, leading, specifyPositive, degrees));}
};

std::pair<double, double> makeVector_Rad(double x, double y) {
    return {sqrt(pow(x, 2) + pow(y, 2)), atan(y / x)};
}
std::pair<double, double> makeVector_Deg(double x, double y) {
    std::pair<double, double> output = makeVector_Rad(x, y);
    return {output.first, output.second * 180 / M_PI};
}

std::pair<double, double> breakVector_Rad(double mag, double dir) {
    return {mag * cos(dir), mag * sin(dir)};
}
std::pair<double, double> breakVector_Deg(double mag, double dir) {
    return breakVector_Rad(mag, dir * M_PI / 180);
}

namespace kmtcs {
    class projectile {
        private:
            double Vi, Vix, Viy, Vf, Vfx, Vfy, Dy, Yi, Yf, My, Dx, Xi, Xf, Ai, Af, Tt, Pt;

        public:
            /**
             * Default Constructor (sets everything to zero)
             */
            projectile() {Vi = Vix = Viy = Vf = Vfx = Vfy = Dy = Yi = Yf = My = Dx = Xi = Xf = Ai = Af = Tt = Pt = 0;}
            /**
             * Paramteric Constructor - Heights only
             */
            projectile(double initialHeight, double finalHeight) {
                Yi = initialHeight;
                Yf = finalHeight;
            }
            /**
             * Parametric Constructor - Most generic projectile motion scenario
             */
            projectile(double launchVelocity, double launchAngle, double initialHeight, double finalHeight) {
                Yi = initialHeight;
                Yf = finalHeight;
            }
    };

    double launchAngle_Rad(double dx, double dy, double mag, bool useMin = true) {
        double angle1 = atan((pow(mag, 2) + sqrt(pow(mag, 4) - Constants.g * (Constants.g * pow(dx, 2) + 2 * dy * pow(mag, 2)))) / (Constants.g * dx));
        double angle2 = atan((pow(mag, 2) - sqrt(pow(mag, 4) - Constants.g * (Constants.g * pow(dx, 2) + 2 * dy * pow(mag, 2)))) / (Constants.g * dx));
        return useMin ? fmin(angle1, angle2) : fmax(angle1, angle2);
    }
    double launchAngle_Deg(double dx, double dy, double mag, bool useMin = true) {
        return launchAngle_Rad(dx, dy, mag, useMin) * 180 / M_PI;
    }
    
    double airTime_Deg(double mag, double angle, double dy) {
        double vy = mag * sin(angle * M_PI / 180);
        double time1 = (-vy + sqrt(pow(vy, 2) - 2 * Constants.g * dy)) / -Constants.g;
        double time2 = (-vy - sqrt(pow(vy, 2) - 2 * Constants.g * dy)) / -Constants.g;
        return fmax(time1, time2);
    }
    double peakTime_Deg(double mag, double angle) {
        return -(mag * sin(angle * M_PI / 180)) / -Constants.g;
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
}

#endif /* PHYSICS */
