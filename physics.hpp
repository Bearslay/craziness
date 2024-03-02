#ifndef PHYSICS
#define PHYSICS

#include "astr.hpp"

#include <string>
#include <utility>
#include <cmath>

/**
 * Collection of constants that can be used in various equations (and changed too in some cases)
 */
struct {
    /**
     * Acceleration due to gravity
     */
    double g = 9.80665;
} Constants;

/**
 * A vector with three dimensions that contains values for both a cartesian and spherical coordinate system
 * This is mainly intended for physics/math-based applications
 * 
 * @tparam ArithType An arithmetic data type (generally and integer or a float)
 */
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
         * 
         * When limited to integral data types, angles get real inaccurate for vectors of increasing magnitudes, so is always a double to account for that
         */
        double Theta = 0;
        /**
         * The azimuth angle (φ) within a spherical coordinate system stored in radians
         * This context is using the physics convention and would instead be theta using the mathematics conventions
         * 
         * When limited to integral data types, angles get real inaccurate for vectors of increasing magnitudes, so is always a double to account for that
         */
        double Phi = 0;

        /**
         * Calculate the cartesian components of the vector based off of its magnitude and directions
         */
        void calcCartesian() {
            if (std::is_integral<ArithType>::value) {
                XAxis = std::round(Magnitude * sin(Phi) * cos(Theta));
                YAxis = std::round(Magnitude * sin(Phi) * sin(Theta));
                ZAxis = std::round(Magnitude * cos(Phi));
                return;
            }

            XAxis = Magnitude * sin(Phi) * cos(Theta);
            YAxis = Magnitude * sin(Phi) * sin(Theta);
            ZAxis = Magnitude * cos(Phi);
        }

        /**
         * Calculate the magnitude and directions of the vector based off of its cartesian components
         */
        void calcSpherical() {
            double mag = sqrt(pow(XAxis, 2) + pow(YAxis, 2) + pow(ZAxis, 2));
            
            if (std::is_integral<ArithType>::value) {Magnitude = round(mag);}
            else {Magnitude = mag;}

            Theta = acos(ZAxis / mag);
            Phi = acos(XAxis / sqrt(pow(XAxis, 2) + pow(YAxis, 2)));
        }

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

            calcSpherical();
        }

        /**
         * Parametric Constructor - Input spherical coordinates
         * 
         * @param magnitude The magnitude of the vector
         * @param theta The polar angle of the vector (angle perpendicular to xy plane)
         * @param phi The azimuth angle of the vector (angle perpendicular to the xz plane)
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
            calcCartesian();
        }

        /**
         * Copy Constructor
         * 
         * @param vector A Vector_3D of a matching data type to be copied from
         */
        Vector_3D(const Vector_3D<ArithType> &vector) {
            XAxis = vector.getX();
            YAxis = vector.getY();
            ZAxis = vector.getZ();
            Magnitude = vector.getMag();
            Theta = vector.getTheta();
            Phi = vector.getPhi();
        }

        ArithType getX() const {return XAxis;}
        ArithType getY() const {return YAxis;}
        ArithType getZ() const {return ZAxis;}
        ArithType getMag() const {return Magnitude;}
        double getTheta() const {return Theta;}
        double getPhi() const {return Phi;}

        Vector_3D<ArithType> setX(ArithType x) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(XAxis, YAxis, ZAxis);
            XAxis = x;
            return output;
        }
        Vector_3D<ArithType> setY(ArithType y) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(XAxis, YAxis, ZAxis);
            YAxis = y;
            return output;
        }
        Vector_3D<ArithType> setZ(ArithType z) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(XAxis, YAxis, ZAxis);
            ZAxis = z;
            return output;
        }
        Vector_3D<ArithType> setMag(ArithType mag) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(Magnitude, Theta, Phi, false);
            Magnitude = mag;
            return output;
        }
        Vector_3D<ArithType> setTheta(double theta, bool degrees = true) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(Magnitude, Theta, Phi, false);
            Theta = theta * (degrees ? M_PI / 180 ; 1);
            return output;
        }
        Vector_3D<ArithType> setPhi(double phi, bool degrees = true) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(Magnitude, Theta, Phi, false);
            Phi = phi * (degrees ? M_PI / 180 ; 1);
            return output;
        }

        std::string rToString(bool specifyPositive = false, short round = 3) const {return "(" + astr::toString(astr::round(XAxis, round), specifyPositive) + ", " + astr::toString(astr::round(YAxis, round), specifyPositive) + ", " + astr::toString(astr::round(ZAxis, round), specifyPositive) + ")";}
        std::string sToString(bool specifyPositive = false, short round = 3, bool degrees = true) const {return "(" + astr::toString(astr::round(Magnitude), specifyPositive) + ", " + astr::toString(astr::round(Theta * (degrees ? 180 / M_PI : 1)), specifyPositive) + ", " + astr::toString(astr::round(Phi * (degrees ? 180 / M_PI : 1)), specifyPositive) + ")";}
        std::string rToString_Places(unsigned int beforeDecimal, unsigned int afterDecimal = 0, bool add = false, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Places(astr::round(XAxis, round), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(YAxis, round), specifyPositive) + ", " + astr::toString_Places(astr::round(ZAxis, round), specifyPositive) + ")";}
        std::string sToString_Places(unsigned int beforeDecimal, unsigned int afterDecimal = 0, bool add = false, bool specifyPositive = false, short round = 3, bool degrees = true) const {return "(" + astr::toString_Places(astr::round(Magnitude), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(Theta * (degrees ? 180 / M_PI : 1)), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(Phi * (degrees ? 180 / M_PI : 1)), beforeDecimal, afterDecimal, add, specifyPositive) + ")";}
        std::string rToString_Length(unsigned int length, bool leading = true, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Length(astr::round(XAxis, round), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(YAxis, round), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(ZAxis, round), length, leading, specifyPositive) + ")";}
        std::string sToString_Length(unsigned int length, bool leading = true, bool specifyPositive = false, short round = 3, bool degrees = true) const {return "(" + astr::toString_Length(astr::round(Magnitude), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(Theta * (degrees ? 180 / M_PI : 1)), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(Phi * (degrees ? 180 / M_PI : 1)), length, leading, specifyPositive) + ")";}
        std::string rToString_Sci(short decimals = 2, short exponentDigits = 1, bool e = true, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Sci(XAxis, decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(YAxis, decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(ZAxis, decimals, exponentDigits, e, specifyPositive) + ")";}
        std::string sToString_Sci(short decimals = 2, short exponentDigits = 1, bool e = true, bool specifyPositive = false, short round = 3, bool degrees = true) const {return "(" + astr::toString_Sci(Magnitude, decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(astr::round(Theta * (degrees ? 180 / M_PI : 1)), decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(astr::round(Phi * (degrees ? 180 / M_PI : 1)), decimals, exponentDigits, e, specifyPositive) + ")";}


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
