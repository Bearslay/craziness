#ifndef PHYSICS
#define PHYSICS

#include "astr.hpp"

#include <string>
#include <utility>
#include <cmath>
#include <vector>

#define COORD_RELATE_COMMON 0
#define COORD_RELATE_XCOORD 1
#define COORD_RELATE_YCOORD 2
#define COORD_RELATE_ZCOORD 3
#define COORD_RELATE_EUCLID 4
#define COORD_RELATE_TAXICB 5

/**
 * Collection of constants that can be used in various equations (and changed too in some cases)
 */
struct {
    /**
     * Acceleration due to gravity
     */
    double g = 9.80665;
} Constants;

template <typename ArithType> class Coord_3D {
    private:
        ArithType X = 0;
        ArithType Y = 0;
        ArithType Z = 0;
        unsigned char RelationMetric = COORD_RELATE_COMMON;

    public:
        Coord_3D() {
            static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        }
        Coord_3D(ArithType x, ArithType y, ArithType z) {
            static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
            X = x;
            Y = y;
            Z = z;
        }
        Coord_3D(const Coord_3D<ArithType> &coord) {
            X = coord.X;
            Y = coord.Y;
            Z = coord.Z;
        }

        ArithType getX() const {return X;}
        ArithType getY() const {return Y;}
        ArithType getZ() const {return Z;}

        Coord_3D<ArithType> setX(ArithType x) {
            Coord_3D<ArithType> output = Coord_3D<ArithType>(X, Y, Z);
            X = x;
            return output;
        }
        Coord_3D<ArithType> setY(ArithType y) {
            Coord_3D<ArithType> output = Coord_3D<ArithType>(X, Y, Z);
            Y = y;
            return output;
        }
        Coord_3D<ArithType> setZ(ArithType z) {
            Coord_3D<ArithType> output = Coord_3D<ArithType>(X, Y, Z);
            Z = z;
            return output;
        }

        std::vector<ArithType> toVector() const {return {X, Y, Z};}
        std::string toString(bool specifyPositive = false, short round = 3) const {return "(" + astr::toString(astr::round(X, round), specifyPositive) + ", " + astr::toString(astr::round(Y, round), specifyPositive) + ", " + astr::toString(astr::round(Z, round), specifyPositive) + ")";}
        std::string toString_Places(unsigned int beforeDecimal, unsigned int afterDecimal = 0, bool add = false, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Places(astr::round(X, round), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(Y, round), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(Z, round), beforeDecimal, afterDecimal, add, specifyPositive) + ")";}
        std::string toString_Length(unsigned int length, bool leading = true, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Length(astr::round(X, round), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(Y, round), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(Z, round), length, leading, specifyPositive) + ")";}
        std::string toString_Sci(short decimals = 2, short exponentDigits = 1, bool e = true, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Sci(X, decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(Y, decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(Z, decimals, exponentDigits, e, specifyPositive) + ")";}

        Coord_3D<ArithType> operator ! () const {return Coord_3D<ArithType>(-X, -Y, -Z);}
        bool operator == (const Coord_3D<ArithType> &coord) {
            switch (RelationMetric) {
                default:
                case COORD_RELATE_COMMON:
                    return X == coord.X && Y == coord.Y && Z == coord.Z;
                case COORD_RELATE_XCOORD:
                    return X == coord.X;
                case COORD_RELATE_YCOORD:
                    return Y == coord.Y;
                case COORD_RELATE_ZCOORD:
                    return Z == coord.Z;
                case COORD_RELATE_EUCLID:
                    return euclideanDistance() == coord.euclideanDistance();
                case COORD_RELATE_TAXICB:
                    return taxicabDistance() == coord.taxicabDistance();
            }
        }
        bool operator != (const Coord_3D<ArithType> &coord) {return !(Coord_3D<ArithType>(X, Y, Z) == coord);}
        bool operator < (const Coord_3D<ArithType> &coord) const {
            switch (RelationMetric) {
                default:
                case COORD_RELATE_COMMON:
                case COORD_RELATE_TAXICB:
                    return taxicabDistance() < coord.taxicabDistance();
                case COORD_RELATE_XCOORD:
                    return X < coord.X;
                case COORD_RELATE_YCOORD:
                    return Y < coord.Y;
                case COORD_RELATE_ZCOORD:
                    return Z < coord.Z;
                case COORD_RELATE_EUCLID:
                    return euclideanDistance() < coord.euclideanDistance();
            }
        }
        bool operator <= (const Coord_3D<ArithType> &coord) const {
            switch (RelationMetric) {
                default:
                case COORD_RELATE_COMMON:
                case COORD_RELATE_TAXICB:
                    return taxicabDistance() <= coord.taxicabDistance();
                case COORD_RELATE_XCOORD:
                    return X <= coord.X;
                case COORD_RELATE_YCOORD:
                    return Y <= coord.Y;
                case COORD_RELATE_ZCOORD:
                    return Z <= coord.Z;
                case COORD_RELATE_EUCLID:
                    return euclideanDistance() <= coord.euclideanDistance();
            }
        }
        bool operator > (const Coord_3D<ArithType> &coord) const {
            switch (RelationMetric) {
                default:
                case COORD_RELATE_TAXICB:
                case COORD_RELATE_TAXICB:
                    return taxicabDistance() > coord.taxicabDistance();
                case COORD_RELATE_XCOORD:
                    return X > coord.X;
                case COORD_RELATE_YCOORD:
                    return Y > coord.Y;
                case COORD_RELATE_ZCOORD:
                    return Z > coord.Z;
                case COORD_RELATE_EUCLID:
                    return euclideanDistance() > coord.euclideanDistance();
            }
        }
        bool operator >= (const Coord_3D<ArithType> &coord) const {
            switch (RelationMetric) {
                default:
                case COORD_RELATE_COMMON:
                case COORD_RELATE_TAXICB:
                    return taxicabDistance() >= coord.taxicabDistance();
                case COORD_RELATE_XCOORD:
                    return X >= coord.X;
                case COORD_RELATE_YCOORD:
                    return Y >= coord.Y;
                case COORD_RELATE_ZCOORD:
                    return Z >= coord.Z;
                case COORD_RELATE_EUCLID:
                    return euclideanDistance() >= coord.euclideanDistance();
            }
        }
        unsigned char setRelationMetric(unsigned char metric = COORD_RELATE_COMMON) {
            unsigned char output = RelationMetric;
            RelationMetric = (6 + metric % 6) % 6;
            return output;
        }
        unsigned char getRelationMetric() const {return RelationMetric;}

        Coord_3D<ArithType> operator += (const Coord_3D<ArithType> &coord) {
            Coord_3D<ArithType> output = Coord_3D<ArithType>(X, Y, Z);
            X += coord.X;
            Y += coord.Y;
            Z += coord.Z;
            return output;
        }
        Coord_3D<ArithType> operator -= (const Coord_3D<ArithType> &coord) {
            Coord_3D<ArithType> output = Coord_3D<ArithType>(X, Y, Z);
            X -= coord.X;
            Y -= coord.Y;
            Z -= coord.Z;
            return output;
        }
        Coord_3D<ArithType> operator *= (const ArithType &scalar) {
            Coord_3D<ArithType> output = Coord_3D<ArithType>(X, Y, Z);
            X *= scalar;
            Y *= scalar;
            Z *= scalar;
            return output;
        }
        Coord_3D<ArithType> operator + (const Coord_3D<ArithType> &coord) const {return Coord_3D<ArithType>(X + coord.X, Y + coord.Y, Z + coord.Z);}
        Coord_3D<ArithType> operator - (const Coord_3D<ArithType> &coord) const {return Coord_3D<ArithType>(X - coord.X, Y - coord.Y, Z - coord.Z);}
        Coord_3D<ArithType> operator * (const ArithType &scalar) const {return Coord_3D<ArithType>(X * scalar, Y * scalar, Z * scalar);}

        ArithType euclideanDistance(const Coord_3D<ArithType> &coord = Coord_3D<ArithType>(0, 0, 0)) const {
            double distance = std::sqrt(std::pow(X - coord.X, 2) + std::pow(Y - coord.Y, 2) + std::pow(Z - coord.Z, 2));
            return std::is_integral<ArithType>::value ? std::round(distance) : distance;
        }
        ArithType taxicabDistance(const Coord_3D<ArithType> &coord = Coord_3D<ArithType>(0, 0, 0)) const {return std::fabs(X - coord.X) + std::fabs(Y - coord.Y) + std::fabs(Z - coord.Z);}
};

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
        ArithType X = 0;
        /**
         * The y-component of the vector with respect to the origin (0, 0, 0)
         */
        ArithType Y = 0;
        /**
         * The z-component of the vector with respect to the origin (0, 0, 0)
         */
        ArithType Z = 0;
        
        /**
         * The magnitude of the vector
         */
        ArithType M = 0;
        /**
         * The polar angle (Θ) within a spherical coordinate system stored in radians
         * This context is using the physics convention and would instead be phi using the mathematics conventions
         * 
         * When limited to integral data types, angles get real inaccurate for vectors of increasing magnitudes, so is always a double to account for that
         */
        double T = 0;
        /**
         * The azimuth angle (φ) within a spherical coordinate system stored in radians
         * This context is using the physics convention and would instead be theta using the mathematics conventions
         * 
         * When limited to integral data types, angles get real inaccurate for vectors of increasing magnitudes, so is always a double to account for that
         */
        double P = 0;

        /**
         * Calculate the cartesian components of the vector based off of its magnitude and directions
         */
        void calcCartesian() {
            if (std::is_integral<ArithType>::value) {
                X = std::round(M * sin(P) * cos(T));
                Y = std::round(M * sin(P) * sin(T));
                Z = std::round(M * cos(P));
                return;
            }

            X = M * sin(P) * cos(T);
            Y = M * sin(P) * sin(T);
            Z = M * cos(P);
        }

        /**
         * Calculate the magnitude and directions of the vector based off of its cartesian components
         */
        void calcSpherical() {
            double m = sqrt(pow(X, 2) + pow(Y, 2) + pow(Z, 2));
            
            if (std::is_integral<ArithType>::value) {M = round(m);}
            else {M = m;}

            T = acos(Z / m);
            P = acos(X / sqrt(pow(X, 2) + pow(Y, 2)));
        }

    public:
        /**
         * Default Constructor (set all values to zero)
         */
        Vector_3D() {
            static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
        }
        /**
         * Parametric Constructor - Input cartesian coordinates
         * 
         * @param xaxis The x-component of the vector
         * @param yaxis The y-component of the vector
         * @param zaxis The z-component of the vector
         */
        Vector_3D(ArithType x, ArithType y, ArithType z) {
            static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");
            X = x;
            Y = y;
            Z = z;
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

            M = magnitude;
            T = theta;
            P = phi;

            if (degrees) {
                T *= M_PI / 180;
                P *= M_PI / 180;
            }
            calcCartesian();
        }
        /**
         * Copy Constructor
         * 
         * @param vector A Vector_3D of a matching data type to be copied from
         */
        Vector_3D(const Vector_3D<ArithType> &vector) {
            X = vector.X;
            Y = vector.Y;
            Z = vector.Z;
            M = vector.M;
            T = vector.T;
            P = vector.P;
        }

        ArithType getX() const {return X;}
        ArithType getY() const {return Y;}
        ArithType getZ() const {return Z;}
        ArithType getMag() const {return M;}
        double getTheta() const {return T;}
        double getPhi() const {return P;}

        Vector_3D<ArithType> setX(ArithType x) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(X, Y, Z);
            X = x;
            return output;
        }
        Vector_3D<ArithType> setY(ArithType y) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(X, Y, Z);
            Y = y;
            return output;
        }
        Vector_3D<ArithType> setZ(ArithType z) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(X, Y, Z);
            Z = z;
            return output;
        }
        Vector_3D<ArithType> setMag(ArithType magnitude) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(M, T, P, false);
            M = magnitude;
            return output;
        }
        Vector_3D<ArithType> setTheta(double theta, bool degrees = true) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(M, T, P, false);
            T = theta * (degrees ? M_PI / 180 : 1);
            return output;
        }
        Vector_3D<ArithType> setPhi(double phi, bool degrees = true) {
            Vector_3D<ArithType> output = Vector_3D<ArithType>(M, T, P, false);
            P = phi * (degrees ? M_PI / 180 : 1);
            return output;
        }

        std::string rToString(bool specifyPositive = false, short round = 3) const {return "(" + astr::toString(astr::round(X, round), specifyPositive) + ", " + astr::toString(astr::round(Y, round), specifyPositive) + ", " + astr::toString(astr::round(Z, round), specifyPositive) + ")";}
        std::string sToString(bool specifyPositive = false, short round = 3, bool degrees = true) const {return "(" + astr::toString(astr::round(M), specifyPositive) + ", " + astr::toString(astr::round(T * (degrees ? 180 / M_PI : 1)), specifyPositive) + ", " + astr::toString(astr::round(P * (degrees ? 180 / M_PI : 1)), specifyPositive) + ")";}
        std::string rToString_Places(unsigned int beforeDecimal, unsigned int afterDecimal = 0, bool add = false, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Places(astr::round(X, round), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(Y, round), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(Z, round), beforeDecimal, afterDecimal, add, specifyPositive) + ")";}
        std::string sToString_Places(unsigned int beforeDecimal, unsigned int afterDecimal = 0, bool add = false, bool specifyPositive = false, short round = 3, bool degrees = true) const {return "(" + astr::toString_Places(astr::round(M), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(T * (degrees ? 180 / M_PI : 1)), beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(astr::round(P * (degrees ? 180 / M_PI : 1)), beforeDecimal, afterDecimal, add, specifyPositive) + ")";}
        std::string rToString_Length(unsigned int length, bool leading = true, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Length(astr::round(X, round), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(Y, round), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(Z, round), length, leading, specifyPositive) + ")";}
        std::string sToString_Length(unsigned int length, bool leading = true, bool specifyPositive = false, short round = 3, bool degrees = true) const {return "(" + astr::toString_Length(astr::round(M), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(T * (degrees ? 180 / M_PI : 1)), length, leading, specifyPositive) + ", " + astr::toString_Length(astr::round(P * (degrees ? 180 / M_PI : 1)), length, leading, specifyPositive) + ")";}
        std::string rToString_Sci(short decimals = 2, short exponentDigits = 1, bool e = true, bool specifyPositive = false, short round = 3) const {return "(" + astr::toString_Sci(X, decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(Y, decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(Z, decimals, exponentDigits, e, specifyPositive) + ")";}
        std::string sToString_Sci(short decimals = 2, short exponentDigits = 1, bool e = true, bool specifyPositive = false, short round = 3, bool degrees = true) const {return "(" + astr::toString_Sci(M, decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(astr::round(T * (degrees ? 180 / M_PI : 1)), decimals, exponentDigits, e, specifyPositive) + ", " + astr::toString_Sci(astr::round(P * (degrees ? 180 / M_PI : 1)), decimals, exponentDigits, e, specifyPositive) + ")";}

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
