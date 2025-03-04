#ifndef HEX
#define HEX

#include "astr.hpp"

// #include <unordered_map>
#include <vector>
#include <string>
#include <cmath>

#define HEX_DIR_000 0
#define HEX_DIR_030 1
#define HEX_DIR_060 2
#define HEX_DIR_090 3
#define HEX_DIR_120 4
#define HEX_DIR_150 5
#define HEX_DIR_180 6
#define HEX_DIR_210 7
#define HEX_DIR_240 8
#define HEX_DIR_270 9
#define HDX_DIR_300 10
#define HDX_DIR_330 11

#define HEX_FLT true
#define HEX_PTY false

#define HEX_RELATE_COMMON 0
#define HEX_RELATE_ICOORD 1
#define HEX_RELATE_JCOORD 2
#define HEX_RELATE_KCOORD 3
#define HEX_RELATE_EUCLID 4
#define HEX_RELATE_TAXICB 5

namespace hex {
    /**
     * Contains the twelve cardinal/intermediate directions for hexagons (seperated by 30 degrees)
     * By default, the directions are arranged for flat hexagons (vertex at 0), but is adaptable for pointy ones with an index offset of 1 (30 degrees)
     */
    const char UnitDirections[12][3] = {
        {2, -1, -1},
        {1, -1, 0},
        {1, -2, 1},
        {0, -1, 1},
        {-1, -1, 2},
        {-1, 0, 1},
        {-2, 1, 1},
        {-1, 1, 0},
        {-1, 2, -1},
        {0, 1, -1},
        {1, 1, -2},
        {1, 0, -1}
    };

    /**
     * A coordinate that can be used to give hexagons in a 2d-plane positions
     * This particular system adapts cubic coordinates into a hexagonal system, so has three-axis that are offset by 60 degrees rather than the more typical 90 degrees
     * 
     * The components of any coordinate should add up to 0 (or close to zero with floating-point data types)
     * 
     * @tparam Type An arithmetic type used for the positioning of the coordinate, usually an integer or float depending on the application of HexCoord
     */
    template <typename Type> class HexCoord {
        private:
            /**
             * The i-component of the HexCoord
             * Using flat hexagons, the i-axis goes top/bottom
             * The value increases going 'right' and decreases going 'left'
             */
            Type I = 0;
            /**
             * The j-component of the HexCoord
             * Using flat hexagons, the j-axis goes top-left/bottom-right
             * The value increases going 'down' or 'left' and decreases going 'up' or 'right'
             */
            Type J = 0;
            /**
             * The k-component of the HexCoord
             * Using flat hexagons, the k-axis goes bottom-left/top-right
             * The value increases going 'up' or 'left' and decreases going 'down' or 'right'
             */
            Type K = 0;
            /*
            * The metric that the HexCoord will use to measure equality and relation
            * This is useful for sorting by different aspects of the HexCoord or just having the ability to use relational operators
            * 
            * HEX_RELATE_COMMON [0] = Most common relationship used for a particular relationship (varies per operator)
            * HEX_RELATE_ICOORD [1] = The value of the HexCoord's i-component
            * HEX_RELATE_JCOORD [2] = The value of the HexCoord's j-component
            * HEX_RELATE_KCOORD [3] = The value of the HexCoord's k-component
            * HEX_RELATE_EUCLID [4] = The euclidean distance of the HexCoord in relation to (0, 0, 0)
            * HEX_RELATE_TAXICB [5] = The taxicab distance of the HexCoord in relation to (0, 0, 0)
            */
            unsigned char RelationMetric = HEX_RELATE_COMMON;

        public:
            /**
             * Parametric constructor
             * 
             * Due to the nature of hexagonal coordinates, if 'i + j + k != 0' then only i and j retain their parameter-defined values while k is automatically calculated
             * This fail-safe is reminiscent of an alternate way to store hexagonal coordinates, where i and j are the only true variables and k is implicitly defined
             * 
             * @param i i-component of the coordinate
             * @param j j-component of the coordinate
             * @param k k-component of the coordinate
             */
            HexCoord(Type i, Type j, Type k) {
                static_assert(std::is_arithmetic<Type>::value, "Type must be an arithmetic type");
                I = i;
                J = j;
                K = I + J + K == 0 ? k : 0 - I - J;
            }
            /**
             * Default constructor - flat hexagon at (0, 0, 0)
             */
            HexCoord() {HexCoord(0, 0, 0);}
            /**
             * Copy constructor
             * 
             * @param coord HexCoord to copy from
             */
            HexCoord(const HexCoord<Type> &coord) {HexCoord(coord.I, coord.J, coord.K);}
            
            template <typename IntType> HexCoord<IntType> toIntHex(HexCoord<Type> &coord) const {
                if (std::is_integral<Type>::value || std::is_floating_point<IntType>::value) {return coord;}

                IntType i = std::round(coord.I);
                IntType j = std::round(coord.J);
                IntType k = std::round(coord.K);

                double iDiff = std::fabs(i - coord.I);
                double jDiff = std::fabs(j - coord.J);
                double kDiff = std::fabs(k - coord.K);

                if (iDiff > jDiff && iDiff > kDiff) {
                    i = -j - k;
                } else if (jDiff > kDiff) {
                    j = -i - k;
                } else {
                    k = -i - j;
                }

                return HexCoord<IntType>(i, j, k);
            }

            /**
             * Get the Hecoord's i-component
             * 
             * @returns The i-component
             */
            Type getI() const {return I;}
            /**
             * Get the Hecoord's j-component
             * 
             * @returns The j-component
             */
            Type getJ() const {return J;}
            /**
             * Get the Hecoord's k-component
             * 
             * @returns The k-component
             */
            Type getK() const {return K;}

            Type getY(Type radius = 1) const {
                double y = sqrt(3) * (I / 2 + J) * radius;
                return std::is_integral<Type>::value ? (y) : y;
            }
            Type getX(Type radius = 1) const {
                double x = I * (3 / 2) * radius;
                return std::is_integral<Type>::value ? std::round(x) : x;
            }

            int getY_npp(int dimy) const {return J * dimy + (I < 0 ? -1 : 1) * (I % 2 != 0 ? dimy / 2 : 0) + I / 2 * dimy;}
            int getX_npp(int dimx, int inset) const {return I * dimx - I * inset;}

            /**
             * Get each component arranged in a vector
             * 
             * @returns A vector containing the three components of the HexCoord {i, j, k}
             */
            std::vector<Type> toVector() const {return {I, J, K};}
            /**
             * Get the HexCoord as a string formatted as (i, j, k)
             * 
             * @returns The HexCoord formatted as a string
             */
            std::string toString(bool specifyPositive = false) const {return "(" + astr::toString(I, specifyPositive) + ", " + astr::toString(J, specifyPositive) + ", " + astr::toString(K, specifyPositive) + ")";}
            /**
             * Get the HexCoord as a string formatted as (i, j, k)
             * This particular overload allows for extra formatting with leading/trailing zeros
             * 
             * @param beforeDecimal The amount of digits to place before the decimal point
             * @param afterDecimal The amount of digits to place after the decimal point (even for integral types)
             * @param add If true, beforeDecimal/afterDecimal switch from amount of digits to amount of leading/trailing zeros (regardless of digits present already)
             * @returns The HexCoord formatted as a string
             */
            std::string toString_Places(unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) const {return "(" + astr::toString_Places(I, beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(J, beforeDecimal, afterDecimal, add, specifyPositive) + ", " + astr::toString_Places(K, beforeDecimal, afterDecimal, add, specifyPositive) + ")";}
            /**
             * Get the HexCoord as a string formatted as (i, j, k)
             * This particular overload allows for extra formatting by specifying a length for each component
             * 
             * @param length The desired width of the output string
             * @param leading Whether to add the extra zeros to the beginning or end of the string
             * @returns The HexCoord formatted as a string
             */
            std::string toString_Length(unsigned long length, bool leading = true, bool specifyPositive = false) const {return "(" + astr::toString_Length(I, length, leading, specifyPositive) + ", " + astr::toString_Length(J, length, leading, specifyPositive) + ", " + astr::toString_Length(K, length, leading, specifyPositive) + ")";}

            /**
             * Get the HexCoord's position if rotated 180 degrees with respect to the origin
             * 
             * @returns The HexCoord's position if rotated 180 degrees with respect to the origin
             */
            HexCoord<Type> operator ! () const {return HexCoord<Type>(-I, -J, -K);}
            /**
             * Checks for equality between two HexCoords
             * 
             * Depending on what metric the HexCoord is using to measure relations, the output will differ
             * In this case, HEX_RELATE_COMMON simply combines HEX_RELATE_ICOORD, HEX_RELATE_JCOORD, and HEX_RELATE_KCOORD
             * 
             * @param coord HexCoord to compare with
             * @returns true for equality or false otherwise
             */
            bool operator == (const HexCoord<Type> &coord) const {
                switch (RelationMetric) {
                    default:
                    case HEX_RELATE_COMMON:
                        return I == coord.I && J == coord.J && K == coord.K;
                    case HEX_RELATE_ICOORD:
                        return I == coord.I;
                    case HEX_RELATE_JCOORD:
                        return J == coord.J;
                    case HEX_RELATE_KCOORD:
                        return K == coord.K;
                    case HEX_RELATE_EUCLID:
                        return euclideanDistance() == coord.euclideanDistance();
                    case HEX_RELATE_TAXICB:
                        return taxicabDistance() == coord.taxicabDistance();
                }
            }
            /**
             * Checks for inequality between two HexCoords
             * 
             * Depending on what metric the HexCoord is using to measure relations, the output will differ
             * In this case, HEX_RELATE_COMMON simply combines HEX_RELATE_ICOORD, HEX_RELATE_JCOORD, and HEX_RELATE_KCOORD
             * 
             * @param coord HexCoord to compare with
             * @returns true for unequality or false otherwise
             */
            bool operator != (const HexCoord<Type> &coord) const {return !(HexCoord<Type>(I, J, K) == coord);}
            /**
             * Checks for whether a HexCoord is less than another HexCoord
             * 
             * Depending on what metric the HexCoord is using to measure relations, the output will differ
             * In this case, HEX_RELATE_COMMON evaluates the same as HEX_RELATE_TAXICB
             * 
             * @param coord HexCoord to compare with
             * @returns true if the HexCoord is less than the other or false otherwise
             */
            bool operator < (const HexCoord<Type> &coord) const {
                switch (RelationMetric) {
                    default:
                    case HEX_RELATE_COMMON:
                    case HEX_RELATE_TAXICB:
                        return taxicabDistance() < coord.taxicabDistance();
                    case HEX_RELATE_ICOORD:
                        return I < coord.I;
                    case HEX_RELATE_JCOORD:
                        return J < coord.J;
                    case HEX_RELATE_KCOORD:
                        return K < coord.K;
                    case HEX_RELATE_EUCLID:
                        return euclideanDistance() < coord.euclideanDistance();
                }
            }
            /**
             * Checks for whether a HexCoord is less than or equal to another HexCoord
             * 
             * Depending on what metric the HexCoord is using to measure relations, the output will differ
             * In this case, HEX_RELATE_COMMON evaluates the same as HEX_RELATE_TAXICB
             * 
             * @param coord HexCoord to compare with
             * @returns true if the HexCoord is less than or equal to the other or false otherwise
             */
            bool operator <= (const HexCoord<Type> &coord) const {
                switch (RelationMetric) {
                    default:
                    case HEX_RELATE_COMMON:
                    case HEX_RELATE_TAXICB:
                        return taxicabDistance() <= coord.taxicabDistance();
                    case HEX_RELATE_ICOORD:
                        return I <= coord.I;
                    case HEX_RELATE_JCOORD:
                        return J <= coord.J;
                    case HEX_RELATE_KCOORD:
                        return K <= coord.K;
                    case HEX_RELATE_EUCLID:
                        return euclideanDistance() <= coord.euclideanDistance();
                }
            }
            /**
             * Checks for whether a HexCoord is greater than another HexCoord
             * 
             * Depending on what metric the HexCoord is using to measure relations, the output will differ
             * In this case, HEX_RELATE_COMMON evaluates the same as HEX_RELATE_TAXICB
             * 
             * @param coord HexCoord to compare with
             * @returns true if the HexCoord is greater than the other or false otherwise
             */
            bool operator > (const HexCoord<Type> &coord) const {
                switch (RelationMetric) {
                    default:
                    case HEX_RELATE_COMMON:
                    case HEX_RELATE_TAXICB:
                        return taxicabDistance() > coord.taxicabDistance();
                    case HEX_RELATE_ICOORD:
                        return I > coord.I;
                    case HEX_RELATE_JCOORD:
                        return J > coord.J;
                    case HEX_RELATE_KCOORD:
                        return K > coord.K;
                    case HEX_RELATE_EUCLID:
                        return euclideanDistance() > coord.euclideanDistance();
                }
            }
            /**
             * Checks for whether a HexCoord is greater than or equal to another HexCoord
             * 
             * Depending on what metric the HexCoord is using to measure relations, the output will differ
             * In this case, HEX_RELATE_COMMON evaluates the same as HEX_RELATE_TAXICB
             * 
             * @param coord HexCoord to compare with
             * @returns true if the HexCoord is greater than or equal to the other or false otherwise
             */
            bool operator >= (const HexCoord<Type> &coord) const {
                switch (RelationMetric) {
                    default:
                    case HEX_RELATE_COMMON:
                    case HEX_RELATE_TAXICB:
                        return taxicabDistance() >= coord.taxicabDistance();
                    case HEX_RELATE_ICOORD:
                        return I >= coord.I;
                    case HEX_RELATE_JCOORD:
                        return J >= coord.J;
                    case HEX_RELATE_KCOORD:
                        return K >= coord.K;
                    case HEX_RELATE_EUCLID:
                        return euclideanDistance() >= coord.euclideanDistance();
                }
            }
            /**
             * Set the metric that the HexCoord will use to measure equality and relation
             * This is useful for sorting by different aspects of the HexCoord or just having the ability to use relational operators
             * 
             * HEX_RELATE_COMMON [0] = Most common relationship used for a particular relationship (varies per operator)
             * HEX_RELATE_ICOORD [1] = The value of the HexCoord's i-component
             * HEX_RELATE_JCOORD [2] = The value of the HexCoord's j-component
             * HEX_RELATE_KCOORD [3] = The value of the HexCoord's k-component
             * HEX_RELATE_EUCLID [4] = The euclidean distance of the HexCoord in relation to (0, 0, 0)
             * HEX_RELATE_TAXICB [5] = The taxicab distance of the HexCoord in relation to (0, 0, 0)
             * 
             * @param metric The new metric to be used (expressed as a number 0-5)
             * @returns The old metric used (expressed as a number 0-5)
             */
            unsigned char setRelationMetric(unsigned char metric = HEX_RELATE_COMMON) {
                unsigned char output = RelationMetric;
                RelationMetric = (6 + metric % 6) % 6;
                return output;
            }
            /**
             * Get the metric that the HexCoord uses to measure equality and relation
             * This is useful for sorting by different aspects of the HexCoord or just having the ability to use relational operators
             * 
             * HEX_RELATE_COMMON [0] = Most common relationship used for a particular relationship (varies per operator)
             * HEX_RELATE_ICOORD [1] = The value of the HexCoord's i-component
             * HEX_RELATE_JCOORD [2] = The value of the HexCoord's j-component
             * HEX_RELATE_KCOORD [3] = The value of the HexCoord's k-component
             * HEX_RELATE_EUCLID [4] = The euclidean distance of the HexCoord in relation to (0, 0, 0)
             * HEX_RELATE_TAXICB [5] = The taxicab distance of the HexCoord in relation to (0, 0, 0)
             * 
             * @returns The old metric used (expressed as a number 0-5)
             */
            unsigned char getRelationMetric() const {return RelationMetric;}

            /**
             * Add the position of a HexCoord to the current one
             * This is essentially vector addition
             * 
             * @param coord HexCoordinate to add to the current one
             * @returns The original position of the HexCoord
             */
            HexCoord<Type> operator += (const HexCoord<Type> &coord) {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);
                I += coord.I;
                J += coord.J;
                K += coord.K;
                return output;
            }
            /**
             * Subtract the position of a HexCoord from the current one
             * This is essentially vector subtraction
             * 
             * @param coord HexCoordinate to subtract from the current one
             * @returns The original position of the HexCoord
             */
            HexCoord<Type> operator -= (const HexCoord<Type> &coord) {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);
                I -= coord.I;
                J -= coord.J;
                K -= coord.K;
                return output;
            }
            /**
             * Multiply the position of the HexCoord by a scalar value
             * This is essentially vector multiplication
             * 
             * @param scalar Scalar value of an arithmetic type matching the HexCoord to multiply the current HexCoord by
             * @returns The original position of the HexCoord
             */
            HexCoord<Type> operator *= (const Type &scalar) {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);
                I *= scalar;
                J *= scalar;
                K *= scalar;
                return output;
            }
            /**
             * Add the position of two HexCoords
             * This is essentially vector addition
             * 
             * @param coord HexCoordinate to add to the current one
             * @returns The position of the resulting HexCoord
             */
            HexCoord<Type> operator + (const HexCoord<Type> &coord) const {return HexCoord<Type>(I + coord.I, J + coord.J, K + coord.K);}
            /**
             * Subtract the position of two HexCoords
             * This is essentially vector subtraction
             * 
             * @param coord HexCoordinate to subtract from the current one
             * @returns The position of the resulting HexCoord
             */
            HexCoord<Type> operator - (const HexCoord<Type> &coord) const {return HexCoord<Type>(I - coord.I, J - coord.J, K - coord.K);}
            /**
             * Multiply the position of a HexCoord by a scalar value
             * This is essentially vector multiplication
             * 
             * @param scalar Scalar value of an arithmetic type matching the HexCoord to multiply the current HexCoord by
             * @returns The position of the resulting HexCoord
             */
            HexCoord<Type> operator * (const Type &scalar) const {return HexCoord<Type>(I * scalar, J * scalar, K * scalar);}

            /**
             * Get the euclidean distance between two HexCoords
             * By default, this function gets the distance from (0, 0, 0)
             * 
             * @param coord HexCoord to find distance between
             * @returns An arithmetic type (from the current HexCoord) with the euclidean distance between the specified coordinates
             */
            Type euclideanDistance(const HexCoord<Type> &coord = HexCoord<Type>(0, 0, 0)) const {
                double distance = sqrt(std::pow(I - coord.I, 2) + std::pow(J - coord.J, 2) + std::pow(K - coord.K, 2));
                return std::is_integral<Type>::value ? std::round(distance) : distance;
            }
            /**
             * Get the taxicab distance between two HexCoords
             * By default, this function gets the distance from (0, 0, 0)
             * 
             * @param coord HexCoord to find distance between
             * @returns An arithmetic type (from the current HexCoord) with the taxicab distance between the specified coordinates
             */
            Type taxicabDistance(const HexCoord<Type> &coord = HexCoord<Type>(0, 0, 0)) const {return (std::fabs(I - coord.I) + std::fabs(J - coord.J) + std::fabs(K - coord.K)) / 2;}

            /**
             * Get the HexCoord located in the direction specified
             * The adjacent directions are considered the six directions that match up with the edges of the hexagon, not the vertices
             * 
             * Out of the twelve directions available, the odd ones should be used (direction * 60 + 30)
             * 
             * @param direction The direction to project in
             * @param offset A scalar used to get the location of further coordinates
             * @returns The HexCoord found at the offset in the specified direction
             */
            HexCoord<Type> getAdjacent(unsigned char direction, Type offset = 1) const {
                direction = std::abs(direction) % 6 * 2 + 1;
                return HexCoord<Type>(I + UnitDirections[direction][0] * std::fabs(offset), J + UnitDirections[direction][1] * std::fabs(offset), K + UnitDirections[direction][2] * std::fabs(offset));
            }
            /**
             * Get the HexCoord located in the direction specified
             * The diagonal directions are considered the six directions that match up with the vertices of the hexagon, not the edges
             * 
             * Out of the twelve directions available, the even ones should be used (direction * 60)
             * 
             * @param direction The direction to project in
             * @param offset A scalar used to get the location of further coordinates (will have a somewhat dramatic effect due to how diagonals work with hexagons)
             * @returns The HexCoord found at the offset in the specified direction
             */
            HexCoord<Type> getDiagonal(unsigned char direction, Type offset = 1) const {
                direction = std::abs(direction) % 6 * 2;
                return HexCoord<Type>(I + UnitDirections[direction][0] * std::fabs(offset), J + UnitDirections[direction][1] * std::fabs(offset), K + UnitDirections[direction][2] * std::fabs(offset));
            }
            /**
             * Get the HexCoord located in the direction specified
             * This function combines getAdjacent() and getDiagonal(), so any of the twelve directions work
             * 
             * @param direction The direction to project in
             * @param offset A scalar used to get the location of further coordinates (will have a somewhat dramatic effect when direction is even due to how diagonals work with hexagons)
             * @returns The HexCoord found at the offset in the specified direction
             */
            HexCoord<Type> getNeighbor(unsigned char direction, Type offset = 1) const {
                direction = std::abs(direction) % 12;
                return direction % 2 == 0 ? getDiagonal(direction, offset) : getAdjacent(direction, offset);
            }

            HexCoord<Type> move(unsigned char direction, Type distance = 1) {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);

                direction = std::abs(direction) % 12 * 2;
                HexCoord<Type> newpos = direction % 2 == 0 ? getDiagonal(direction, distance) : getAdjacent(direction, distance);
                I = newpos.getI();
                J = newpos.getJ();
                K = newpos.getK();

                return output;
            }
            HexCoord<Type> rotate(int steps = 1) {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);

                if (abs(steps) % 6 == 0) {return output;}

                HexCoord<Type> newpos = HexCoord<Type>(I, J, K);
                for (int i = 0; i < std::abs(steps); i++) {
                    newpos = steps > 0 ? HexCoord<Type>(-newpos.getK(), -newpos.getI(), -newpos.getJ()) : HexCoord<Type>(-newpos.getJ(), -newpos.getK(), -newpos.getI());
                }
                I = newpos.getI();
                J = newpos.getJ();
                K = newpos.getK();

                return output;
            }
            HexCoord<Type> reflectI() {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);

                Type temp = J;
                J = K;
                K = temp;

                return output;
            }
            HexCoord<Type> reflectJ() {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);

                Type temp = I;
                I = K;
                K = temp;

                return output;
            }
            HexCoord<Type> reflectK() {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);

                Type temp = I;
                I = J;
                J = temp;

                return output;
            }
            HexCoord<Type> reflect(bool i, bool j, bool k) {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);

                if (i) {reflectI();}
                if (j) {reflectJ();}
                if (k) {reflectK();}

                return output;
            }
    };
}

#endif /* HEX */
