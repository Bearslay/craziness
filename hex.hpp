#ifndef HEX
#define HEX

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

/**
 * Contain custom toString() and toWideString() functions that mainly feature improved leading/trailing zero formatting
 */
namespace str {
    /**
     * Get an arithmetic data type as a string
     * 
     * @returns The inputted data type formatted as a string
     */
    template <typename Type> const std::string toString(Type input, bool specifyPositive = false) {
        // std::to_string() only works for arithmetic data types
        static_assert(std::is_arithmetic<Type>::value, "Type must be an arithmetic type");

        std::string output = std::to_string(input);
        // Remove any extra decimal places on the end
        if (std::is_floating_point<decltype(input)>::value) {
            for (unsigned long i = output.length() - 1; i >= 0; i--) {
                if (output[i] == '.') {
                    output.erase(i, 1);
                    break;
                } else if (output[i] != '0') {
                    break;
                }
                output.erase(i, 1);
            }
        }
        return (specifyPositive && input > 0 ? "+" : "") + output;
    }
    /**
     * Get an arithmetic data type as a string
     * This particular overload allows for extra formatting with leading/trailing zeros
     * 
     * @param beforeDecimal The amount of digits to place before the decimal point
     * @param afterDecimal The amount of digits to place after the decimal point (even for integral types)
     * @param add If true, beforeDecimal/afterDecimal switch from amount of digits to amount of leading/trailing zeros (regardless of digits present already)
     * @returns The inputted data type formatted as a string
     */
    template <typename Type> const std::string toString_Places(Type input, unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) {
        std::string output = str::toString(input, specifyPositive);
        unsigned long originalLength = output.length();

        if (beforeDecimal == 0 && afterDecimal == 0) {return output;}

        unsigned long beforeDecimalCount = 0;
        while (beforeDecimalCount < output.length()) {
            if (output[beforeDecimalCount] == '.') {break;}
            beforeDecimalCount++;
        }

        // Add leading zeros
        if (beforeDecimal > 0) {
            for (unsigned long i = add ? 0 : beforeDecimalCount; i < beforeDecimal; i++) {
                output.insert(input < 0 || specifyPositive ? 1 : 0, 1, '0');
            }

            if (afterDecimal == 0) {return output;}
        }

        unsigned long afterDecimalCount = output.find('.') == std::string::npos ? 0 : originalLength - beforeDecimalCount - 1;

        // Add trailing zeros
        if (afterDecimalCount == 0) {output += ".";}
        for (unsigned long i = add ? 0 : afterDecimalCount; i < afterDecimal; i++) {
            output += "0";
        }

        return output;
    }
    /**
     * Get an arithmetic data type as a string
     * This particular overload allows for extra formatting by specifying a length for each component
     * 
     * @param length The desired width of the output string
     * @param leading Whether to add the extra zeros to the beginning or end of the string
     * @returns The inputted data type formatted as a string
     */
    template <typename Type> const std::string toString_Length(Type input, unsigned long length, bool leading = true, bool specifyPositive = false) {
        std::string output = str::toString(input, specifyPositive);
        bool hasDecimal = output.find('.') != std::string::npos;

        unsigned long beforeDecimal = 0;
        while (beforeDecimal < output.length()) {
            if (output[beforeDecimal] == '.') {break;}
            beforeDecimal++;
        }
        unsigned long afterDecimal = hasDecimal ? output.length() - beforeDecimal - 1 : 0;

        if (leading) {return str::toString_Places(input, length - afterDecimal - (hasDecimal ? 1 : 0), 0, false, specifyPositive);}
        return str::toString_Places(input, 0, length - beforeDecimal - 1, false, specifyPositive);
    }
    /**
     * Get a string as a wide string
     * 
     * @returns The inputted string as a wide string
     */
    const std::wstring toWideString(std::string input) {
        std::wstring output;
        for (unsigned long i = 0; i < input.length(); i++) {
            output += input[i];
        }
        return output;
    }
    /**
     * Get an arithmetic data type as a wide string
     * 
     * @returns The inputted data type formatted as a wide string
     */
    template <typename Type> const std::wstring toWideString(Type input, bool specifyPositive = false) {return str::toWideString(str::toString(input, specifyPositive));}
    /**
     * Get an arithmetic data type as a wide string
     * This particular overload allows for extra formatting with leading/trailing zeros
     * 
     * @param beforeDecimal The amount of digits to place before the decimal point
     * @param afterDecimal The amount of digits to place after the decimal point (even for integral types)
     * @param add If true, beforeDecimal/afterDecimal switch from amount of digits to amount of leading/trailing zeros (regardless of digits present already)
     * @returns The inputted data type formatted as a wide string
     */
    template <typename Type> const std::wstring toWideString_Places(Type input, unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) {return str::toWideString(str::toString_Places(input, beforeDecimal, afterDecimal, add, specifyPositive));}
    /**
     * Get an arithmetic data type as a wide string
     * This particular overload allows for extra formatting by specifying a length for each component
     * 
     * @param length The desired width of the output wide string
     * @param leading Whether to add the extra zeros to the beginning or end of the wide string
     * @returns The inputted data type formatted as a wide string
     */
    template <typename Type> const std::wstring toWideString_Length(Type input, unsigned long length, bool leading = true, bool specifyPositive = false) {return str::toWideString(str::toString_Length(input, length, leading, specifyPositive));}
}

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
            
            template <typename IntType> HexCoord<IntType> round(HexCoord<Type> &coord) {
                if (std::is_integral<Type>::value || std::is_floating_point<IntType>::value) {return coord;}

                IntType i = round(coord.I);
                IntType j = round(coord.J);
                IntType k = round(coord.K);

                double iDiff = fabs(i - coord.I);
                double jDiff = fabs(j - coord.J);
                double kDiff = fabs(k - coord.K);

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
            const Type getI() {return I;}
            /**
             * Get the Hecoord's j-component
             * 
             * @returns The j-component
             */
            const Type getJ() {return J;}
            /**
             * Get the Hecoord's k-component
             * 
             * @returns The k-component
             */
            const Type getK() {return K;}

            const Type getY(Type radius = 1) {
                double y = sqrt(3) * (I / 2 + J) * radius;
                return std::is_integral<Type>::value ? (y) : y;
            }
            const Type getX(Type radius = 1) {
                double x = I * (3 / 2) * radius;
                return std::is_integral<Type>::value ? round(x) : x;
            }

            const int getY_npp(int dimy) {return J * dimy + (I < 0 ? -1 : 1) * (I % 2 != 0 ? dimy / 2 : 0) + I / 2 * dimy;}
            const int getX_npp(int dimx, int inset) {return I * dimx - I * inset;}

            /**
             * Get each component arranged in a vector
             * 
             * @returns A vector containing the three components of the HexCoord {i, j, k}
             */
            const std::vector<Type> toVector() {return {I, J, K};}
            /**
             * Get the HexCoord as a string formatted as (i, j, k)
             * 
             * @returns The HexCoord formatted as a string
             */
            const std::string toString(bool specifyPositive = false) {return "(" + str::toString(I, specifyPositive) + ", " + str::toString(J, specifyPositive) + ", " + str::toString(K, specifyPositive) + ")";}
            /**
             * Get the HexCoord as a string formatted as (i, j, k)
             * This particular overload allows for extra formatting with leading/trailing zeros
             * 
             * @param beforeDecimal The amount of digits to place before the decimal point
             * @param afterDecimal The amount of digits to place after the decimal point (even for integral types)
             * @param add If true, beforeDecimal/afterDecimal switch from amount of digits to amount of leading/trailing zeros (regardless of digits present already)
             * @returns The HexCoord formatted as a string
             */
            const std::string toString_Places(unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) {return "(" + str::toString_Places(I, beforeDecimal, afterDecimal, add, specifyPositive) + ", " + str::toString_Places(J, beforeDecimal, afterDecimal, add, specifyPositive) + ", " + str::toString_Places(K, beforeDecimal, afterDecimal, add, specifyPositive) + ")";}
            /**
             * Get the HexCoord as a string formatted as (i, j, k)
             * This particular overload allows for extra formatting by specifying a length for each component
             * 
             * @param length The desired width of the output string
             * @param leading Whether to add the extra zeros to the beginning or end of the string
             * @returns The HexCoord formatted as a string
             */
            const std::string toString_Length(unsigned long length, bool leading = true, bool specifyPositive = false) {return "(" + str::toString_Length(I, length, leading, specifyPositive) + ", " + str::toString_Length(J, length, leading, specifyPositive) + ", " + str::toString_Length(K, length, leading, specifyPositive) + ")";}
            /**
             * Get the HexCoord as a wide string formatted as (i, j, k)
             * 
             * @returns The HexCoord formatted as a wide string
             */
            const std::wstring toWideString(bool specifyPositive = false) {return str::toWideString(toString(specifyPositive));}
            /**
             * Get the HexCoord as a wide string formatted as (i, j, k)
             * This particular overload allows for extra formatting with leading/trailing zeros
             * 
             * @param beforeDecimal The amount of digits to place before the decimal point
             * @param afterDecimal The amount of digits to place after the decimal point (even for integral types)
             * @param add If true, beforeDecimal/afterDecimal switch from amount of digits to amount of leading/trailing zeros (regardless of digits present already)
             * @returns The HexCoord formatted as a wide string
             */
            const std::wstring toWideString_Places(unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) {return str::toWideString(toString_Places(beforeDecimal, afterDecimal, add, specifyPositive));}
            /**
             * Get the HexCoord as a wide string formatted as (i, j, k)
             * This particular overload allows for extra formatting by specifying a length for each component
             * 
             * @param length The desired width of the output wide string
             * @param leading Whether to add the extra zeros to the beginning or end of the wide string
             * @returns The HexCoord formatted as a wide string
             */
            const std::wstring toWideString_Length(unsigned long length, bool leading = true, bool specifyPositive = false) {return str::toWideString(toString_Length(length, leading, specifyPositive));}

            /**
             * Get the HexCoord's position if rotated 180 degrees with respect to the origin
             * 
             * @returns The HexCoord's position if rotated 180 degrees with respect to the origin
             */
            const HexCoord<Type> operator ! () {return HexCoord<Type>(-I, -J, -K);}
            /**
             * Checks for equality between two HexCoords
             * 
             * Depending on what metric the HexCoord is using to measure relations, the output will differ
             * In this case, HEX_RELATE_COMMON simply combines HEX_RELATE_ICOORD, HEX_RELATE_JCOORD, and HEX_RELATE_KCOORD
             * 
             * @param coord HexCoord to compare with
             * @returns true for equality or false otherwise
             */
            const bool operator == (const HexCoord<Type> &coord) {
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
            const bool operator != (const HexCoord<Type> &coord) {return !(this == coord);}
            /**
             * Checks for whether a HexCoord is less than another HexCoord
             * 
             * Depending on what metric the HexCoord is using to measure relations, the output will differ
             * In this case, HEX_RELATE_COMMON evaluates the same as HEX_RELATE_TAXICB
             * 
             * @param coord HexCoord to compare with
             * @returns true if the HexCoord is less than the other or false otherwise
             */
            const bool operator < (const HexCoord<Type> &coord) {
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
            const bool operator <= (const HexCoord<Type> &coord) {
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
            const bool operator > (const HexCoord<Type> &coord) {
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
            const bool operator >= (const HexCoord<Type> &coord) {
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
            const unsigned char getRelationMetric() {return RelationMetric;}

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
            HexCoord<Type> operator + (const HexCoord<Type> &coord) {return HexCoord<Type>(I + coord.I, J + coord.J, K + coord.K);}
            /**
             * Subtract the position of two HexCoords
             * This is essentially vector subtraction
             * 
             * @param coord HexCoordinate to subtract from the current one
             * @returns The position of the resulting HexCoord
             */
            HexCoord<Type> operator - (const HexCoord<Type> &coord) {return HexCoord<Type>(I - coord.I, J - coord.J, K - coord.K);}
            /**
             * Multiply the position of a HexCoord by a scalar value
             * This is essentially vector multiplication
             * 
             * @param scalar Scalar value of an arithmetic type matching the HexCoord to multiply the current HexCoord by
             * @returns The position of the resulting HexCoord
             */
            HexCoord<Type> operator * (const Type &scalar) {return HexCoord<Type>(I * scalar, J * scalar, K * scalar);}

            /**
             * Get the euclidean distance between two HexCoords
             * By default, this function gets the distance from (0, 0, 0)
             * 
             * @param coord HexCoord to find distance between
             * @returns An arithmetic type (from the current HexCoord) with the euclidean distance between the specified coordinates
             */
            const Type euclideanDistance(const HexCoord<Type> &coord = HexCoord<Type>(0, 0, 0)) {
                double distance = sqrt(pow(I - coord.I, 2) + pow(J - coord.J, 2) + pow(K - coord.K, 2));
                return (std::is_integral<Type>::value) ? round(distance) : distance;
            }
            /**
             * Get the taxicab distance between two HexCoords
             * By default, this function gets the distance from (0, 0, 0)
             * 
             * @param coord HexCoord to find distance between
             * @returns An arithmetic type (from the current HexCoord) with the taxicab distance between the specified coordinates
             */
            const Type taxicabDistance(const HexCoord<Type> &coord = HexCoord<Type>(0, 0, 0)) {return (fabs(I - coord.I) + fabs(J - coord.J) + fabs(K - coord.K)) / 2;}

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
            HexCoord<Type> getAdjacent(unsigned char direction, Type offset = 1) {
                direction = abs(direction) % 6 * 2 + 1;
                return HexCoord<Type>(I + UnitDirections[direction][0] * fabs(offset), J + UnitDirections[direction][1] * fabs(offset), K + UnitDirections[direction][2] * fabs(offset));
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
            HexCoord<Type> getDiagonal(unsigned char direction, Type offset = 1) {
                direction = abs(direction) % 6 * 2;
                return HexCoord<Type>(I + UnitDirections[direction][0] * fabs(offset), J + UnitDirections[direction][1] * fabs(offset), K + UnitDirections[direction][2] * fabs(offset));
            }
            /**
             * Get the HexCoord located in the direction specified
             * This function combines getAdjacent() and getDiagonal(), so any of the twelve directions work
             * 
             * @param direction The direction to project in
             * @param offset A scalar used to get the location of further coordinates (will have a somewhat dramatic effect when direction is even due to how diagonals work with hexagons)
             * @returns The HexCoord found at the offset in the specified direction
             */
            HexCoord<Type> getNeighbor(unsigned char direction, Type offset = 1) {
                direction = abs(direction) % 12;
                return direction % 2 == 0 ? getDiagonal(direction, offset) : getAdjacent(direction, offset);
            }

            HexCoord<Type> move(unsigned char direction, Type distance = 1) {
                HexCoord<Type> output = HexCoord<Type>(I, J, K);

                direction = abs(direction) % 12 * 2;
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
                for (int i = 0; i < abs(steps); i++) {
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
    };
}

#endif /* EVERYTHING */
