#ifndef EVERYTHING
#define EVERYTHING

#include <unordered_map>
#include <vector>

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

template <typename Type> class HexCoord {
    private:
        Type I = 0, J = 0, K = 0;

        /**
         * Contains the twelve cardinal/intermediate directions for hexagons (seperated by 30 degrees)
         * By default, the directions are arranged for flat hexagons (vertex at 0), but is adaptable for pointy ones with an index offset of 1 (30 degrees)
         */
        const std::vector<HexCoord<Type>> UnitDirections = {
            HexCoord<Type>(2, -1, -1),
            HexCoord<Type>(1, -1, 0),
            HexCoord<Type>(1, -2, 1),
            HexCoord<Type>(0, -1, 1),
            HexCoord<Type>(-1, -1, 2),
            HexCoord<Type>(-1, 0, 1),
            HexCoord<Type>(-2, 1, 1),
            HexCoord<Type>(-1, 1, 0),
            HexCoord<Type>(-1, 2, -1),
            HexCoord<Type>(0, 1, -1),
            HexCoord<Type>(1, 1, -2),
            HexCoord<Type>(1, 0, -1)
        };
        unsigned char LastDirectionMoved = HEX_DIR_000;

    public:
        /**
         * Parametric constructor
         * 
         * @param i i-component of the coordinate
         * @param j j-component of the coordinate
         * @param k k-component of the coordinate
         */
        HexCoord(Type i, Type j, Type k) {
            static_assert(std::is_arithmetic<Type>::value, "Type must be an arithmetic type");
            I = i;
            J = j;
            K = k;
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
        /**
         * Deconstructor
         */
        ~HexCoord() {}

        /**
         * Checks for equality between two HexCoordinates
         * 
         * @param coord HexCoord to compare with
         * @returns true for equal positions or false for unequal positions
        */
        const bool operator == (const HexCoord<Type> &coord) {return I == coord.I && J == coord.J && K == coord.K;}
        /**
         * Checks for inequality between two HexCoordinates
         * 
         * @param coord HexCoord to compare with
         * @returns true for unequal positions or false for equal positions
        */
        const bool operator != (const HexCoord<Type> &coord) {return !(this == coord);}

        HexCoord<Type> operator += (const HexCoord<Type> &coord) {
            HexCoord<Type> output = this;
            I += coord.I;
            J += coord.J;
            K += coord.K;
            return output;
        }
        HexCoord<Type> operator -= (const HexCoord<Type> &coord) {
            HexCoord<Type> output = this;
            I -= coord.I;
            J -= coord.J;
            K -= coord.K;
            return output;
        }
        HexCoord<Type> operator *= (const Type &scalar) {
            HexCoord<Type> output = this;
            I *= scalar;
            J *= scalar;
            K *= scalar;
            return output;
        }
        HexCoord<Type> operator + (const HexCoord<Type> &coord) {return HexCoord<Type>(I + coord.I, J + coord.J, K + coord.K);}
        HexCoord<Type> operator - (const HexCoord<Type> &coord) {return HexCoord<Type>(I - coord.I, J - coord.J, K - coord.K);}
        HexCoord<Type> operator * (const Type &scalar) {return HexCoord<Type>(I * scalar, J * scalar, K * scalar);}
};

template <typename Type> class HexPlane {
    private:
        std::unordered_map<HexCoord, Type> List;

    public:

};

#endif /* EVERYTHING */
