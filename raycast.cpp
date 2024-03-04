#include "ncursespp.hpp"
#include "physics.hpp"
#include <iostream>
#include <unordered_map>

template <typename ArithType> Coord_3D<ArithType> projectCoord(const Coord_3D<ArithType> &coord, double focalLength) {
    static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");

    // Protect from a floating-point exception, even though a projection doesn't occur
    if (focalLength + coord.getX() == 0) {return coord;}

    Coord_3D<ArithType> output;

    output.setY((focalLength * coord.getY()) / (focalLength + coord.getX()));
    output.setZ((focalLength * coord.getZ()) / (focalLength + coord.getX()));

    return output;
}

template <typename ArithType> std::vector<Coord_3D<int>> extractLine(Coord_3D<ArithType> &point1, Coord_3D<ArithType> &point2) {
    static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");

    // Extract the x and y values out of the input coordinates for ease
    int x1 = std::round(point1.getY()), x2 = std::round(point2.getY()), y1 = std::round(point1.getZ()), y2 = std::round(point2.getZ());

    // Calculate the differences
    int dx = x2 - x1;
    int dy = y2 - y1;

    // Swap x and y if the slope is greater than 1
    bool steep = std::abs(dy) > std::abs(dx);
    if (steep) {
        int temp = x1;
        x1 = y1;
        y1 = temp;
        temp = x2;
        x2 = y2;
        y2 = temp;
    }

    // Ensure that the second coordinate is to the right of the first one
    bool swapLR = false;
    if (x1 > x2) {
        int temp = x1;
        x1 = x2;
        x2 = temp;
        temp = y1;
        y1 = y2;
        y2 = temp;
        swapLR = true;
    }

    // Recalculate
    dx = x2 - x1;
    dy = y2 - y1;

    // Reflect the line across the x-axis if the slope is less than 0
    bool flipX = y1 > y2;
    char ystep = flipX ? -1 : 1;

    int y = y1;
    int error = (int)(dx / 2);
    std::vector<Coord_3D<int>> output;
    for (int x = x1; x <= x2; x++) {
        if (steep) {output.emplace_back(Coord_3D<int>(0, y, x));}
        else {output.emplace_back(Coord_3D<int>(0, x, y));}

        error -= std::abs(dy);
        if (error < 0) {
            y += ystep;
            error += dx;
        }
    }

    // Reverse the order of the output vector if the original x1 > x2
    if (swapLR) {
        for (int i = 0; i < output.size() / 2; i++) {
            Coord_3D<int> temp = output[i];
            output[i] = output[output.size() - 1 - i];
            output[output.size() - 1 - i] = temp;
        }
    }

    return output;
}

template <typename ArithType> Coord_3D<ArithType> rotateCoord_3D(const Vector_3D<ArithType> &axis, const Coord_3D<ArithType> &coord, double theta, bool degrees = true) {
    static_assert(std::is_arithmetic<ArithType>::value, "ArithType must be an arithmetic type");

    // Convert the angle to radians if needed
    if (degrees) {theta *= M_PI / 180;}

    ArithType x = axis.getX();
    ArithType y = axis.getY();
    ArithType z = axis.getZ();

    const std::vector<std::vector<double>> R = {
        {std::cos(theta) + std::pow(x, 2) * (1 - std::cos(theta)), x * y * (1 - std::cos(theta)) - z * std::sin(theta), x * z * (1 - std::cos(theta)) + y * std::sin(theta)},
        {y * x * (1 - std::cos(theta)) + z * std::sin(theta), std::cos(theta) + std::pow(y, 2) * (1 - std::cos(theta)), y * z * (1 - std::cos(theta)) - x * std::sin(theta)},
        {z * x * (1 - std::cos(theta)) - y * std::sin(theta), z * y * (1 - std::cos(theta)) + x * std::sin(theta), std::cos(theta) + std::pow(z, 2) * (1 - std::cos(theta))}
    };

    std::vector<ArithType> original = coord.toVector();
    std::vector<double> rotated = {0, 0, 0};
    for (unsigned char i = 0; i < 3; i++) {
        for (unsigned char j = 0; j < 3; j++) {
            rotated[i] += original[j] * R[i][j];
        }
    }

    if (std::is_floating_point<ArithType>::value) {
        return Coord_3D<ArithType>(rotated[0], rotated[1], rotated[2]);
    }
    
    return Coord_3D<ArithType>(std::round(rotated[0]), std::round(rotated[1]), std::round(rotated[2]));
}

int main() {
    npp::init();
   
    std::vector<std::pair<unsigned char, unsigned char>> pairs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };
    std::vector<Coord_3D<int>> points = {
        Coord_3D<int>(-10, -10,  10),
        Coord_3D<int>(-10,  10,  10),
        Coord_3D<int>( 10,  10,  10),
        Coord_3D<int>( 10, -10,  10),
        Coord_3D<int>(-10, -10, -10),
        Coord_3D<int>(-10,  10, -10),
        Coord_3D<int>( 10,  10, -10),
        Coord_3D<int>( 10, -10, -10)
    };
    std::vector<Coord_3D<int>> incremental = {
        Coord_3D<int>(-1, -1,  1),
        Coord_3D<int>(-1,  1,  1),
        Coord_3D<int>( 1,  1,  1),
        Coord_3D<int>( 1, -1,  1),
        Coord_3D<int>(-1, -1, -1),
        Coord_3D<int>(-1,  1, -1),
        Coord_3D<int>( 1,  1, -1),
        Coord_3D<int>( 1, -1, -1)    
    };

    std::vector<Coord_3D<int>> pointsNow, projected;
    for (unsigned char i = 0; i < points.size(); i++) {
        pointsNow.emplace_back(points[i]);
        projected.emplace_back();
    }

    std::vector<double> angles;
    for (unsigned char i = 0; i <= 24; i++) {
        angles.emplace_back(i * 15);
    }

    double focalLength = 50;
    double theta = -45;
    double phi = 90;

    npp::mwin.dhline(npp::mwin.gdimy() / 2, 0, npp::mwin.gdimx(), false, {DOUBLED_HORIZONTAL, DASHED_NONE});
    npp::mwin.dvline(0, npp::mwin.gdimx() / 2, npp::mwin.gdimy());
    npp::mwin.dvline(0, npp::mwin.gdimx() / 2 + 1, npp::mwin.gdimy());
    npp::mwin.dbox();

    int ch;
    char ax, ay, az;
    bool goodCycle;
    while (true) {
        goodCycle = true;
        npp::mwin.reset();
        npp::mwin.dhline(npp::mwin.gdimy() / 2, 0, npp::mwin.gdimx(), false, {DOUBLED_HORIZONTAL, DASHED_NONE});
        npp::mwin.dvline(0, npp::mwin.gdimx() / 2, npp::mwin.gdimy());
        npp::mwin.dvline(0, npp::mwin.gdimx() / 2 + 1, npp::mwin.gdimy());
        npp::mwin.dbox();

        switch (ch) {
            case 'i':
                ay = (ay == -23) ? 0 : ay - 1;
                break;
            case 'k':
                ay = (ay == 23) ? 0 : ay + 1;
                break;
            case 'j':
                az = (az == 23) ? 0 : az + 1;
                break;
            case 'l':
                az = (az == -23) ? 0 : az - 1;
                break;
            case 'u':
                ax = (ax == -23) ? 0 : ax - 1;
                break;
            case 'o':
                ax = (ax == -23) ? 0 : ax + 1;
                break;
            case 'p':
                focalLength -= 5;
                break;
            case ';':
                focalLength += 5;
                break;
            case 'n':
                for (unsigned char i = 0; i < points.size(); i++) {
                    points[i] -= incremental[i];
                }
                break;
            case 'm':
                for (unsigned char i = 0; i < points.size(); i++) {
                    points[i] += incremental[i];
                }
                break;
        }

        for (unsigned char i = 0; i < points.size(); i++) {
            pointsNow[i] = rotateCoord_3D(Vector_3D<int>(0, 1, 0), points[i], 15 * ay);
            pointsNow[i] = rotateCoord_3D(Vector_3D<int>(0, 0, 1), pointsNow[i], 15 * az);
            pointsNow[i] = rotateCoord_3D(Vector_3D<int>(1, 0, 0), pointsNow[i], 15 * ax);
        }

        bool revertChange = false;
        for (unsigned char i = 0; i < points.size(); i++) {
            
        }

        for (unsigned char i = 0; i < points.size() i++) {
            projected[i] = projectCoord<int>(pointsNow[i], focalLength);
        }

        for (unsigned char i = 0; i < pairs.size(); i++) {
            std::vector<Coord_3D<int>> line = extractLine(projected[pairs[i].first], projected[pairs[i].second]);

            for (unsigned char j = 0; j < line.size(); j++) {
                npp::mwin.wstr(npp::mwin.gdimy() / 2 + line[j].getZ(), npp::mwin.gdimx() / 2 + line[j].getY() * 2, L"▒▒");
            }
        }
        for (unsigned char i = 0; i < projected.size(); i++) {
            npp::mwin.wstr(npp::mwin.gdimy() / 2 + projected[i].getZ(), npp::mwin.gdimx() / 2 + projected[i].getY() * 2, L"██", 5 + i);
        }

        if ((ch = npp::mwin.gchar()) == 'q') {break;}
    }

    return npp::end();
}
