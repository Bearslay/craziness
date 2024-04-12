#include "physics.hpp"
#include <iostream>

int main() {
    bool useRadians = false;
    double y0 = 10;
    double yf = 20;
    double x0 = 0;
    double xf = 100;

    double dx = xf - x0;
    double dy = yf - y0;
    double mag = 34;

    double a = kmtcs::launchAngle<double>(dx, dy, mag, true, useRadians);
    std::pair<double, double> landing = kmtcs::landingVector<double>(mag, a, dy, useRadians);
    if (isnanf(a)) {
        std::cout << "Can't reach target with given launch velocity\n";
        return -1;
    }

    std::cout << " Change in Position (dx): " << dx << " m\n";
    std::cout << "       Final Height (yf): " << yf << " m\n";
    std::cout << "   Final X-Position (xf): " << xf << " m\n\n";

    std::cout << "   Change in Height (dy): " << dy << " m\n";
    std::cout << "     Initial Height (y0): " << y0 << " m\n";
    std::cout << " Initial X-Position (x0): " << x0 << " m\n\n";

    std::cout << "    Launch Velocity (V0): " << mag << " m/s\n";
    std::cout << " Launch Velocity X (V0x): " << breakVector<double>(mag, a, useRadians).first << " m/s\n";
    std::cout << " Launch Velocity Y (V0y): " << breakVector<double>(mag, a, useRadians).second << " m/s\n";
    std::cout << "            Launch Angle: " << a << " deg\n\n";

    std::cout << "   Landing Velocity (Vf): " << landing.first << " m/s\n";
    std::cout << "Landing Velocity X (Vfx): " << breakVector<double>(landing.first, landing.second, useRadians).first << " m/s\n";
    std::cout << "Landing Velocity Y (Vfy): " << breakVector<double>(landing.first, landing.second, useRadians).second << " m/s\n";
    std::cout << "           Landing Angle: " << landing.second << " deg\n\n";

    std::cout << "          Maximum Height: " << kmtcs::maxHeight<double>(mag, a, y0, useRadians) << " m\n";
    std::cout << "Time to Reach Max Height: " << kmtcs:: peakTime<double>(mag, a,     useRadians) << " sec\n";
    std::cout << "           Total Airtime: " << kmtcs::  airTime<double>(mag, a, dy, useRadians) << " sec\n\n";

    return 0;
}
