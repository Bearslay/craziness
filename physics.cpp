#include "physics.hpp"
#include <iostream>

int main() {
    double y0 = 10;
    double yf = 20;
    double x0 = 0;
    double xf = 100;

    double dx = xf - x0;
    double dy = yf - y0;
    double mag = 34;

    double a = kmtcs::launchAngle_Deg(dx, dy, mag);
    if (isnanf(a)) {
        std::cout << "Can't reach target with given launch velocity\n";
        return -1;
    }

    std::cout << "Change in Height (dy):    " << dy << " m\n";
    std::cout << "Initial Height (y0):      " << y0 << " m\n";
    std::cout << "Initial X-Position (x0):  " << x0 << " m\n\n";

    std::cout << "Launch Velocity (V0):     " << mag << " m/s\n";
    std::cout << "Launch Velocity X (V0x):  " << breakVector_Deg(mag, a).first << " m/s\n";
    std::cout << "Launch Velocity Y (V0y):  " << breakVector_Deg(mag, a).second << " m/s\n";
    std::cout << "Launch Angle:             " << a << " deg\n\n";

    std::cout << "Change in Position (dx):  " << dx << " m\n";
    std::cout << "Final Height (yf):        " << yf << " m\n";
    std::cout << "Final X-Position (xf):    " << xf << " m\n\n";

    std::cout << "Landing Velocity (Vf):    " << kmtcs::landingVelocity_Deg(mag, a, dy) << " m/s\n";
    std::cout << "Landing Velocity X (Vfx): " << breakVector_Deg(kmtcs::landingVelocity_Deg(mag, a, dy), kmtcs::landingAngle_Deg(mag, a, dy)).first << " m/s\n";
    std::cout << "Landing Velocity Y (Vfy): " << breakVector_Deg(kmtcs::landingVelocity_Deg(mag, a, dy), kmtcs::landingAngle_Deg(mag, a, dy)).second << " m/s\n";
    std::cout << "Landing Angle:            " << kmtcs::landingAngle_Deg(mag, a, dy) << " deg\n\n";

    std::cout << "Maximum Height:           " << kmtcs::maxHeight_Deg(mag, a, y0) << " m\n";
    std::cout << "Time to Reach Max Height: " << kmtcs::peakTime_Deg(mag, a) << " sec\n";
    std::cout << "Total Airtime:            " << kmtcs::airTime_Deg(mag, a, dy) << " sec\n\n";

    Vector_3D<double> vr = Vector_3D<double>(5, 4, 3);
    std::cout << vr.rToString_Length(5) << "\n";
    std::cout << vr.sToString_Length(5) << "\n";
    Vector_3D<double> vs = Vector_3D<double>(10, 30, 45, true);
    std::cout << vs.sToString_Length(5) << "\n";
    std::cout << vs.rToString_Length(5) << "\n\n";

    Coord_3D<int> c1 = Coord_3D<int>(4, 7, 2);
    std::cout << c1.toString_Places(2, 2) << "\n\n";

    double num = 2124.24532;
    double num1 = -0.0425334;
    long double avagadro = 602000000000000000000000.0;
    int num2 = 12;
    std::cout << astr::toString_Places(num, 4, 4) << "\n";
    std::cout << astr::toString_Places(astr::round(num, 3), 4, 4) << "\n\n";
    std::cout << astr::toString(avagadro) << "\n";
    std::cout << astr::toString_Sci(avagadro) << "\n";
    std::cout << astr::toString(num1) << "\n";
    std::cout << astr::toString_Sci(num1) << "\n";
    std::cout << astr::toString_Sci(num1, 2, 2, false) << "\n";
    std::cout << astr::toString(num2) << "\n";
    std::cout << astr::toString_Sci(num2) << "\n";

    return 0;
}
