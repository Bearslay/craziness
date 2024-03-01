#ifndef PHYSICS
#define PHYSICS

#include <utility>
#include <cmath>

struct {
    double g = 9.80665;
} Constants;

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
