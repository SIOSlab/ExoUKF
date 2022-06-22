#ifndef ANGLES_HPP
#define ANGLES_HPP

constexpr double pi = 3.14159265358979323846;

inline double rad2deg(double x) {
    return x * (180 / pi);
}

inline double deg2rad(double x) {
    return x * (pi / 180);
}

inline double rad2min(double x) {
    return x * (180 * 60 / pi);
}

inline double min2rad(double x) {
    return x * (pi / (180 * 60));
}

inline double sec2rad(double x) {
    return x * (pi / (180 * 60 * 60));
}

inline double rad2sec(double x) {
    return x * ((180 * 60 * 60) / pi);
}

#endif
