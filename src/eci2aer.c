// Much of this code is inspired/ported from the matlab/octave coordinate
// conversions here: https://github.com/geospace-code/matmap3d

#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

void eci2aer(double jdate, double *eci, double *lla, double *aer) {

    double ecef[3] = {0, 0, 0};

    // jdate to greenwhich sidereal time
    double tUT1 = (jdate - 2451545) / 36525;
    double tUT1_power2 = (0.093104 * tUT1) * tUT1;
    double tUT1_power3 = (6.2e-6 * tUT1) * tUT1 * tUT1;
    double gmst_sec = 67310.54841 + (3164400184.812866) * tUT1 + tUT1_power2 - tUT1_power3;
    // 1/86400 and %(2*pi) implied by units of radians
    double gst = fmod(gmst_sec * (2 * PI) / 86400.0, 2 * PI);
    // printf("greenwhich stime: %f\n", gst * (180.0/PI));

    // Satellite ECI to ECEF
    ecef[0] = cos(gst)*eci[0]*1000.0 + sin(gst)*eci[1]*1000.0;
    ecef[1] = -sin(gst)*eci[0]*1000.0 + cos(gst)*eci[1]*1000.0;
    ecef[2] = eci[2]*1000.0;

    // LLA2ECEF (Observer)
    float rad_lat = lla[0] * (PI/180.0);
    float rad_lon = lla[1] * (PI/180.0);

    float a = 6378137.0;
    float finv = 298.257223563;
    float f = 1 / finv;
    float e2 = 1 - (1 - f) * (1 - f);
    float v = a / sqrt(1 - e2 * sin(rad_lat) * sin(rad_lat));
    float x0 = (v + lla[2]) * cos(rad_lat) * cos(rad_lon);
    float y0 = (v + lla[2]) * cos(rad_lat) * sin(rad_lon);
    float z0 = (v * (1 - e2) + lla[2]) * sin(rad_lat);
    // printf("X: %f, Y: %f Z:%f\n", x0, y0, z0);

    // ECEF2ENU
    float U = ecef[0] - x0;
    float V = ecef[1] - y0;
    float W = ecef[2] - z0;

    float t = cos(rad_lon) * U + sin(rad_lon) * V;
    float e = -sin(rad_lon) * U + cos(rad_lon) * V;
    float up = cos(rad_lat) * t + sin(rad_lat) * W;
    float n = -sin(rad_lat) * t + cos(rad_lat) * W;

    if (fabs(e) < 0.001) { e = 0; }
    if (fabs(n) < 0.001) { n = 0; }
    if (fabs(up) < 0.001) { up = 0; }
    // printf("E: %f, N: %f U:%f\n", e, n, up);

    float r = hypot(e, n);
    float slantrange = hypot(r, up);
    // % radians
    float elev = atan2(up, r);
    float az = fmod(atan2(e, n), 2*PI);
    // printf("Az: %f, El: %f Range:%f\n", az * (180.0/PI), elev * (180.0/PI), slantrange);
    aer[0] = az * (180.0/PI);
    aer[1] = elev * (180.0/PI);
    aer[2] = slantrange;
}