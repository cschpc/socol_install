#ifndef GRID_CONVERT_H
#define GRID_CONVERT_H

#include <cmath>

#ifdef HAVE_SINCOS

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

void sincos(double x, double *sin, double *cos);
static inline void
gcLLtoXYZ(double lon, double lat, double *xyz)
{
  double sinlon, coslon, sinlat, coslat;
  sincos(lon, &sinlon, &coslon);
  sincos(lat, &sinlat, &coslat);
  xyz[0] = coslat * coslon;
  xyz[1] = coslat * sinlon;
  xyz[2] = sinlat;
}

#else

static inline void
gcLLtoXYZ(double lon, double lat, double *xyz)
{
  const auto cos_lat = std::cos(lat);
  xyz[0] = cos_lat * std::cos(lon);
  xyz[1] = cos_lat * std::sin(lon);
  xyz[2] = std::sin(lat);
}

#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288  // pi
#endif

#ifndef RAD_CONVERT
#define RAD_CONVERT
constexpr double RAD2DEG = 180.0 / M_PI;  // conversion for rad to deg
constexpr double DEG2RAD = M_PI / 180.0;  // conversion for deg to rad
#endif

static inline void
gcLLtoXYZdeg(double lon, double lat, double *xyz)
{
  gcLLtoXYZ(lon * DEG2RAD, lat * DEG2RAD, xyz);
}

#endif
