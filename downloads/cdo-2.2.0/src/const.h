#ifndef CONST_H
#define CONST_H

constexpr double MIN_PS = 20000.;
constexpr double MAX_PS = 120000.;

constexpr double MIN_FIS = -100000.;
constexpr double MAX_FIS = 100000.;

constexpr double MIN_T = 150.;
constexpr double MAX_T = 400.;

constexpr double MIN_Q = 0.0;
constexpr double MAX_Q = 0.1;

#ifndef M_LN10
#define M_LN10 2.30258509299404568401799145468436421 /* log_e 10 */
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 /* pi */
#endif

#ifndef RAD_CONVERT
#define RAD_CONVERT
constexpr double RAD2DEG = 180. / M_PI;  // conversion for rad to deg
constexpr double DEG2RAD = M_PI / 180.;  // conversion for deg to rad
#endif

#endif /* CONST_H */
