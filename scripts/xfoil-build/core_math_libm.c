/* Linked into every build of the reference (scripts/xfoil-build.sh).
 *
 * XFOIL's transcendentals (EXP, LOG, LOG10, real-exponent **, TANH, SIN, COS, ATAN, ATAN2,
 * ASIN) compile to calls of the C library's exp, log, log10, pow, tanh, sin, cos, atan, atan2
 * and asin. Where GCC fuses a SIN and a COS of the same argument it calls sincos (glibc targets) or,
 * on a target without one (Darwin), cexp of the imaginary argument: both are here too, cexp built
 * from cr_exp, cr_cos and cr_sin, which for the zero real part GCC passes gives exactly
 * (cr_cos(x), cr_sin(x)), as separate SIN and COS calls would. Defined here, in the executable, they take precedence over the host libm's and
 * forward to CORE-MATH's correctly rounded functions (third-party/core-math), the same sources
 * yFoil's build.rs compiles. A correctly rounded result is defined by its argument alone, so the
 * reference and yFoil get the same bits on every host, where the host libms (Apple libSystem,
 * glibc) differ by one ULP on a fraction of inputs (CLAUDE.md Rule 1).
 *
 * MDES's other complex functions (cabs, clog) stay with the host: inverse design is out of scope.
 * xfoil-build.sh fails a build in which any other maths function is left to the host libm.
 */
#include <complex.h>

double cr_exp(double);
double cr_log(double);
double cr_log10(double);
double cr_pow(double, double);
double cr_tanh(double);
double cr_sin(double);
double cr_cos(double);
double cr_atan(double);
double cr_atan2(double, double);
double cr_asin(double);
void cr_sincos(double, double *, double *);

double exp(double x) { return cr_exp(x); }
double log(double x) { return cr_log(x); }
double log10(double x) { return cr_log10(x); }
double pow(double x, double y) { return cr_pow(x, y); }
double tanh(double x) { return cr_tanh(x); }
double sin(double x) { return cr_sin(x); }
double cos(double x) { return cr_cos(x); }
double atan(double x) { return cr_atan(x); }
double atan2(double y, double x) { return cr_atan2(y, x); }
double asin(double x) { return cr_asin(x); }
void sincos(double x, double *s, double *c) { cr_sincos(x, s, c); }
double complex cexp(double complex z) {
    double e = cr_exp(creal(z));
    return CMPLX(e * cr_cos(cimag(z)), e * cr_sin(cimag(z)));
}
