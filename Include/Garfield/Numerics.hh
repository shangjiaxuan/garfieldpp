#ifndef G_NUMERICS_H
#define G_NUMERICS_H

#include <functional>
#include <complex>
#include <vector>

namespace Garfield {

/// Collection of numerical routines.
namespace Numerics {

/// Functions for performing numerical integration (quadrature). 
/// Reimplemented from %QUADPACK.
///  
/// R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner,
/// %QUADPACK, a Subroutine Package for Automatic Integration,
/// Springer, 1983
namespace QUADPACK {

/// Estimates an integral over a semi-infinite or infinite interval.
/// Calculates an approximation to an integral
/// \f[ I = \int_{a}^{\infty} f\left(x\right) dx, \f] 
/// or 
/// \f[ I = \int_{-\infty}^{a} f\left(x\right) dx, \f] 
/// or 
/// \f[ I = \int_{-\infty}^{\infty} f\left(x\right) dx \f] 
/// hopefully satisfying
/// \f[
///     \left| I - \mathrm{RESULT} \right| \leq
///     \max(\varepsilon_{\mathrm{abs}}, \varepsilon_{\mathrm{rel}} \left|I\right|).
/// \f]
/// \param f function to be integrated.
/// \param bound value of the finite endpoint of the integration range (if any).
/// \param inf indicates the type of integration range.
///   -  1: (    bound, +Infinity)
///   - -1: (-Infinity,     bound)
///   -  2: (-Infinity, +Infinity)
/// \param epsabs requested absolute accuracy
/// \param epsrel requested relative accuracy
/// \param result the estimated value of the integral 
/// \param abserr estimated error
/// \param status error flag
///   -   0: normal and reliable termination, requested accuracy
///           has been achieved.
///   - > 0: abnormal termination, estimates for result and error
///           are less reliable.
///     - 1: maximum number of subdivisions reached.
///     - 2: occurance of roundoff error prevents the requested 
///          tolerance from being achieved. Error may be underestimated.
///     - 3: extremely bad integrand behaviour at some points of the 
///          integration interval.
///     - 4: algorithm does not converge, roundoff error is detected in
///          the extrapolation table. It is assumed that the requested
///          tolerance cannot be achieved and that the returned result
///          is the best that can be obtained.
///     - 5: integral is probably divergent, or slowly convergent.
///     - 6: invalid input.  
void qagi(std::function<double(double)> f, double bound, const int inf, 
          const double epsabs, const double epsrel, 
          double& result, double& abserr, unsigned int& status);

/// 15-point Gauss-Kronrod integration with (semi-)infinite integration range.
/// \param f function to be integrated.
/// \param bound finite bound of original integration range (0 if inf = 2).
/// \param inf indicates the type of integration range.
/// \param a lower limit for integration over subrange of (0, 1).
/// \param b upper limit for integration over subrange of (0, 1).
/// \param result approximation to the integral.
/// \param abserr estimate of the modulus of the absolute error.
/// \param resabs approximation to the integral over \f$\left|f\right|\f$.
/// \param resasc approximation to the integral of 
///               \f$\left|f - I / (b-a)\right|\f$ over \f$(a,b)\f$.
void qk15i(std::function<double(double)> f, double bound, const int inf,
           const double a, const double b, double& result, double& abserr,
           double& resabs, double& resasc);

/// 15-point Gauss-Kronrod integration with finite integration range.
void qk15(std::function<double(double)> f, const double a, const double b,
          double& result, double& abserr, double& resabs, double& resasc);

}

// Linear algebra routines from CERNLIB
/// Replace B by the solution X of A*X=B, after which A is undefined.
void Deqn(const int n, std::vector<std::vector<double> >& a,
          int& ifail, std::vector<double>& b);
void Dfact(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir, int& ifail, double& det, int& jfail);
void Dfeqn(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir, std::vector<double>& b);
void Dfinv(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir);
/// Replace b by the solution x of A x = b, and replace A by its inverse.
void Deqinv(const int n, std::vector<std::vector<double> >& a, int& ifail,
            std::vector<double>& b);

void Cfact(const int n, std::vector<std::vector<std::complex<double> > >& a,
           std::vector<int>& ir, int& ifail, std::complex<double>& det,
           int& jfail);
void Cfinv(const int n, std::vector<std::vector<std::complex<double> > >& a,
           std::vector<int>& ir);
/// Replace A by its inverse.
void Cinv(const int n, std::vector<std::vector<std::complex<double> > >& a,
          int& ifail);


/// Modified Bessel functions.
/// Series expansions from Abramowitz and Stegun.
inline double BesselI0S(const double xx) {
  const double y = xx / 3.75;
  const double y2 = y * y;
  return 1. + 3.5156229 * y2 + 3.0899424 * y2 * y2 +
         1.2067492 * pow(y2, 3) + 0.2659732 * pow(y2, 4) +
         0.0360768 * pow(y2, 5) + 0.0045813 * pow(y2, 6);
}

inline double BesselI1S(const double xx) {
  const double y = xx / 3.75;
  const double y2 = y * y;
  return xx *
         (0.5 + 0.87890594 * y2 + 0.51498869 * y2 * y2 + 
          0.15084934 * pow(y2, 3) + 0.02658733 * pow(y2, 4) + 
          0.00301532 * pow(y2, 5) + 0.00032411 * pow(y2, 6));
}

inline double BesselK0S(const double xx) {
  const double y = xx / 2.;
  const double y2 = y * y;
  return -log(y) * BesselI0S(xx) - 0.57721566 +
         0.42278420 * y2 + 0.23069756 * y2 * y2 +
         0.03488590 * pow(y2, 3) + 0.00262698 * pow(y2, 4) +
         0.00010750 * pow(y2, 5) + 0.00000740 * pow(y2, 6);
}

inline double BesselK0L(const double xx) {
  const double y = 2. / xx;
  return (exp(-xx) / sqrt(xx)) *
         (1.25331414 - 0.07832358 * y + 0.02189568 * y * y -
          0.01062446 * pow(y, 3) + 0.00587872 * pow(y, 4) -
          0.00251540 * pow(y, 5) + 0.00053208 * pow(y, 6));
}

inline double BesselK1S(const double xx) {
  const double y = xx / 2.;
  const double y2 = y * y;
  return log(y) * BesselI1S(xx) +
         (1. / xx) *
             (1. + 0.15443144 * y2 - 0.67278579 * y2 * y2 -
              0.18156897 * pow(y2, 3) - 0.01919402 * pow(y2, 4) -
              0.00110404 * pow(y2, 5) - 0.00004686 * pow(y2, 6));
}

inline double BesselK1L(const double xx) {
  const double y = 2. / xx;
  return (exp(-xx) / sqrt(xx)) *
         (1.25331414 + 0.23498619 * y - 0.03655620 * y * y +
          0.01504268 * pow(y, 3) - 0.00780353 * pow(y, 4) +
          0.00325614 * pow(y, 5) - 0.00068245 * pow(y, 6));
}

/// C++ version of DIVDIF (CERN program library E105) which performs
/// tabular interpolation using symmetrically placed argument points.
double Divdif(const std::vector<double>& f, const std::vector<double>& a,
              int nn, double x, int mm);

/// Interpolation of order 1 and 2 in an irregular rectangular
/// two-dimensional grid.
bool Boxin2(const std::vector<std::vector<double> >& value,
            const std::vector<double>& xAxis, const std::vector<double>& yAxis,
            const int nx, const int ny, const double xx, const double yy,
            double& f, const int iOrder);
/// Interpolation of order 1 and 2 in an irregular rectangular 
/// three-dimensional grid.
bool Boxin3(const std::vector<std::vector<std::vector<double> > >& value,
            const std::vector<double>& xAxis, const std::vector<double>& yAxis,
            const std::vector<double>& zAxis, const int nx, const int ny,
            const int nz, const double xx, const double yy, const double zz,
            double& f, const int iOrder);
}

}

#endif
