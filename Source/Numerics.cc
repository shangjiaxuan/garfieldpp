#include <cmath>
#include <iostream>
#include <algorithm>
#include <array>

#include "Garfield/Numerics.hh"

namespace Garfield {

namespace Numerics {

void Dfact(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir, int& ifail, double& det, int& jfail) {
  constexpr double g1 = 1.e-19;
  constexpr double g2 = 1.e-19;

  double t;
  int k;

  ifail = jfail = 0;

  int nxch = 0;
  det = 1.;

  for (int j = 1; j <= n; ++j) {
    k = j;
    double p = fabs(a[j - 1][j - 1]);
    if (j == n) {
      if (p <= 0.) {
        det = 0.;
        ifail = -1;
        jfail = 0;
        return;
      }
      det *= a[j - 1][j - 1];
      a[j - 1][j - 1] = 1. / a[j - 1][j - 1];
      t = fabs(det);
      if (t < g1) {
        det = 0.;
        if (jfail == 0) jfail = -1;
      } else if (t > g2) {
        det = 1.;
        if (jfail == 0) jfail = +1;
      }
      continue;
    }
    for (int i = j + 1; i <= n; ++i) {
      double q = fabs(a[i - 1][j - 1]);
      if (q <= p) continue;
      k = i;
      p = q;
    }
    if (k != j) {
      for (int l = 1; l <= n; ++l) {
        double tf = a[j - 1][l - 1];
        a[j - 1][l - 1] = a[k - 1][l - 1];
        a[k - 1][l - 1] = tf;
      }
      ++nxch;
      ir[nxch - 1] = j * 4096 + k;
    } else if (p <= 0.) {
      det = 0.;
      ifail = -1;
      jfail = 0;
      return;
    }
    det *= a[j - 1][j - 1];
    a[j - 1][j - 1] = 1. / a[j - 1][j - 1];
    t = fabs(det);
    if (t < g1) {
      det = 0.;
      if (jfail == 0) jfail = -1;
    } else if (t > g2) {
      det = 1.;
      if (jfail == 0) jfail = +1;
    }
    for (k = j + 1; k <= n; ++k) {
      double s11 = -a[j - 1][k - 1];
      double s12 = -a[k - 1][j];
      if (j == 1) {
        a[j - 1][k - 1] = -s11 * a[j - 1][j - 1];
        a[k - 1][j] = -(s12 + a[j - 1][j] * a[k - 1][j - 1]);
        continue;
      }
      for (int i = 1; i <= j - 1; ++i) {
        s11 += a[i - 1][k - 1] * a[j - 1][i - 1];
        s12 += a[i - 1][j] * a[k - 1][i - 1];
      }
      a[j - 1][k - 1] = -s11 * a[j - 1][j - 1];
      a[k - 1][j] = -a[j - 1][j] * a[k - 1][j - 1] - s12;
    }
  }

  if (nxch % 2 != 0) det = -det;
  if (jfail != 0) det = 0.;
  ir[n - 1] = nxch;
}

void Dfeqn(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir, std::vector<double>& b) {
  if (n <= 0) return;

  int nxch = ir[n - 1];
  if (nxch != 0) {
    for (int m = 1; m <= nxch; ++m) {
      int ij = ir[m - 1];
      int i = ij / 4096;
      int j = ij % 4096;
      double te = b[i - 1];
      b[i - 1] = b[j - 1];
      b[j - 1] = te;
    }
  }

  b[0] *= a[0][0];
  if (n == 1) return;

  for (int i = 2; i <= n; ++i) {
    double s21 = -b[i - 1];
    for (int j = 1; j <= i - 1; ++j) {
      s21 += a[i - 1][j - 1] * b[j - 1];
    }
    b[i - 1] = -a[i - 1][i - 1] * s21;
  }

  for (int i = 1; i <= n - 1; ++i) {
    double s22 = -b[n - i - 1];
    for (int j = 1; j <= i; ++j) {
      s22 += a[n - i - 1][n - j] * b[n - j];
    }
    b[n - i - 1] = -s22;
  }
}

void Dfinv(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir) {

  if (n <= 1) return;
  a[1][0] = -a[1][1] * a[0][0] * a[1][0];
  a[0][1] = -a[0][1];
  if (n > 2) {
    for (int i = 3; i <= n; ++i) {
      for (int j = 1; j <= i - 2; ++j) {
        double s31 = 0.;
        double s32 = a[j - 1][i - 1];
        for (int k = j; k <= i - 2; ++k) {
          s31 += a[k - 1][j - 1] * a[i - 1][k - 1];
          s32 += a[j - 1][k] * a[k][i - 1];
        }
        a[i - 1][j - 1] =
            -a[i - 1][i - 1] * (s31 + a[i - 2][j - 1] * a[i - 1][i - 2]);
        a[j - 1][i - 1] = -s32;
      }
      a[i - 1][i - 2] = -a[i - 1][i - 1] * a[i - 2][i - 2] * a[i - 1][i - 2];
      a[i - 2][i - 1] = -a[i - 2][i - 1];
    }
  }

  for (int i = 1; i <= n - 1; ++i) {
    for (int j = 1; j <= i; ++j) {
      double s33 = a[i - 1][j - 1];
      for (int k = 1; k <= n - i; ++k) {
        s33 += a[i + k - 1][j - 1] * a[i - 1][i + k - 1];
      }
      a[i - 1][j - 1] = s33;
    }
    for (int j = 1; j <= n - i; ++j) {
      double s34 = 0.;
      for (int k = j; k <= n - i; ++k) {
        s34 += a[i + k - 1][i + j - 1] * a[i - 1][i + k - 1];
      }
      a[i - 1][i + j - 1] = s34;
    }
  }

  int nxch = ir[n - 1];
  if (nxch == 0) return;

  for (int m = 1; m <= nxch; ++m) {
    int k = nxch - m + 1;
    int ij = ir[k - 1];
    int i = ij / 4096;
    int j = ij % 4096;
    for (k = 1; k <= n; ++k) {
      double ti = a[k - 1][i - 1];
      a[k - 1][i - 1] = a[k - 1][j - 1];
      a[k - 1][j - 1] = ti;
    }
  }
}

//   ******************************************************************
//
//   REPLACES B BY THE SOLUTION X OF A*X=B, AND REPLACES A BY ITS IN-
//   VERSE.
//
//   n            ORDER OF THE SQUARE MATRIX IN ARRAY A.
//   A            (DOUBLE PRECISION) TWO-DIMENSIONAL ARRAY CONTAINING
//                AN n BY n MATRIX.
//
//   IFAIL        OUTPUT PARAMETER.   IFAIL= 0 ... NORMAL EXIT.
//                                    IFAIL=-1 ... SINGULAR MATRIX.
//
//   B            (DOUBLE PRECISION) ONE-DIMENSIONAL ARRAY
//
//   CALLS ... DFACT, DFINV.
//
//   ******************************************************************

void Deqinv(const int n, std::vector<std::vector<double> >& a, int& ifail,
            std::vector<double>& b) {

  // Test for parameter errors.
  if (n < 1) {
    ifail = +1;
    return;
  }

  ifail = 0;
  int jfail = 0;

  double det = 0.;
  if (n > 3) {
    // n > 3 cases. Factorize matrix, invert and solve system.
    std::vector<int> ir(n, 0);
    Dfact(n, a, ir, ifail, det, jfail);
    if (ifail != 0) return;
    Dfeqn(n, a, ir, b);
    Dfinv(n, a, ir);
  } else if (n == 3) {
    // n = 3 case. Compute cofactors.
    const double c11 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    const double c12 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    const double c13 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    const double c21 = a[2][1] * a[0][2] - a[2][2] * a[0][1];
    const double c22 = a[2][2] * a[0][0] - a[2][0] * a[0][2];
    const double c23 = a[2][0] * a[0][1] - a[2][1] * a[0][0];
    const double c31 = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    const double c32 = a[0][2] * a[1][0] - a[0][0] * a[1][2];
    const double c33 = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    const double t1 = fabs(a[0][0]);
    const double t2 = fabs(a[1][0]);
    const double t3 = fabs(a[2][0]);

    // Set temp = pivot and det = pivot * det.
    double temp = 0.;
    if (t1 >= t2) {
      if (t3 >= t1) {
        // Pivot is A31
        temp = a[2][0];
        det = c23 * c12 - c22 * c13;
      } else {
        // Pivot is A11
        temp = a[0][0];
        det = c22 * c33 - c23 * c32;
      }
    } else {
      if (t3 >= t2) {
        // Pivot is A31
        temp = a[2][0];
        det = c23 * c12 - c22 * c13;
      } else {
        // Pivot is A21
        temp = a[1][0];
        det = c13 * c32 - c12 * c33;
      }
    }

    // Set elements of inverse in A.
    if (det == 0.) {
      ifail = -1;
      return;
    }
    const double s = temp / det;
    a[0][0] = s * c11;
    a[0][1] = s * c21;
    a[0][2] = s * c31;
    a[1][0] = s * c12;
    a[1][1] = s * c22;
    a[1][2] = s * c32;
    a[2][0] = s * c13;
    a[2][1] = s * c23;
    a[2][2] = s * c33;

    // Replace b by Ainv * b.
    const double b1 = b[0];
    const double b2 = b[1];
    b[0] = a[0][0] * b1 + a[0][1] * b2 + a[0][2] * b[2];
    b[1] = a[1][0] * b1 + a[1][1] * b2 + a[1][2] * b[2];
    b[2] = a[2][0] * b1 + a[2][1] * b2 + a[2][2] * b[2];
  } else if (n == 2) {
    // n = 2 case by Cramer's rule.
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if (det == 0.) {
      ifail = -1;
      return;
    }
    const double s = 1. / det;
    const double c11 = s * a[1][1];
    a[0][1] = -s * a[0][1];
    a[1][0] = -s * a[1][0];
    a[1][1] = s * a[0][0];
    a[0][0] = c11;

    const double b1 = b[0];
    b[0] = c11 * b1 + a[0][1] * b[1];
    b[1] = a[1][0] * b1 + a[1][1] * b[1];
  } else {
    // n = 1 case.
    if (a[0][0] == 0.) {
      ifail = -1;
      return;
    }
    a[0][0] = 1. / a[0][0];
    b[0] = a[0][0] * b[0];
  }
}

void Cfact(const int n, std::vector<std::vector<std::complex<double> > >& a,
           std::vector<int>& ir, int& ifail, std::complex<double>& det,
           int& jfail) {
  constexpr double g1 = 1.e-19;
  constexpr double g2 = 1.e-19;

  ifail = jfail = 0;

  int nxch = 0;
  det = std::complex<double>(1., 0.);

  for (int j = 1; j <= n; ++j) {
    int k = j;
    double p = std::max(fabs(real(a[j - 1][j - 1])), fabs(imag(a[j - 1][j - 1])));
    if (j == n) {
      if (p <= 0.) {
        det = std::complex<double>(0., 0.);
        ifail = -1;
        jfail = 0;
        return;
      }
      det *= a[j - 1][j - 1];
      a[j - 1][j - 1] = std::complex<double>(1., 0.) / a[j - 1][j - 1];
      const double t = std::max(fabs(real(det)), fabs(imag(det)));
      if (t < g1) {
        det = std::complex<double>(0., 0.);
        if (jfail == 0) jfail = -1;
      } else if (t > g2) {
        det = std::complex<double>(1., 0.);
        if (jfail == 0) jfail = +1;
      }
      continue;
    }
    for (int i = j + 1; i <= n; ++i) {
      double q = std::max(fabs(real(a[i - 1][j - 1])), fabs(imag(a[i - 1][j - 1])));
      if (q <= p) continue;
      k = i;
      p = q;
    }
    if (k != j) {
      for (int l = 1; l <= n; ++l) {
        const auto tf = a[j - 1][l - 1];
        a[j - 1][l - 1] = a[k - 1][l - 1];
        a[k - 1][l - 1] = tf;
      }
      ++nxch;
      ir[nxch - 1] = j * 4096 + k;
    } else if (p <= 0.) {
      det = std::complex<double>(0., 0.);
      ifail = -1;
      jfail = 0;
      return;
    }
    det *= a[j - 1][j - 1];
    a[j - 1][j - 1] = 1. / a[j - 1][j - 1];
    const double t = std::max(fabs(real(det)), fabs(imag(det)));
    if (t < g1) {
      det = std::complex<double>(0., 0.);
      if (jfail == 0) jfail = -1;
    } else if (t > g2) {
      det = std::complex<double>(1., 0.);
      if (jfail == 0) jfail = +1;
    }
    for (k = j + 1; k <= n; ++k) {
      auto s11 = -a[j - 1][k - 1];
      auto s12 = -a[k - 1][j];
      if (j == 1) {
        a[j - 1][k - 1] = -s11 * a[j - 1][j - 1];
        a[k - 1][j] = -(s12 + a[j - 1][j] * a[k - 1][j - 1]);
        continue;
      }
      for (int i = 1; i <= j - 1; ++i) {
        s11 += a[i - 1][k - 1] * a[j - 1][i - 1];
        s12 += a[i - 1][j] * a[k - 1][i - 1];
      }
      a[j - 1][k - 1] = -s11 * a[j - 1][j - 1];
      a[k - 1][j] = -a[j - 1][j] * a[k - 1][j - 1] - s12;
    }
  }

  if (nxch % 2 != 0) det = -det;
  if (jfail != 0) det = std::complex<double>(0., 0.);
  ir[n - 1] = nxch;
}

void Cfinv(const int n, std::vector<std::vector<std::complex<double> > >& a,
           std::vector<int>& ir) {

  if (n <= 1) return;
  a[1][0] = -a[1][1] * a[0][0] * a[1][0];
  a[0][1] = -a[0][1];
  if (n > 2) {
    for (int i = 3; i <= n; ++i) {
      for (int j = 1; j <= i - 2; ++j) {
        auto s31 = std::complex<double>(0., 0.);
        auto s32 = a[j - 1][i - 1];
        for (int k = j; k <= i - 2; ++k) {
          s31 += a[k - 1][j - 1] * a[i - 1][k - 1];
          s32 += a[j - 1][k] * a[k][i - 1];
        }
        a[i - 1][j - 1] =
            -a[i - 1][i - 1] * (s31 + a[i - 2][j - 1] * a[i - 1][i - 2]);
        a[j - 1][i - 1] = -s32;
      }
      a[i - 1][i - 2] = -a[i - 1][i - 1] * a[i - 2][i - 2] * a[i - 1][i - 2];
      a[i - 2][i - 1] = -a[i - 2][i - 1];
    }
  }

  for (int i = 1; i <= n - 1; ++i) {
    for (int j = 1; j <= i; ++j) {
      auto s33 = a[i - 1][j - 1];
      for (int k = 1; k <= n - i; ++k) {
        s33 += a[i + k - 1][j - 1] * a[i - 1][i + k - 1];
      }
      a[i - 1][j - 1] = s33;
    }
    for (int j = 1; j <= n - i; ++j) {
      std::complex<double> s34(0., 0.);
      for (int k = j; k <= n - i; ++k) {
        s34 += a[i + k - 1][i + j - 1] * a[i - 1][i + k - 1];
      }
      a[i - 1][i + j - 1] = s34;
    }
  }

  int nxch = ir[n - 1];
  if (nxch == 0) return;

  for (int m = 1; m <= nxch; ++m) {
    int k = nxch - m + 1;
    int ij = ir[k - 1];
    int i = ij / 4096;
    int j = ij % 4096;
    for (k = 1; k <= n; ++k) {
      const auto ti = a[k - 1][i - 1];
      a[k - 1][i - 1] = a[k - 1][j - 1];
      a[k - 1][j - 1] = ti;
    }
  }
}

//    ******************************************************************
//
//     REPLACES A BY ITS INVERSE.
//
//     (PARAMETERS AS FOR CEQINV.)
//
//     CALLS ... CFACT, CFINV.
//
//     ******************************************************************

void Cinv(const int n, std::vector<std::vector<std::complex<double> > >& a,
          int& ifail) {

  std::complex<double> det(0., 0.);
  // Test for parameter errors.
  if (n < 1) {
    ifail = 1;
    return;
  }

  ifail = 0;
  int jfail = 0;

  if (n > 3) {
    // n > 3 cases. Factorize matrix and invert.
    std::vector<int> ir(n, 0);
    Cfact(n, a, ir, ifail, det, jfail);
    if (ifail != 0) return;
    Cfinv(n, a, ir);
  } else if (n == 3) {
    // n = 3 case. Compute cofactors.
    const auto c11 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    const auto c12 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    const auto c13 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    const auto c21 = a[2][1] * a[0][2] - a[2][2] * a[0][1];
    const auto c22 = a[2][2] * a[0][0] - a[2][0] * a[0][2];
    const auto c23 = a[2][0] * a[0][1] - a[2][1] * a[0][0];
    const auto c31 = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    const auto c32 = a[0][2] * a[1][0] - a[0][0] * a[1][2];
    const auto c33 = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    const double t1 = fabs(real(a[0][0])) + fabs(imag(a[0][0]));
    const double t2 = fabs(real(a[1][0])) + fabs(imag(a[1][0]));
    const double t3 = fabs(real(a[2][0])) + fabs(imag(a[2][0]));

    // Set temp = pivot and det = pivot * det.
    std::complex<double> temp(0., 0.);
    if (t1 >= t2) {
      if (t3 >= t1) {
        // Pivot is A31
        temp = a[2][0];
        det = c23 * c12 - c22 * c13;
      } else {
        // Pivot is A11
        temp = a[0][0];
        det = c22 * c33 - c23 * c32;
      }
    } else {
      if (t3 >= t2) {
        // Pivot is A31
        temp = a[2][0];
        det = c23 * c12 - c22 * c13;
      } else {
        // Pivot is A21
        temp = a[1][0];
        det = c13 * c32 - c12 * c33;
      }
    }
    // Set elements of inverse in A.
    if (real(det) == 0. && imag(det) == 0.) {
      ifail = -1;
      return;
    }
    const auto s = temp / det;
    a[0][0] = s * c11;
    a[0][1] = s * c21;
    a[0][2] = s * c31;
    a[1][0] = s * c12;
    a[1][1] = s * c22;
    a[1][2] = s * c32;
    a[2][0] = s * c13;
    a[2][1] = s * c23;
    a[2][2] = s * c33;
  } else if (n == 2) {
    // n=2 case by Cramer's rule.
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if (real(det) == 0. && imag(det) == 0.) {
      ifail = -1;
      return;
    }
    const auto s = std::complex<double>(1., 0.) / det;
    const auto c11 = s * a[1][1];
    a[0][1] = -s * a[0][1];
    a[1][0] = -s * a[1][0];
    a[1][1] = s * a[0][0];
    a[0][0] = c11;
  } else {
    // n = 1 case.
    if (real(a[0][0]) == 0. && imag(a[0][0]) == 0.) {
      ifail = -1;
      return;
    }
    a[0][0] = std::complex<double>(1., 0.) / a[0][0];
  }
}

// Numerical integration using 15-point Gauss-Kronrod algorithm
// Origin: QUADPACK
double GaussKronrod15(double (*f)(const double), const double a,
                      const double b) {
  // Abscissae of the 15-point Kronrod rule
  // xGK[1], xGK[3], ... abscissae of the 7-point Gauss rule
  // xGK[0], xGK[2], ... abscissae which are optimally added
  //                     to the 7-point Gauss rule
  constexpr double xGK[8] = {9.914553711208126e-01, 9.491079123427585e-01,
                             8.648644233597691e-01, 7.415311855993944e-01,
                             5.860872354676911e-01, 4.058451513773972e-01,
                             2.077849550078985e-01, 0.0e+00};
  // Weights of the 15-point Kronrod rule
  constexpr double wGK[8] = {2.293532201052922e-02, 6.309209262997855e-02,
                             1.047900103222502e-01, 1.406532597155259e-01,
                             1.690047266392679e-01, 1.903505780647854e-01,
                             2.044329400752989e-01, 2.094821410847278e-01};
  // Weights of the 7-point Gauss rule
  constexpr double wG[4] = {1.294849661688697e-01, 2.797053914892767e-01,
                            3.818300505051189e-01, 4.179591836734694e-01};

  // Mid-point of the interval
  const double center = 0.5 * (a + b);
  // Half-length of the interval
  const double halfLength = 0.5 * (b - a);

  double fC = f(center);
  // Result of the 7-point Gauss formula
  double resG = fC * wG[3];
  // Result of the 15-point Kronrod formula
  double resK = fC * wGK[7];

  for (int j = 0; j < 3; ++j) {
    const int i = j * 2 + 1;
    // Abscissa
    const double x = halfLength * xGK[i];
    // Function value
    const double fSum = f(center - x) + f(center + x);
    resG += wG[j] * fSum;
    resK += wGK[i] * fSum;
  }

  for (int j = 0; j < 4; ++j) {
    const int i = j * 2;
    const double x = halfLength * xGK[i];
    const double fSum = f(center - x) + f(center + x);
    resK += wGK[i] * fSum;
  }

  return resK * halfLength;
}

double Divdif(const std::vector<double>& f, const std::vector<double>& a,
              const int nn, const double x, const int mm) {

  double t[20], d[20];

  // Check the arguments.
  if (nn < 2) {
    std::cerr << "Divdif: Array length < 2.\n";
    return 0.;
  }
  if (mm < 1) {
    std::cerr << "Divdif: Interpolation order < 1.\n";
    return 0.;
  }

  // Deal with the case that X is located at first or last point.
  const double tol = 1.e-6 * (fabs(a[0]) + fabs(a[nn - 1]));
  if (fabs(x - a[0]) < tol) return f[0];
  if (fabs(x - a[nn - 1]) < tol) return f[nn - 1];

  // Find subscript IX of X in array A.
  constexpr int mmax = 10;
  const int m = std::min({mm, mmax, nn - 1});
  const int mplus = m + 1;
  int ix = 0;
  int iy = nn + 1;
  if (a[0] > a[nn - 1]) {
    // Search decreasing arguments.
    do {
      const int mid = (ix + iy) / 2;
      if (x > a[mid - 1]) {
        iy = mid;
      } else {
        ix = mid;
      }
    } while (iy - ix > 1);
  } else {
    // Search increasing arguments.
    do {
      const int mid = (ix + iy) / 2;
      if (x < a[mid - 1]) {
        iy = mid;
      } else {
        ix = mid;
      }
    } while (iy - ix > 1);
  }
  //  Copy reordered interpolation points into (T[I],D[I]), setting
  //  EXTRA to True if M+2 points to be used.
  int npts = m + 2 - (m % 2);
  int ip = 0;
  int l = 0;
  do {
    const int isub = ix + l;
    if ((1 > isub) || (isub > nn)) {
      // Skip point.
      npts = mplus;
    } else {
      // Insert point.
      ip++;
      t[ip - 1] = a[isub - 1];
      d[ip - 1] = f[isub - 1];
    }
    if (ip < npts) {
      l = -l;
      if (l >= 0) ++l;
    }
  } while (ip < npts);

  const bool extra = npts != mplus;
  // Replace d by the leading diagonal of a divided-difference table,
  // supplemented by an extra line if EXTRA is True.
  for (l = 1; l <= m; l++) {
    if (extra) {
      const int isub = mplus - l;
      d[m + 1] = (d[m + 1] - d[m - 1]) / (t[m + 1] - t[isub - 1]);
    }
    int i = mplus;
    for (int j = l; j <= m; j++) {
      const int isub = i - l;
      d[i - 1] = (d[i - 1] - d[i - 2]) / (t[i - 1] - t[isub - 1]);
      i--;
    }
  }
  // Evaluate the Newton interpolation formula at X, averaging two values
  // of last difference if EXTRA is True.
  double sum = d[mplus - 1];
  if (extra) {
    sum = 0.5 * (sum + d[m + 1]);
  }
  int j = m;
  for (l = 1; l <= m; l++) {
    sum = d[j - 1] + (x - t[j - 1]) * sum;
    j--;
  }
  return sum;
}

bool Boxin2(const std::vector<std::vector<double> >& value,
            const std::vector<double>& xAxis, const std::vector<double>& yAxis,
            const int nx, const int ny, const double x, const double y,
            double& f, const int iOrder) {
  //-----------------------------------------------------------------------
  //   BOXIN2 - Interpolation of order 1 and 2 in an irregular rectangular
  //            2-dimensional grid.
  //-----------------------------------------------------------------------
  int iX0 = 0, iX1 = 0;
  int iY0 = 0, iY1 = 0;
  std::array<double, 3> fX;
  std::array<double, 3> fY;
  f = 0.;
  // Ensure we are in the grid.
  if ((xAxis[nx - 1] - x) * (x - xAxis[0]) < 0 ||
      (yAxis[ny - 1] - y) * (y - yAxis[0]) < 0) {
    std::cerr << "Boxin2: Point not in the grid; no interpolation.\n";
    return false;
  }
  // Make sure we have enough points.
  if (iOrder < 0 || iOrder > 2) {
    std::cerr << "Boxin2: Incorrect order; no interpolation.\n";
    return false;
  } else if (nx < 1 || ny < 1) {
    std::cerr << "Boxin2: Incorrect number of points; no interpolation.\n";
    return false;
  }
  if (iOrder == 0 || nx <= 1) {
    // Zeroth order interpolation in x.
    // Find the nearest node.
    double dist = fabs(x - xAxis[0]);
    int iNode = 0;
    for (int i = 1; i < nx; i++) {
      if (fabs(x - xAxis[i]) < dist) {
        dist = fabs(x - xAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iX0 = iNode;
    iX1 = iNode;
    // Establish the shape functions.
    fX = {1., 0., 0.};
  } else if (iOrder == 1 || nx <= 2) {
    // First order interpolation in x.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nx; i++) {
      if ((xAxis[i - 1] - x) * (x - xAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    const double x0 = xAxis[iGrid - 1];
    const double x1 = xAxis[iGrid];
    // Ensure there won't be divisions by zero.
    if (x1 == x0) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double xL = (x - x0) / (x1 - x0);
    // Set the summing range.
    iX0 = iGrid - 1;
    iX1 = iGrid;
    // Set the shape functions.
    fX = {1. - xL, xL, 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in x.
    // Find the nearest node and the grid segment.
    double dist = fabs(x - xAxis[0]);
    int iNode = 0;
    for (int i = 1; i < nx; ++i) {
      if (fabs(x - xAxis[i]) < dist) {
        dist = fabs(x - xAxis[i]);
        iNode = i;
      }
    }
    // Find the nearest fitting 2x2 matrix.
    int iGrid = std::max(1, std::min(nx - 2, iNode));
    const double x0 = xAxis[iGrid - 1];
    const double x1 = xAxis[iGrid];
    const double x2 = xAxis[iGrid + 1];
    // Ensure there won't be divisions by zero.
    if (x2 == x0) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute the alpha and local coordinate for this grid segment.
    const double xAlpha = (x1 - x0) / (x2 - x0);
    const double xL = (x - x0) / (x2 - x0);
    // Ensure there won't be divisions by zero.
    if (xAlpha <= 0 || xAlpha >= 1) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Set the summing range.
    iX0 = iGrid - 1;
    iX1 = iGrid + 1;
    // Set the shape functions.
    const double xL2 = xL * xL;
    fX[0] = xL2 / xAlpha - xL * (1. + xAlpha) / xAlpha + 1.;
    fX[1] = (xL2 - xL) / (xAlpha * xAlpha - xAlpha);
    fX[2] = (xL2 - xL * xAlpha) / (1. - xAlpha);
  }
  if (iOrder == 0 || ny <= 1) {
    // Zeroth order interpolation in y.
    // Find the nearest node.
    double dist = fabs(y - yAxis[0]);
    int iNode = 0;
    for (int i = 1; i < ny; i++) {
      if (fabs(y - yAxis[i]) < dist) {
        dist = fabs(y - yAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iY0 = iNode;
    iY1 = iNode;
    // Establish the shape functions.
    fY = {1., 0., 0.};
  } else if (iOrder == 1 || ny <= 2) {
    // First order interpolation in y.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < ny; ++i) {
      if ((yAxis[i - 1] - y) * (y - yAxis[i]) >= 0) {
        iGrid = i;
      }
    }
    const double y0 = yAxis[iGrid - 1];
    const double y1 = yAxis[iGrid];
    // Ensure there won't be divisions by zero.
    if (y1 == y0) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double yL = (y - y0) / (y1 - y0);
    // Set the summing range.
    iY0 = iGrid - 1;
    iY1 = iGrid;
    // Set the shape functions.
    fY = {1. - yL, yL, 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in y.
    // Find the nearest node and the grid segment.
    double dist = fabs(y - yAxis[0]);
    int iNode = 0;
    for (int i = 1; i < ny; ++i) {
      if (fabs(y - yAxis[i]) < dist) {
        dist = fabs(y - yAxis[i]);
        iNode = i;
      }
    }
    // Find the nearest fitting 2x2 matrix.
    int iGrid = std::max(1, std::min(ny - 2, iNode));
    const double y0 = yAxis[iGrid - 1];
    const double y1 = yAxis[iGrid];
    const double y2 = yAxis[iGrid + 1];
    // Ensure there won't be divisions by zero.
    if (y2 == y0) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute the alpha and local coordinate for this grid segment.
    const double yAlpha = (y1 - y0) / (y2 - y0);
    const double yL = (y - y0) / (y2 - y0);
    // Ensure there won't be divisions by zero.
    if (yAlpha <= 0 || yAlpha >= 1) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Set the summing range.
    iY0 = iGrid - 1;
    iY1 = iGrid + 1;
    // Set the shape functions.
    const double yL2 = yL * yL;
    fY[0] = yL2 / yAlpha - yL * (1. + yAlpha) / yAlpha + 1.;
    fY[1] = (yL2 - yL) / (yAlpha * yAlpha - yAlpha);
    fY[2] = (yL2 - yL * yAlpha) / (1. - yAlpha);
  }

  // Sum the shape functions.
  for (int i = iX0; i <= iX1; ++i) {
    for (int j = iY0; j <= iY1; ++j) {
      f += value[i][j] * fX[i - iX0] * fY[j - iY0];
    }
  }
  return true;
}

bool Boxin3(const std::vector<std::vector<std::vector<double> > >& value,
            const std::vector<double>& xAxis, const std::vector<double>& yAxis,
            const std::vector<double>& zAxis, const int nx, const int ny,
            const int nz, const double xx, const double yy, const double zz,
            double& f, const int iOrder) {
  //-----------------------------------------------------------------------
  //   BOXIN3 - interpolation of order 1 and 2 in an irregular rectangular
  //            3-dimensional grid.
  //-----------------------------------------------------------------------
  int iX0 = 0, iX1 = 0;
  int iY0 = 0, iY1 = 0;
  int iZ0 = 0, iZ1 = 0;
  std::array<double, 4> fX;
  std::array<double, 4> fY;
  std::array<double, 4> fZ;

  f = 0.;
  // Ensure we are in the grid.
  const double x = std::min(std::max(xx, std::min(xAxis[0], xAxis[nx - 1])),
                            std::max(xAxis[0], xAxis[nx - 1]));
  const double y = std::min(std::max(yy, std::min(yAxis[0], yAxis[ny - 1])),
                            std::max(yAxis[0], yAxis[ny - 1]));
  const double z = std::min(std::max(zz, std::min(zAxis[0], zAxis[nz - 1])),
                            std::max(zAxis[0], zAxis[nz - 1]));

  // Make sure we have enough points.
  if (iOrder < 0 || iOrder > 2) {
    std::cerr << "Boxin3: Incorrect order; no interpolation.\n";
    return false;
  } else if (nx < 1 || ny < 1 || nz < 1) {
    std::cerr << "Boxin3: Incorrect number of points; no interpolation.\n";
    return false;
  }
  if (iOrder == 0 || nx == 1) {
    // Zeroth order interpolation in x.
    // Find the nearest node.
    double dist = fabs(x - xAxis[0]);
    int iNode = 0;
    for (int i = 1; i < nx; i++) {
      if (fabs(x - xAxis[i]) < dist) {
        dist = fabs(x - xAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iX0 = iNode;
    iX1 = iNode;
    // Establish the shape functions.
    fX = {1., 0., 0., 0.};
  } else if (iOrder == 1 || nx == 2) {
    // First order interpolation in x.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nx; i++) {
      if ((xAxis[i - 1] - x) * (x - xAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    const double x0 = xAxis[iGrid - 1];
    const double x1 = xAxis[iGrid];
    // Ensure there won't be divisions by zero.
    if (x1 == x0) {
      std::cerr << "Boxin3: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double xL = (x - x0) / (x1 - x0);
    // Set the summing range.
    iX0 = iGrid - 1;
    iX1 = iGrid;
    // Set the shape functions.
    fX = {1. - xL, xL, 0., 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in x.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nx; i++) {
      if ((xAxis[i - 1] - x) * (x - xAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    if (iGrid == 1) {
      iX0 = iGrid - 1;
      iX1 = iGrid + 1;
      const double x0 = xAxis[iX0];
      const double x1 = xAxis[iX0 + 1];
      const double x2 = xAxis[iX0 + 2];
      if (x0 == x1 || x0 == x2 || x1 == x2) {
        std::cerr << "Boxin3: One or more grid points in x coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fX[0] = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
      fX[1] = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
      fX[2] = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
    } else if (iGrid == nx - 1) {
      iX0 = iGrid - 2;
      iX1 = iGrid;
      const double x0 = xAxis[iX0];
      const double x1 = xAxis[iX0 + 1];
      const double x2 = xAxis[iX0 + 2];
      if (x0 == x1 || x0 == x2 || x1 == x2) {
        std::cerr << "Boxin3: One or more grid points in x coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fX[0] = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
      fX[1] = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
      fX[2] = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
    } else {
      iX0 = iGrid - 2;
      iX1 = iGrid + 1;
      const double x0 = xAxis[iX0];
      const double x1 = xAxis[iX0 + 1];
      const double x2 = xAxis[iX0 + 2];
      const double x3 = xAxis[iX0 + 3];
      if (x0 == x1 || x0 == x2 || x0 == x3 || 
          x1 == x2 || x1 == x3 || x2 == x3) {
        std::cerr << "Boxin3: One or more grid points in x coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      // Compute the local coordinate for this grid segment.
      const double xL = (x - x1) / (x2 - x1);
      fX[0] = ((x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))) * (1. - xL);
      fX[1] = ((x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))) * (1. - xL) +
              ((x - x2) * (x - x3) / ((x1 - x2) * (x1 - x3))) * xL;
      fX[2] = ((x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1))) * (1. - xL) + 
              ((x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3))) * xL;
      fX[3] = ((x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2))) * xL;
    }
  }

  if (iOrder == 0 || ny == 1) {
    // Zeroth order interpolation in y.
    // Find the nearest node.
    double dist = fabs(y - yAxis[0]);
    int iNode = 0;
    for (int i = 1; i < ny; i++) {
      if (fabs(y - yAxis[i]) < dist) {
        dist = fabs(y - yAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iY0 = iNode;
    iY1 = iNode;
    // Establish the shape functions.
    fY = {1., 0., 0., 0.};
  } else if (iOrder == 1 || ny == 2) {
    // First order interpolation in y.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < ny; i++) {
      if ((yAxis[i - 1] - y) * (y - yAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    // Ensure there won't be divisions by zero.
    const double y0 = yAxis[iGrid - 1];
    const double y1 = yAxis[iGrid];
    if (y1 == y0) {
      std::cerr << "Boxin3: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double yL = (y - y0) / (y1 - y0);
    // Set the summing range.
    iY0 = iGrid - 1;
    iY1 = iGrid;
    // Set the shape functions.
    fY = {1. - yL, yL, 0., 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in y.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < ny; i++) {
      if ((yAxis[i - 1] - y) * (y - yAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    if (iGrid == 1) {
      iY0 = iGrid - 1;
      iY1 = iGrid + 1;
      const double y0 = yAxis[iY0];
      const double y1 = yAxis[iY0 + 1];
      const double y2 = yAxis[iY0 + 2];
      if (y0 == y1 || y0 == y2 || y1 == y2) {
        std::cerr << "Boxin3: One or more grid points in y coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fY[0] = (y - y1) * (y - y2) / ((y0 - y1) * (y0 - y2));
      fY[1] = (y - y0) * (y - y2) / ((y1 - y0) * (y1 - y2));
      fY[2] = (y - y0) * (y - y1) / ((y2 - y0) * (y2 - y1));
    } else if (iGrid == ny - 1) {
      iY0 = iGrid - 2;
      iY1 = iGrid;
      const double y0 = yAxis[iY0];
      const double y1 = yAxis[iY0 + 1];
      const double y2 = yAxis[iY0 + 2];
      if (y0 == y1 || y0 == y2 || y1 == y2) {
        std::cerr << "Boxin3: One or more grid points in y coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fY[0] = (y - y1) * (y - y2) / ((y0 - y1) * (y0 - y2));
      fY[1] = (y - y0) * (y - y2) / ((y1 - y0) * (y1 - y2));
      fY[2] = (y - y0) * (y - y1) / ((y2 - y0) * (y2 - y1));
    } else {
      iY0 = iGrid - 2;
      iY1 = iGrid + 1;
      const double y0 = yAxis[iY0];
      const double y1 = yAxis[iY0 + 1];
      const double y2 = yAxis[iY0 + 2];
      const double y3 = yAxis[iY0 + 3];
      if (y0 == y1 || y0 == y2 || y0 == y3 || 
          y1 == y2 || y1 == y3 || y2 == y3) {
        std::cerr << "Boxin3: One or more grid points in y coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      // Compute the local coordinate for this grid segment.
      const double yL = (y - y1) / (y2 - y1);
      fY[0] = ((y - y1) * (y - y2) / ((y0 - y1) * (y0 - y2))) * (1. - yL);
      fY[1] = ((y - y0) * (y - y2) / ((y1 - y0) * (y1 - y2))) * (1. - yL) +
              ((y - y2) * (y - y3) / ((y1 - y2) * (y1 - y3))) * yL;
      fY[2] = ((y - y0) * (y - y1) / ((y2 - y0) * (y2 - y1))) * (1. - yL) +
              ((y - y1) * (y - y3) / ((y2 - y1) * (y2 - y3))) * yL;
      fY[3] = ((y - y1) * (y - y2) / ((y3 - y1) * (y3 - y2))) * yL;
    }
  }

  if (iOrder == 0 || nz == 1) {
    // Zeroth order interpolation in z.
    // Find the nearest node.
    double dist = fabs(z - zAxis[0]);
    int iNode = 0;
    for (int i = 1; i < nz; i++) {
      if (fabs(z - zAxis[i]) < dist) {
        dist = fabs(z - zAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iZ0 = iNode;
    iZ1 = iNode;
    // Establish the shape functions.
    fZ = {1., 0., 0., 0.};
  } else if (iOrder == 1 || nz == 2) {
    // First order interpolation in z.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nz; i++) {
      if ((zAxis[i - 1] - z) * (z - zAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    const double z0 = zAxis[iGrid - 1];
    const double z1 = zAxis[iGrid];
    // Ensure there won't be divisions by zero.
    if (z1 == z0) {
      std::cerr << "Boxin3: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double zL = (z - z0) / (z1 - z0);
    // Set the summing range.
    iZ0 = iGrid - 1;
    iZ1 = iGrid;
    // Set the shape functions.
    fZ = {1. - zL, zL, 0., 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in z.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nz; i++) {
      if ((zAxis[i - 1] - z) * (z - zAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    if (iGrid == 1) {
      iZ0 = iGrid - 1;
      iZ1 = iGrid + 1;
      const double z0 = zAxis[iZ0];
      const double z1 = zAxis[iZ0 + 1];
      const double z2 = zAxis[iZ0 + 2];
      if (z0 == z1 || z0 == z2 || z1 == z2) {
        std::cerr << "Boxin3: One or more grid points in z coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fZ[0] = (z - z1) * (z - z2) / ((z0 - z1) * (z0 - z2));
      fZ[1] = (z - z0) * (z - z2) / ((z1 - z0) * (z1 - z2));
      fZ[2] = (z - z0) * (z - z1) / ((z2 - z0) * (z2 - z1));
    } else if (iGrid == nz - 1) {
      iZ0 = iGrid - 2;
      iZ1 = iGrid;
      const double z0 = zAxis[iZ0];
      const double z1 = zAxis[iZ0 + 1];
      const double z2 = zAxis[iZ0 + 2];
      if (z0 == z1 || z0 == z2 || z1 == z2) {
        std::cerr << "Boxin3: One or more grid points in z coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fZ[0] = (z - z1) * (z - z2) / ((z0 - z1) * (z0 - z2));
      fZ[1] = (z - z0) * (z - z2) / ((z1 - z0) * (z1 - z2));
      fZ[2] = (z - z0) * (z - z1) / ((z2 - z0) * (z2 - z1));
    } else {
      iZ0 = iGrid - 2;
      iZ1 = iGrid + 1;
      const double z0 = zAxis[iZ0];
      const double z1 = zAxis[iZ0 + 1];
      const double z2 = zAxis[iZ0 + 2];
      const double z3 = zAxis[iZ0 + 3];
      if (z0 == z1 || z0 == z2 || z0 == z3 || 
          z1 == z2 || z1 == z3 || z2 == z3) {
        std::cerr << "Boxin3: One or more grid points in z coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      // Compute the local coordinate for this grid segment.
      const double zL = (z - z1) / (z2 - z1);
      fZ[0] = ((z - z1) * (z - z2) / ((z0 - z1) * (z0 - z2))) * (1. - zL);
      fZ[1] = ((z - z0) * (z - z2) / ((z1 - z0) * (z1 - z2))) * (1. - zL) +
              ((z - z2) * (z - z3) / ((z1 - z2) * (z1 - z3))) * zL;
      fZ[2] = ((z - z0) * (z - z1) / ((z2 - z0) * (z2 - z1))) * (1. - zL) +
              ((z - z1) * (z - z3) / ((z2 - z1) * (z2 - z3))) * zL;
      fZ[3] = ((z - z1) * (z - z2) / ((z3 - z1) * (z3 - z2))) * zL;
    }
  }

  for (int i = iX0; i <= iX1; ++i) {
    for (int j = iY0; j <= iY1; ++j) {
      for (int k = iZ0; k <= iZ1; ++k) {
        f += value[i][j][k] * fX[i - iX0] * fY[j - iY0] * fZ[k - iZ0];
      }
    }
  }
  return true;
}
}
}
