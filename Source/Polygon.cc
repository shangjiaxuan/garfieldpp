#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

#include "Garfield/Polygon.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"

namespace {

unsigned int NextPoint(const unsigned int i, const unsigned int n) {
  const unsigned int j = i + 1;
  return j < n ? j : 0;
}

unsigned int PrevPoint(const unsigned int i, const unsigned int n) {
  return i > 0 ? i - 1 : n - 1;
}

double Mag(const std::array<double, 3>& a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

std::array<double, 3> UnitVector(const std::array<double, 3>& a) {

  const double mag = Mag(a);
  if (mag < 1.e-12) return a;
  const std::array<double, 3> b = {a[0] / mag, a[1] / mag, a[2] / mag};
  return b;
}
std::array<double, 3> CrossProduct(const std::array<double, 3>& u,
                                   const std::array<double, 3>& v) {

  const std::array<double, 3> w = {u[1] * v[2] - u[2] * v[1],
                                   u[2] * v[0] - u[0] * v[2],
                                   u[0] * v[1] - u[1] * v[0]};
  return w;
}

/// Compute lambda for a point on a line (0 = start, 1 = end).
double Lambda(const double x1, const double x0, const double x2,
              const double y1, const double y0, const double y2) {
  // Segment of zero length.
  if ((x1 - x2) == 0. && (y1 - y2) == 0.) {
    std::cerr << "Lambda: Zero length segment.\n";
    return 2.;
  }

  double xl = 0.;
  const double dx1 = x0 - x1;
  const double dy1 = y0 - y1;
  const double dx2 = x0 - x2;
  const double dy2 = y0 - y2;
  if (dx1 * dx1 + dy1 * dy1 < dx2 * dx2 + dy2 * dy2) {
    // Point nearer to (x1, y1).
    if (fabs(y1 - y2) > fabs(x1 - x2)) {
      xl = dy1 / (y2 - y1);
    } else {
      xl = dx1 / (x2 - x1);
    }
  } else {
    // Point nearer to (x2, y2).
    if (fabs(y1 - y2) > fabs(x1 - x2)) {
      xl = 1. - dy2 / (y1 - y2);
    } else {
      xl = 1. - dx2 / (x1 - x2);
    }
  }
  return xl;
}

/// Determine whether a point (u, v) lies on a straight line
/// (x1, y1) to (x2, y2).
bool OnLine(const double x1, const double y1, const double x2, const double y2,
            const double u, const double v) {
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v)});
  epsx = std::max(1.e-10, epsx);
  epsy = std::max(1.e-10, epsy);

  if ((fabs(x1 - u) <= epsx && fabs(y1 - v) <= epsy) ||
      (fabs(x2 - u) <= epsx && fabs(y2 - v) <= epsy)) {
    // Point to be examined coincides with start or end.
    return true;
  } else if (fabs(x1 - x2) <= epsx && fabs(y1 - y2) <= epsy) {
    // The line (x1, y1) to (x2, y2) is in fact a point.
    return false;
  }
  double xc = 0., yc = 0.;
  if (fabs(u - x1) + fabs(v - y1) < fabs(u - x2) + fabs(v - y2)) {
    // (u, v) is nearer to (x1, y1).
    const double dx = (x2 - x1);
    const double dy = (y2 - y1);
    const double xl = ((u - x1) * dx + (v - y1) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x1;
      yc = y1;
    } else if (xl > 1.) {
      xc = x2;
      yc = y2;
    } else {
      xc = x1 + xl * dx;
      yc = y1 + xl * dy;
    }
  } else {
    // (u, v) is nearer to (x2, y2).
    const double dx = (x1 - x2);
    const double dy = (y1 - y2);
    const double xl = ((u - x2) * dx + (v - y2) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x2;
      yc = y2;
    } else if (xl > 1.) {
      xc = x1;
      yc = y1;
    } else {
      xc = x2 + xl * dx;
      yc = y2 + xl * dy;
    }
  }
  // See whether the point is on the line.
  if (fabs(u - xc) < epsx && fabs(v - yc) < epsy) {
    return true;
  }
  return false;
}

/// Determine whether the 2 straight lines (x1, y1) to (x2, y2)
/// and (u1, v1) to (u2, v2) cross at an intermediate point for both lines.
bool Crossing(const double x1, const double y1, const double x2,
              const double y2, const double u1, const double v1,
              const double u2, const double v2, double& xc, double& yc) {
  /// Matrix to compute the crossing point.
  std::array<std::array<double, 2>, 2> a;
  a[0][0] = y2 - y1;
  a[0][1] = v2 - v1;
  a[1][0] = x1 - x2;
  a[1][1] = u1 - u2;
  const double det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  // Initial values.
  xc = 0.;
  yc = 0.;
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u1), fabs(u2)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v1), fabs(v2)});
  epsx = std::max(epsx, 1.e-10);
  epsy = std::max(epsy, 1.e-10);
  // Check for a point of one line located on the other line.
  if (OnLine(x1, y1, x2, y2, u1, v1)) {
    xc = u1;
    yc = v1;
    return true;
  } else if (OnLine(x1, y1, x2, y2, u2, v2)) {
    xc = u2;
    yc = v2;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x1, y1)) {
    xc = x1;
    yc = y1;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x2, y2)) {
    xc = x2;
    yc = y2;
    return true;
  } else if (fabs(det) < epsx * epsy) {
    // Parallel, non-touching.
    return false;
  }
  // Crossing, non-trivial lines: solve crossing equations.
  const double aux = a[1][1];
  a[1][1] = a[0][0] / det;
  a[0][0] = aux / det;
  a[1][0] = -a[1][0] / det;
  a[0][1] = -a[0][1] / det;
  // Compute crossing point.
  xc = a[0][0] * (x1 * y2 - x2 * y1) + a[1][0] * (u1 * v2 - u2 * v1);
  yc = a[0][1] * (x1 * y2 - x2 * y1) + a[1][1] * (u1 * v2 - u2 * v1);
  // See whether the crossing point is on both lines.
  if (OnLine(x1, y1, x2, y2, xc, yc) && OnLine(u1, v1, u2, v2, xc, yc)) {
    // Intersecting lines.
    return true;
  }
  // Crossing point not on both lines.
  return false;
}

}

namespace Garfield {

namespace Polygon {

void Inside(const std::vector<double>& xpl, const std::vector<double>& ypl,
            const double x, const double y, bool& inside, bool& edge) {
  // Initial settings.
  inside = false;
  edge = false;
  const unsigned int npl = xpl.size();
  if (ypl.size() != npl) return;
  // Special treatment for few points.
  if (npl < 2) {
    return;
  } else if (npl == 2) {
    edge = OnLine(xpl[0], ypl[0], xpl[1], ypl[1], x, y);
    return;
  }
  // Determine the range of the data.
  const double xmin = *std::min_element(std::begin(xpl), std::end(xpl));
  const double xmax = *std::max_element(std::begin(xpl), std::end(xpl));
  const double ymin = *std::min_element(std::begin(ypl), std::end(ypl));
  const double ymax = *std::max_element(std::begin(ypl), std::end(ypl));

  // Set tolerances.
  double epsx = 1.e-8 * std::max(fabs(xmin), fabs(xmax));
  double epsy = 1.e-8 * std::max(fabs(ymin), fabs(ymax));
  epsx = std::max(epsx, 1.e-8);
  epsy = std::max(epsy, 1.e-8);

  // Ensure that we have a range.
  if (fabs(xmax - xmin) <= epsx) {
    if (y >= ymin - epsy && y <= ymax + epsy &&
        fabs(xmax + xmin - 2 * x) <= epsx) {
      edge = true;
    } else {
      edge = false;
    }
  } else if (fabs(ymax - ymin) <= epsy) {
    if (x >= xmin - epsx && x <= xmax + epsx &&
        fabs(ymax + ymin - 2 * y) <= epsy) {
      edge = true;
    } else {
      edge = false;
    }
  }
  // Choose a point at "infinity".
  double xinf = xmin - fabs(xmax - xmin);
  double yinf = ymin - fabs(ymax - ymin);

  unsigned int nIter = 0;
  bool ok = false;
  while (!ok && nIter < 100) {
    ok = true;
    // Loop over the edges counting intersections.
    unsigned int nCross = 0;
    for (unsigned int j = 0; j < npl; ++j) {
      const unsigned int jj = NextPoint(j, npl);
      // Flag points located on one of the edges.
      if (OnLine(xpl[j], ypl[j], xpl[jj], ypl[jj], x, y)) {
        edge = true;
        return;
      }
      // Count mid-line intersects.
      double xc = 0., yc = 0.;
      if (Crossing(x, y, xinf, yinf, xpl[j], ypl[j], xpl[jj], ypl[jj], xc,
                   yc)) {
        ++nCross;
      }
      // Ensure that the testing line doesn't cross a corner.
      if (OnLine(x, y, xinf, yinf, xpl[j], ypl[j])) {
        xinf = xmin - ::Garfield::RndmUniform() * fabs(xmax - xinf);
        yinf = ymin - ::Garfield::RndmUniform() * fabs(ymax - yinf);
        ok = false;
        break;
      }
    }
    if (ok) {
      // Set the INSIDE flag.
      if (nCross != 2 * (nCross / 2)) inside = true;
      return;
    }
    ++nIter;
  }

  std::cerr << "Polygon::Inside:\n    Warning. Unable to verify "
            << "whether a point is internal; setting to edge.\n";
  inside = false;
  edge = true;
}

double Area(const std::vector<double>& xp, const std::vector<double>& yp) {

  double f = 0.;
  const unsigned int n = xp.size();
  for (unsigned int i = 0; i < n; ++i) {
    const unsigned int ii = i < n - 1 ? i + 1 : 0;
    f += xp[i] * yp[ii] - xp[ii] * yp[i];  
  }
  return 0.5 * f; 
}

}

}
