#include <cmath>
#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Polygon.hh"
#include "Garfield/SolidTube.hh"

namespace Garfield {

SolidTube::SolidTube(const double cx, const double cy, const double cz,
                     const double rmin, const double rmax, const double lz)
    : Solid(cx, cy, cz, "SolidTube"),
      m_rMin(rmin),
      m_rMax(rmax),
      m_lZ(lz) {}

SolidTube::SolidTube(const double cx, const double cy, const double cz,
                     const double rmin, const double rmax, const double lz,
                     const double dx, const double dy, const double dz)
    : SolidTube(cx, cy, cz, rmin, rmax, lz) {
  SetDirection(dx, dy, dz);
}

SolidTube::SolidTube(const double cx, const double cy, const double cz,
                     const double r, const double lz)
    : SolidTube(cx, cy, cz, 0., r, lz) {}

SolidTube::SolidTube(const double cx, const double cy, const double cz,
                     const double r, const double lz, const double dx,
                     const double dy, const double dz)
    : SolidTube(cx, cy, cz, 0., r, lz, dx, dy, dz) {}

bool SolidTube::IsInside(const double x, const double y, const double z,
                         const bool exact) const {
  // Transform the point to local coordinates.
  double u = x, v = y, w = z;
  ToLocal(x, y, z, u, v, w);

  if (fabs(w) > m_lZ) return false;

  // TODO!
  const double r = sqrt(u * u + v * v);
  if (r < m_rMin || r > m_rMax) return false;
  return true;
}

bool SolidTube::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                               double& xmax, double& ymax, double& zmax) const {
  if (m_cTheta == 1. && m_cPhi == 1.) {
    xmin = m_cX - m_rMax;
    xmax = m_cX + m_rMax;
    ymin = m_cY - m_rMax;
    ymax = m_cY + m_rMax;
    zmin = m_cZ - m_lZ;
    zmax = m_cZ + m_lZ;
    return true;
  }

  const double dd = sqrt(m_rMax * m_rMax + m_lZ * m_lZ);
  xmin = m_cX - dd;
  xmax = m_cX + dd;
  ymin = m_cY - dd;
  ymax = m_cY + dd;
  zmin = m_cZ - dd;
  zmax = m_cZ + dd;
  return true;
}

void SolidTube::SetInnerRadius(const double rmin) {
  if (rmin <= 0.) {
    std::cerr << "SolidTube::SetInnerRadius: Radius must be > 0.\n";
    return;
  }
  if (rmin >= m_rMax) {
    std::cerr << "SolidTube::SetInnerRadius:\n"
              << "    Inner radius must be smaller than outer radius.\n";
    return;
  }
  m_rMin = rmin;
}

void SolidTube::SetOuterRadius(const double rmax) {
  if (rmax <= 0.) {
    std::cerr << "SolidTube::SetOuterRadius: Radius must be > 0.\n";
    return;
  }
  if (rmax <= m_rMin) {
    std::cerr << "SolidTube::SetOuterRadius:\n"
              << "    Outer radius must be greater than inner radius.\n";
    return;
  }
  m_rMax = rmax;
}

void SolidTube::SetRadius(const double r) {
  if (r <= 0.) {
    std::cerr << "SolidTube::SetRadius: Radius must be > 0.\n";
    return;
  }
  m_rMax = r;
  m_rMin = 0.;
}

void SolidTube::SetHalfLength(const double lz) {
  if (lz <= 0.) {
    std::cerr << "SolidTube::SetHalfLength: Half-length must be > 0.\n";
    return;
  }
  m_lZ = lz;
}

void SolidTube::SetSectors(const unsigned int n) {
  if (n < 1) {
    std::cerr << "SolidTube::SetSectors: Number must be > 0.\n";
    return;
  }
  m_n = n;
}

bool SolidTube::SolidPanels(std::vector<Panel>& panels) {

  const auto id = GetId();
  const unsigned int nPanels = panels.size();
  // Direction vector.
  const double fnorm = sqrt(m_dX * m_dX + m_dY * m_dY + m_dZ * m_dZ);
  if (fnorm <= 0) {
    std::cerr << "SolidTube::SolidPanels:\n"
              << "    Zero norm direction vector; no panels generated.\n";
    return false;
  }

  // Set the mean or the outer radius.
  double r = m_rMax;
  if (m_average) {
    const double alpha = Pi / (4. * (m_n - 1.));
    r = 2 * m_rMax / (1. + asinh(tan(alpha)) * cos(alpha) / tan(alpha));
  }

  const unsigned int nPoints = 4 * (m_n - 1);
  // Create the top lid.
  if (m_toplid) {
    std::vector<double> xv;
    std::vector<double> yv;
    std::vector<double> zv;
    for (unsigned int i = 1; i <= nPoints; i++) {
      const double alpha = m_rot + HalfPi * (i - 1.) / (m_n - 1.);
      // Rotate into place.
      double x, y, z;
      ToGlobal(r * cos(alpha), r * sin(alpha), m_lZ, x, y, z);
      xv.push_back(x);
      yv.push_back(y);
      zv.push_back(z);
    }
    Panel panel;
    panel.a = m_cPhi * m_sTheta;
    panel.b = m_sPhi * m_sTheta;
    panel.c = m_cTheta;
    panel.xv = xv;
    panel.yv = yv;
    panel.zv = zv;
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // Create the bottom lid.
  if (m_botlid) {
    std::vector<double> xv;
    std::vector<double> yv;
    std::vector<double> zv;
    for (unsigned int i = 1; i <= nPoints; i++) {
      const double alpha = m_rot + HalfPi * (i - 1.) / (m_n - 1.);
      // Rotate into place.
      double x, y, z;
      ToGlobal(r * cos(alpha), r * sin(alpha), -m_lZ, x, y, z);
      xv.push_back(x);
      yv.push_back(y);
      zv.push_back(z);
    }
    Panel panel;
    panel.a = -m_cPhi * m_sTheta;
    panel.b = -m_sPhi * m_sTheta;
    panel.c = -m_cTheta;
    panel.xv = xv;
    panel.yv = yv;
    panel.zv = zv;
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // Create the side panels.
  double u = r * cos(m_rot);
  double v = r * sin(m_rot);
  // Rotate into place.
  double xv0, yv0, zv0;
  ToGlobal(u, v, -m_lZ, xv0, yv0, zv0);
  double xv1, yv1, zv1;
  ToGlobal(u, v, +m_lZ, xv1, yv1, zv1);
  // Go around the cylinder.
  for (unsigned int i = 2; i <= nPoints + 1; i++) {
    // Bottom and top of the line along the axis of the cylinder.
    double alpha = m_rot + HalfPi * (i - 1.) / (m_n - 1.);
    u = r * cos(alpha);
    v = r * sin(alpha);
    // Rotate into place.
    double xv2, yv2, zv2;
    ToGlobal(u, v, +m_lZ, xv2, yv2, zv2);
    double xv3, yv3, zv3;
    ToGlobal(u, v, -m_lZ, xv3, yv3, zv3);
    // Store the plane.
    Panel panel;
    alpha = m_rot + HalfPi * (i - 1.5) / (m_n - 1.);
    const double cAlpha = cos(alpha);
    const double sAlpha = sin(alpha);
    panel.a = m_cPhi * m_cTheta * cAlpha - m_sPhi * sAlpha;
    panel.b = m_sPhi * m_cTheta * cAlpha + m_cPhi * sAlpha;
    panel.c = -m_sTheta * cAlpha;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
    // Shift the points.
    xv0 = xv3;
    yv0 = yv3;
    zv0 = zv3;
    xv1 = xv2;
    yv1 = yv2;
    zv1 = zv2;
  }
  std::cout << "SolidTube::SolidPanels: " << panels.size() - nPanels
            << " panels.\n";
  return true;
}

double SolidTube::GetDiscretisationLevel(const Panel& panel) {

  // Transform the normal vector to local coordinates.
  double u = 0., v = 0., w = 0.;
  VectorToLocal(panel.a, panel.b, panel.c, u, v, w);
  // Identify the vector.
  if (w > std::max(std::abs(u), std::abs(v))) {
    return m_dis[0];
  } else if (w < -std::max(std::abs(u), std::abs(v))) {
    return m_dis[1];
  }
  return m_dis[2];
} 

void SolidTube::Cut(const double x0, const double y0, const double z0,
                    const double xn, const double yn, const double zn,
                    std::vector<Panel>& panels) {

  // -----------------------------------------------------------------------
  //    PLACYC - Cuts cylinder with a plane.
  // -----------------------------------------------------------------------
  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  // Set the mean or the outer radius.
  double r = m_rMax;
  if (m_average) {
    const double alpha = Pi / (4. * (m_n - 1.));
    r = 2 * m_rMax / (1. + asinh(tan(alpha)) * cos(alpha) / tan(alpha));
  }
  const unsigned int nPoints = 4 * (m_n - 1);
  const double dphi = HalfPi / (m_n - 1.);
  // Go through the lines of the top and bottom lids.
  for (const auto zLid : {-m_lZ, +m_lZ}) {
    double x1, y1, z1;
    ToGlobal(r * cos(m_rot), r * sin(m_rot), zLid, x1, y1, z1);
    // Loop over the points.
    for (unsigned int i = 2; i <= nPoints + 1; ++i) {
      const double phi = m_rot + (i - 1.) * dphi;
      double x2, y2, z2;
      ToGlobal(r * cos(phi), r * sin(phi), zLid, x2, y2, z2);
      // Cut with the plane.
      double xc, yc, zc;
      if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, 
                    xn, yn, zn, xc, yc, zc)) {
        xv.push_back(xc);
        yv.push_back(yc);
        zv.push_back(zc);
      }
      // Shift the coordinates.
      x1 = x2;
      y1 = y2;
      z1 = z2;
    }
  }
  // Go through the ribs.
  for (unsigned int i = 2; i <= nPoints + 1; ++i) {
    // Bottom and top of the line along the axis of the cylinder.
    const double phi = m_rot + (i - 1.) * dphi;
    const double u = r * cos(phi);
    const double v = r * sin(phi);
    double x1, y1, z1;
    ToGlobal(u, v, +m_lZ, x1, y1, z1);
    double x2, y2, z2;
    ToGlobal(u, v, -m_lZ, x2, y2, z2);
    // Cut with the plane.
    double xc, yc, zc;
    if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
      xv.push_back(xc);
      yv.push_back(yc);
      zv.push_back(zc);
    }
  }
  // Get rid of butterflies.
  Polygon::EliminateButterflies(xv, yv, zv);

  if (xv.size() >= 3) {
    Panel panel;
    panel.a = xn;
    panel.b = yn;
    panel.c = zn;
    panel.xv = xv;
    panel.yv = yv;
    panel.zv = zv;
    panel.colour = m_colour;
    panel.volume = GetId();
    panels.push_back(std::move(panel));
  }
}

}
