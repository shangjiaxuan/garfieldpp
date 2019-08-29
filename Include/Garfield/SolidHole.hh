#ifndef G_SOLID_HOLE_H
#define G_SOLID_HOLE_H

#include "Solid.hh"

namespace Garfield {

/// Box with a cylindrical hole.

class SolidHole : public Solid {
 public:
  /// Constructor from centre, upper/lower radii, half-lengths of the box. 
  SolidHole(const double cx, const double cy, const double cz,
            const double rup, const double rlow, 
            const double lx, const double ly, const double lz);
  /// Constructor from centre, upper/lower radii, half-lengths of the box
  /// and orientation. 
  SolidHole(const double cx, const double cy, const double cz,
            const double rup, const double rlow, 
            const double lx, const double ly, const double lz,
            const double dx, const double dy, const double dz);
  /// Destructor
  ~SolidHole() {}

  bool IsInside(const double x, const double y, const double z) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsHole() const override { return true; }

  /// Set the half-length of the box along x.
  void SetHalfLengthX(const double lx);
  /// Set the half-length of the box along y.
  void SetHalfLengthY(const double ly);
  /// Set the half-length of the box along z.
  void SetHalfLengthZ(const double lz);
  /// Set the radius at z = +lz.
  void SetUpperRadius(const double r);
  /// Set the radius at z = -lz.
  void SetLowerRadius(const double r);

  double GetHalfLengthX() const override { return m_lX; }
  double GetHalfLengthY() const override { return m_lY; }
  double GetHalfLengthZ() const override { return m_lZ; }
  double GetUpperRadius() const override { return m_rUp; }
  double GetLowerRadius() const override { return m_rLow; }

  /// When calculating the surface panels, the hole is
  /// approximated as a polygon with a finite number of panels.
  /// The number of corners of the polygon equals \f$4(n - 1)\f$.
  /// Thus, \f$n = 2\f$ will produce a square, \f$n = 3\f$ an octagon etc.
  void SetSectors(const unsigned int n);
  /// By default, the polygon used for approximating the hole when
  /// calculating surface panels is inscribed in a circle
  /// of the specified radius. If the "average-radius" flag is activated,
  /// then the radius will be interpreted as the mean radius of the polygon
  /// that approximates the cylinder.
  void SetAverageRadius(const bool average) { m_average = average; }

  /// Return the order of the approximating polygon.
  unsigned int GetSectors() const { return m_n; }
  /// Return the state of the "average-radius" flag. 
  bool GetAverage() const { return m_average; }
 
  bool SolidPanels(std::vector<Panel>& panels) override;

 private:
  // Upper and lower radius
  double m_rUp;
  double m_rLow;
  // Half-lengths
  double m_lX;
  double m_lY;
  double m_lZ;

  /// Number of sectors
  unsigned int m_n = 2;
  /// Average chord over the sectors.
  bool m_average = false;
};
}

#endif