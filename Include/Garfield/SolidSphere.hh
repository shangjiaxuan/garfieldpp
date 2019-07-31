#ifndef G_SOLID_SPHERE_H
#define G_SOLID_SPHERE_H

#include "Solid.hh"

namespace Garfield {

/// Sphere.

class SolidSphere : public Solid {
 public:
  /// Constructor
  SolidSphere(const double cx, const double cy, const double cz,
              const double rmin, const double rmax);
  /// Destructor
  ~SolidSphere() {}

  bool IsInside(const double x, const double y, const double z) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsSphere() const override { return true; }

  /// Set the inner radius of the sphere.
  void SetInnerRadius(const double rmin);
  /// Set the outer radius of the sphere.
  void SetOuterRadius(const double rmax);

  double GetInnerRadius() const override { return m_rMin; }
  double GetOuterRadius() const override { return m_rMax; }

  /// When calculating surface panels, the sphere is approximated by a set of
  /// parallelograms, much the same way maps are drawn. N specifies the number
  /// of meridians and also the number of parallels.
  void SetMeridians(const unsigned int n);

  bool SolidPanels(std::vector<Panel>& panels) override;

 private:
  /// Inner radius.
  double m_rMin = 0.;
  /// Outer radius
  double m_rMax;

  /// Number of meridians.
  unsigned int m_n = 10;
};
}

#endif
