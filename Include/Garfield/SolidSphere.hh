#ifndef G_SOLID_SPHERE_H
#define G_SOLID_SPHERE_H

#include "Solid.hh"

namespace Garfield {

/// Sphere.

class SolidSphere : public Solid {
 public:
  /// Constructor
  SolidSphere(const double cx, const double cy, const double cz,
              const double r);
  /// Destructor
  ~SolidSphere() {}

  bool IsInside(const double x, const double y, const double z) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsSphere() const override { return true; }

  /// Set the radius of the sphere.
  void SetRadius(const double r);

  double GetRadius() const override { return m_r; }

  /// When calculating surface panels, the sphere is approximated by a set of
  /// parallelograms, much the same way maps are drawn. N specifies the number
  /// of meridians and also the number of parallels.
  void SetMeridians(const unsigned int n);

  bool SolidPanels(std::vector<Panel>& panels) override;
  double GetDiscretisationLevel(const Panel& panel) override;

 private:
  /// Radius
  double m_r = 1.;

  /// Number of meridians.
  unsigned int m_n = 10;

  /// Discretisation level.
  double m_dis = -1.;
};
}

#endif
