#ifndef G_MEDIUM_GAAS_H
#define G_MEDIUM_GAAS_H

#include "Medium.hh"

namespace Garfield {

/// Gallium-Arsenide.

class MediumGaAs : public Medium {
 public:
  /// Constructor
  MediumGaAs();
  /// Destructor
  virtual ~MediumGaAs() {}

  bool IsSemiconductor() const override { return true; }

  void GetComponent(const unsigned int i, std::string& label, 
                    double& f) override;

  // Electron transport parameters
  bool ElectronVelocity(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        double& vx, double& vy, double& vz) override;
  bool ElectronTownsend(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        double& alpha) override;
  bool ElectronAttachment(const double ex, const double ey, const double ez,
                          const double bx, const double by, const double bz,
                          double& eta) override;
  // Hole transport parameters
  bool HoleVelocity(const double ex, const double ey, const double ez,
                    const double bx, const double by, const double bz,
                    double& vx, double& vy, double& vz) override;
  bool HoleTownsend(const double ex, const double ey, const double ez,
                    const double bx, const double by, const double bz,
                    double& alpha) override;
  bool HoleAttachment(const double ex, const double ey, const double ez,
                      const double bx, const double by, const double bz,
                      double& eta) override;

  void SetLowFieldMobility(const double mue, const double muh);

 private:
  // Band-gap energy [eV]
  // double m_bandGap = 1.42;
  // Low-field mobility
  double m_eMobility = 8.8e-6;
  double m_hMobility = 3.2e-6;
  // Hall factor
  double m_eHallFactor = 1.05;
  double m_hHallFactor = 1.25;

  // Models
  bool m_hasUserMobility = false;
};
}

#endif
