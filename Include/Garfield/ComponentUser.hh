#ifndef G_COMPONENT_USER_H
#define G_COMPONENT_USER_H

#include <functional>

#include "ComponentBase.hh"

namespace Garfield {

/// Simple component with electric field given by a user function.

class ComponentUser : public ComponentBase {
 public:
  /// Constructor
  ComponentUser();
  /// Destructor
  ~ComponentUser() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  bool GetVoltageRange(double& vmin, double& vmax) override;
  void MagneticField(const double x, const double y, const double z, double& bx,
                     double& by, double& bz, int& status) override;
  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;
  void DelayedWeightingField(const double x, const double y, const double z,
                             const double t, double& wx, double& wy, double& wz,
                             const std::string& label) override;

  /// Set the function to be called for calculating the electric field.
  void SetElectricField(
    std::function<void(const double, const double, const double,
                       double&, double&, double&)>);
  /// Set the function to be called for calculating the potential.
  void SetPotential(
    std::function<void(const double, const double, const double, double&)>);
  /// Set the function to be called for calculating the weighting field.
  void SetWeightingField(
    std::function<void(const double, const double, const double,
                       double&, double&, double&, const std::string&)>);
  /// Set the function to be called for calculating the weighting potential.
  void SetWeightingPotential(
    std::function<void(const double, const double, const double, 
                       double&, const std::string&)>);
  /// Set the function to be called for calculating the delayed weighting field.
  void SetDelayedWeightingField(
    std::function<void(const double, const double, const double, const double,
                       double&, double&, double&, const std::string&)>);
  /// Set the function to be called for calculating the magnetic field.
  void SetMagneticField(
    std::function<void(const double, const double, const double,
                       double&, double&, double&)>);

 private:
  /// Electric field function
  std::function<void(const double, const double, const double,
                     double&, double&, double&)> m_efield;
  /// Potential function
  std::function<void(const double, const double, const double,
                     double&)> m_potential;

  /// Weighting field function
  std::function<void(const double, const double, const double, 
                     double&, double&, double&, const std::string&)> m_wfield;

  /// Weighting potential function
  std::function<void(const double, const double, const double, 
                     double&, const std::string&)> m_wpot;

  /// Delayed weighting field function
  std::function<void(const double, const double, const double, const double,
                     double&, double&, double&, const std::string&)> m_dwfield;

  /// Magnetic field function
  std::function<void(const double, const double, const double, 
                     double&, double&, double&)> m_bfield;

  /// Reset the component
  void Reset() override;
  // Verify periodicities
  void UpdatePeriodicity() override;
};
}
#endif
