#ifndef G_DRIFTLINE_RKF_H
#define G_DRIFTLINE_RKF_H

#include <string>
#include <vector>
#include <array>

#include "GeometryBase.hh"
#include "Medium.hh"
#include "Sensor.hh"
#include "ViewDrift.hh"

namespace Garfield {

/// Calculation of drift lines based on macroscopic transport coefficients
/// using Runge-Kutta-Fehlberg integration.

class DriftLineRKF {
 public:
  /// Constructor
  DriftLineRKF();
  /// Destructor
  ~DriftLineRKF() {}

  void SetSensor(Sensor* s);

  void EnablePlotting(ViewDrift* view);
  void DisablePlotting();

  /// Set the accuracy of the Runge Kutta Fehlberg drift line integration.
  void SetIntegrationAccuracy(const double a);
  /// Set the maximum step size that is allowed. By default, there is no limit. 
  void SetMaximumStepSize(const double ms);
  /// Do not apply an upper limit to the step size that is allowed.
  void UnsetMaximumStepSize() { m_useStepSizeLimit = false; }
  /// Request (or not) the drift line calculation to be aborted if the 
  /// drift line makes a bend sharper than 90 degrees.
  void EnableRejectKinks(const bool on = true) { m_rejectKinks = on; }

  /// Set multiplication factor for the signal induced by electrons.
  void SetElectronSignalScalingFactor(const double scale) { m_scaleE = scale; }
  /// Set multiplication factor for the signal induced by holes.
  void SetHoleSignalScalingFactor(const double scale) { m_scaleH = scale; }
  /// Set multiplication factor for the signal induced by ions.
  void SetIonSignalScalingFactor(const double scale) { m_scaleI = scale; }

  bool DriftElectron(const double x0, const double y0, const double z0,
                     const double t0);
  bool DriftHole(const double x0, const double y0, const double z0,
                 const double t0);
  bool DriftIon(const double x0, const double y0, const double z0,
                const double t0);

  void GetEndPoint(double& x, double& y, double& z, double& t, int& st) const;
  unsigned int GetNumberOfDriftLinePoints() const { return m_path.size(); }
  void GetDriftLinePoint(const unsigned int i, double& x, double& y, double& z,
                         double& t) const;

  double GetArrivalTimeSpread(const double eps = 1.e-4);
  /// Compute the multiplication factor for the current drift line.
  double GetGain(const double eps = 1.e-4);
  /// Compute the attachment loss factor for the current drift line.
  double GetLoss(const double eps = 1.e-4);
  double GetDriftTime() const { return m_path.empty() ? 0. : m_path.back().t; }

  void EnableDebugging() { m_debug = true; }
  void DisableDebugging() { m_debug = false; }

  void EnableVerbose(const bool on = true) { m_verbose = on; }

 private:
  enum class Particle { Electron = 0, Ion, Hole };

  std::string m_className = "DriftLineRKF";

  Sensor* m_sensor = nullptr;

  Particle m_particleType;
  double m_maxStepSize = 0.;
  double m_accuracy = 1.e-8;
  bool m_rejectKinks = true;
  bool m_useStepSizeLimit = false;

  ViewDrift* m_view = nullptr;

  struct DriftPoint {
    // Position
    std::array<double, 3> x;
    // Time
    double t;
    // Integrated Townsend coefficient
    double alphaint;
  };
  std::vector<DriftPoint> m_path;
  int m_status = 0;

  double m_scaleE = 1.;
  double m_scaleH = 1.;
  double m_scaleI = 1.;

  bool m_debug = false;
  bool m_verbose = false;

  // Calculate a drift line starting at a given position.
  bool DriftLine(const double x0, const double y0, const double z0,
                 const double t0, const Particle particle,
                 std::vector<DriftPoint>& path, int& status);
  // Compute electric and magnetic field at a given position.
  int GetField(const std::array<double, 3>& x,
               double& ex, double& ey, double& ez,
               double& bx, double& by, double& bz,
               Medium*& medium) const;
  // Calculate transport parameters for the respective particle type.
  bool GetVelocity(const std::array<double, 3>& x, const Particle particle,
                   std::array<double, 3>& v, int& status);
  bool GetVelocity(const double ex, const double ey, const double ez,
                   const double bx, const double by, const double bz,
                   Medium* medium, const Particle particle,
                   std::array<double, 3>& v) const;
  bool GetDiffusion(const double ex, const double ey, const double ez,
                    const double bx, const double by, const double bz,
                    Medium* medium, const Particle particle,
                    double& dl, double& dt) const;
  bool GetAlpha(const double ex, const double ey, const double ez,
                const double bx, const double by, const double bz,
                Medium* medium, const Particle particle, double& alpha) const;
  bool GetEta(const double ex, const double ey, const double ez,
              const double bx, const double by, const double bz,
              Medium* medium, const Particle particle, double& eta) const;
  // Add a point to a drift line.
  void AddPoint(const std::array<double, 3>& x, const double t,
                std::vector<DriftPoint>& path);
  // Terminate a drift line at the edge of a boundary.
  bool Terminate(const std::array<double, 3>& xx0,
                 const std::array<double, 3>& xx1,
                 const Particle particle, std::vector<DriftPoint>& path);
  // Drift a particle to a wire
  bool DriftToWire(const double xw, const double yw, const double rw,
                   const Particle particle, std::vector<DriftPoint>& path);
  // Integrate the longitudinal diffusion over a step.
  double IntegrateDiffusion(const std::array<double, 3>& xi,
                            const std::array<double, 3>& xe,
                            const Particle particle, const double tol);
  // Integrate the Townsend coefficient over a step.
  double IntegrateAlpha(const std::array<double, 3>& xi,
                        const std::array<double, 3>& xe, 
                        const Particle particle, const double tol);
  // Integrate the attachment coefficient over a step.
  double IntegrateEta(const std::array<double, 3>& xi,
                      const std::array<double, 3>& xe, 
                      const Particle particle, const double tol);

  // Calculate the signal for the current drift line.
  void ComputeSignal(const double q, const double scale) const;
};
}

#endif
