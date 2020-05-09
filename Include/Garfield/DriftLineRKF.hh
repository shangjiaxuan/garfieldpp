#ifndef G_DRIFTLINE_RKF_H
#define G_DRIFTLINE_RKF_H

#include <string>
#include <vector>
#include <array>

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

  /// Set the sensor.
  void SetSensor(Sensor* s);

  /// Switch on drift line plotting.
  void EnablePlotting(ViewDrift* view);
  /// Switch off drift line plotting.
  void DisablePlotting();

  /// Switch on/off calculation of induced currents (default: enabled).
  void EnableSignalCalculation(const bool on = true) { m_doSignal = on; }
  /// Set the number of points to be used when averaging the delayed 
  /// signal vector over a time bin in the Sensor class.
  /// The averaging is done with a \f$2\times navg + 1\f$ point 
  /// Newton-Raphson integration. Default: 2. 
  void SetSignalAveragingOrder(const unsigned int navg) { m_navg = navg; }

  /// Set the accuracy of the Runge Kutta Fehlberg drift line integration.
  void SetIntegrationAccuracy(const double eps);
  /// Set the maximum step size that is allowed. By default, there is no limit. 
  void SetMaximumStepSize(const double ms);
  /// Do not apply an upper limit to the step size that is allowed.
  void UnsetMaximumStepSize() { m_useStepSizeLimit = false; }
  /// Request (or not) the drift line calculation to be aborted if the 
  /// drift line makes a bend sharper than 90 degrees.
  void RejectKinks(const bool on = true) { m_rejectKinks = on; }

  /// Set multiplication factor for the signal induced by electrons.
  void SetElectronSignalScalingFactor(const double scale) { m_scaleE = scale; }
  /// Set multiplication factor for the signal induced by holes.
  void SetHoleSignalScalingFactor(const double scale) { m_scaleH = scale; }
  /// Set multiplication factor for the signal induced by ions.
  void SetIonSignalScalingFactor(const double scale) { m_scaleI = scale; }

  /// Enable/disable simulation electron multiplication (default: on).
  void EnableAvalanche(const bool on = true) { m_doAvalanche = on; }
  /// Enable/disable simulation of the ion tail (default: off).
  void EnableIonTail(const bool on = true) { m_doIonTail = on; }
  /// Do not randomize the avalanche size.
  void SetGainFluctuationsFixed(const double gain = -1.);
  /// Sample the avalanche size from a Polya distribution with 
  /// shape parameter theta.
  void SetGainFluctuationsPolya(const double theta, const double mean = -1.);

  /// Simulate the drift line of an electron with a given starting point.
  bool DriftElectron(const double x0, const double y0, const double z0,
                     const double t0);
  /// Simulate the drift line of a hole with a given starting point.
  bool DriftHole(const double x0, const double y0, const double z0,
                 const double t0);
  /// Simulate the drift line of an ion with a given starting point.
  bool DriftIon(const double x0, const double y0, const double z0,
                const double t0);
  /// Simulate the drift line of an electron with a given starting point,
  /// assuming that it has positive charge.
  bool DriftPositron(const double x0, const double y0, const double z0,
                     const double t0);
  /// Simulate the drift line of an ion with a given starting point,
  /// assuming that it has negative charge.
  bool DriftNegativeIon(const double x0, const double y0, const double z0,
                        const double t0);

  /// Print the trajectory of the most recent drift line.
  void PrintDriftLine() const;
  /// Get the end point and status flag of the most recent drift line.
  void GetEndPoint(double& x, double& y, double& z, double& t, int& st) const;
  /// Get the number of points of the most recent drift line.
  unsigned int GetNumberOfDriftLinePoints() const { return m_x.size(); }
  /// Get the coordinates and time of a point along the most recent drift line.
  void GetDriftLinePoint(const unsigned int i, double& x, double& y, double& z,
                         double& t) const;

  /// Compute the sigma of the arrival time distribution for the current 
  /// drift line by integrating the longitudinal diffusion coefficient.
  double GetArrivalTimeSpread(const double eps = 1.e-4);
  /// Compute the multiplication factor for the current drift line.
  double GetGain(const double eps = 1.e-4);
  /// Compute the attachment loss factor for the current drift line.
  double GetLoss(const double eps = 1.e-4);
  double GetDriftTime() const { return m_t.empty() ? 0. : m_t.back(); }

  void GetAvalancheSize(double& ne, double& ni) const { ne = m_nE; ni = m_nI; }

  /** Compute an electric field line.
    * \param xi,yi,zi starting point
    * \param xf points along the field line
    * \param electron flag to set the direction in which to follow the field
    */ 
  bool FieldLine(const double xi, const double yi, const double zi,
                 std::vector<std::array<float, 3> >& xl,
                 const bool electron = true);

  void EnableDebugging() { m_debug = true; }
  void DisableDebugging() { m_debug = false; }

 private:
  std::string m_className = "DriftLineRKF";

  // Pointer to sensor.
  Sensor* m_sensor = nullptr;

  enum class Particle { Electron = 0, Ion, Hole, Positron, NegativeIon };
  // Type of particle (of the most current drift line).
  Particle m_particle = Particle::Electron;

  // Maximum allowed step size.
  double m_maxStepSize = 0.;
  // Precision of the stepping algorithm.
  double m_accuracy = 1.e-8;
  // Flag to reject bends > 90 degrees or not.
  bool m_rejectKinks = true;
  // Flag to apply a cut on the maximum allowed step size or not.
  bool m_useStepSizeLimit = false;

  // Pointer to the drift viewer.
  ViewDrift* m_view = nullptr;

  // Points along the current drift line.
  std::vector<std::array<double, 3> > m_x;
  // Times corresponding to the points along the current drift line.
  std::vector<double> m_t;
  // Status flag of the current drift line.
  int m_status = 0;

  // Flag whether to calculate induced signals or not.
  bool m_doSignal = true;
  // Averaging order used when projecting the signal on the time bins.
  unsigned int m_navg = 2;

  // Scaling factor for electron signals.
  double m_scaleE = 1.;
  // Scaling factor for hole signals.
  double m_scaleH = 1.;
  // Scaling factor for ion signals.
  double m_scaleI = 1.;

  // Flag wether to simulate electron multiplication or not.
  bool m_doAvalanche = true;
  enum class GainFluctuations { None = 0, Polya };
  // Model to be used for randomizing the avalanche size.
  GainFluctuations m_gainFluctuations = GainFluctuations::None; 
  // Polya shape parameter.
  double m_theta = 0.;
  // Mean avalanche size (only used if < 1).
  double m_gain = -1.;

  // Flag whether to simulate the ion tail or not.
  bool m_doIonTail = false;
  // Avalanche size.
  double m_nE = 0., m_nI = 0.;

  bool m_debug = false;

  // Calculate a drift line starting at a given position.
  bool DriftLine(const double x0, const double y0, const double z0,
                 const double t0, const Particle particle,
                 std::vector<double>& ts, 
                 std::vector<std::array<double, 3> >& xs, int& status);
  bool Avalanche(const Particle particle, 
                 const std::vector<std::array<double, 3> >& xs,
                 std::vector<double>& ne, std::vector<double>& ni, 
                 double& scale);
  bool AddIonTail(const std::vector<double>& te,
                  const std::vector<std::array<double, 3> >& xe,
                  const std::vector<double>& ni, const double scale);
  // Compute electric and magnetic field at a given position.
  int GetField(const std::array<double, 3>& x,
               double& ex, double& ey, double& ez,
               double& bx, double& by, double& bz,
               Medium*& medium) const;
  // Calculate transport parameters for the respective particle type.
  bool GetVelocity(const std::array<double, 3>& x, const Particle particle,
                   std::array<double, 3>& v, int& status) const;
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

  // Terminate a drift line at the edge of a boundary.
  bool Terminate(const std::array<double, 3>& xx0,
                 const std::array<double, 3>& xx1,
                 const Particle particle, 
                 std::vector<double>& ts, 
                 std::vector<std::array<double, 3> >& xs);

  // Drift a particle to a wire
  bool DriftToWire(const double xw, const double yw, const double rw,
                   const Particle particle, 
                   std::vector<double>& ts,
                   std::vector<std::array<double, 3> >& xs, int& stat);

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
  void ComputeSignal(const Particle particle, const double scale,
                     const std::vector<double>& ts,
                     const std::vector<std::array<double, 3> >& xs,
                     const std::vector<double>& ne) const;

  // Terminate a field line.
  void Terminate(const std::array<double, 3>& xx0,
                 const std::array<double, 3>& xx1,
                 std::vector<std::array<float, 3> >& xs);
};
}

#endif
