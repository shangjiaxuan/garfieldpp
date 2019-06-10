#ifndef G_SHAPER_H
#define G_SHAPER_H

#include <vector>
#include <time.h>
#include <TRandom3.h>

namespace Garfield {

/// Class for signal processing

class Shaper {
 public:
  /// Constructor
  Shaper() = delete;
  Shaper(const int n, const double tau, 
         const std::string& shaperType){
    m_n = n;
    m_tau = tau;
    m_type = shaperType;
    m_rand = new TRandom3(time(NULL));
  }

  /// Destructor
  ~Shaper() {}
  
  // Get transfer function value
  double Shape(double t);
  // Transfer function for a bipolar shaper
  double UnipolarShaper(double t) const;
  // Transfer function for a unipolar shaper
  double BipolarShaper(double t);
  // Heaviside function, t0
  double Heaviside(double t, double t0) const;

  // Add noise

  // Calculate a current to add to the detector signal for a given output ENC.
  double WhiteNoise(int enc, double tStep);
  // Calculate the number of delta pulses 
  double NDeltaPulses(double q_enc, double q0, double tStep, unsigned int nTimeBins);
  double TransferFuncSq(double tStep, unsigned int nTimeBins);

  private:
  std::string m_className = "Shaper";

  static double m_signalConversion;

  // Transfer function
  bool m_unipolarShaper = false;
  bool m_bipolarShaper  = false;
  double (*m_fTransfer)(double t) = nullptr;
  std::vector<double> m_transferFunctionTimes;
  std::vector<double> m_transferFunctionValues;

  // Standard transfer function parameters.
  int m_n;
  double m_tau;
  std::string m_type;

  // Noise.
  TRandom3 * m_rand = nullptr;

  // Switch on/off debugging messages
  bool m_debug = false;

  // Return the current shaper size
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax);

  double InterpolateTransferFunctionTable(const double t) const;
};
}

#endif
