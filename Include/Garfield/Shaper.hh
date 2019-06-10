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
  Shaper(const int n, const double tau, const double g, 
         const std::string& shaperType){
    m_n = n;
    m_tau = tau;
    m_g = g;
    m_type = shaperType;
    m_rand = new TRandom3(time(NULL));
    m_transfer_func_sq = -1.;
  }

  /// Destructor
  ~Shaper() {}
  
  // Get transfer function value
  double Shape(double t);
  // Transfer function for a bipolar shaper
  double UnipolarShaper(double t);
  // Transfer function for a unipolar shaper
  double BipolarShaper(double t);
  // Heaviside function, t0
  double Heaviside(double t, double t0);
  // Time @ peak.
  double PeakingTime();


  // Calculate a noise current to add to the detector signal for a given ENC.
  double WhiteNoise(int enc, double tStep);
  // Calculate the number of delta pulses required for a given ENC.
  double NDeltaPulses(double q_enc, double q0, double tStep);
  // Calculate the integral of the transfer function squared.
  void CalculateTransferFuncSq(double tStep, unsigned int nTimeBins);

  private:
  std::string m_className = "Shaper";

  // Time window for signals                                                                    
  static double m_signalConversion;

  // Standard transfer function parameters.
  int m_n;
  double m_tau;
  double m_g;  
  std::string m_type;
  double m_transfer_func_sq;

  // Noise.
  TRandom3 * m_rand = nullptr;

  // Switch on/off debugging messages
  bool m_debug = false;

};
}

#endif
