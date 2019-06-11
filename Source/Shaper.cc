#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Numerics.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Shaper.hh"

namespace Garfield {

double Shaper::m_signalConversion = ElementaryCharge;

double Shaper::Shape(double t) {
  if (m_type == "unipolar")
    return UnipolarShaper(t);
  else if (m_type == "bipolar")
    return BipolarShaper(t);
  else {
    std::cerr << m_className << "::Shape: Shaper type not known.\n";
    return 0;
  }
}

double Shaper::UnipolarShaper(double t){
  double f = exp(m_n) * pow( (t/(m_n*m_tau)), m_n) * exp(-t/m_tau) * Heaviside(t, 0.);
  return m_g * f;
}

double Shaper::BipolarShaper(double t){
  double r = m_n - sqrt(m_n);
  double f = exp(r)  / sqrt(m_n) * (m_n - t / m_tau) * pow(t/(r*t),m_n-1) * exp(-t/m_tau) * Heaviside(t, 0.);
  return m_g * f;
}

double Shaper::Heaviside(double t, double t0){
  if (t < t0)
    return 0;
  else if (fabs(t - t0) < pow(10,-20.))
    return 0.5;
  else
    return 1;
}

double Shaper::PeakingTime() {
  if (m_type == "unipolar")
    return m_tau * m_n;
  else if (m_type == "bipolar")
    return m_tau * (m_n - sqrt(m_n));
  return 0.;
}

double Shaper::WhiteNoise(int enc, double tStep){
  if (m_transfer_func_sq < 0) {
    std::cerr << m_className << "::WhiteNoise: Transfer function integral not calculated. \n";
  }
  // Convert the desired output sigma
  //double sigv = m_g * enc;

  // Convert the ENC --> charge (fC)
  double q_enc = enc * m_signalConversion;
  // Calculate the number of delta pulses we need to add in this bin
  int npulse = NDeltaPulses(q_enc, m_signalConversion, tStep);
  return npulse;
// Convert to current
//   return npulse * m_signalConversion / tStep; 
}

double Shaper::NDeltaPulses(double q_enc, double q0, double tStep){
  // Mean frequency of random delta pulses to model noise.
  double nu = pow(q_enc,2) / pow(q0,2) * (1/m_transfer_func_sq); 
  return m_rand->Poisson(nu * tStep) - (nu * tStep);
}

void Shaper::CalculateTransferFuncSq(double tStep, unsigned int nTimeBins){
  double integral = 0.;
  double t_temp = 0.;
  for (unsigned int it = 0; it < nTimeBins; it++){
    t_temp = (it + 0.5) * tStep;
    // Normalize transfer function to unity
    integral += pow(Shape(t_temp)/m_g,2.) * tStep;
  }
  m_transfer_func_sq = integral;
}
}
