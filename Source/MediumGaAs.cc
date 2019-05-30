#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MediumGaAs.hh"
#include "Garfield/Random.hh"

namespace Garfield {

MediumGaAs::MediumGaAs() : Medium() {
  m_className = "MediumGaAs";
  m_name = "GaAs";

  SetTemperature(300.);
  SetDielectricConstant(12.9);
  SetAtomicNumber(32);
  SetAtomicWeight(72.32);
  SetMassDensity(5.317);

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = false;

  m_w = 4.35;
  m_fano = 0.1;
}

void MediumGaAs::GetComponent(const unsigned int i, std::string& label,
                              double& f) {
  if (i == 0) {
    label = "Ga";
    f = 0.5;
  } else if (i == 1) {
    label = "As";
    f = 0.5;
  }
}

bool MediumGaAs::ElectronVelocity(const double ex, const double ey,
                                  const double ez, const double bx,
                                  const double by, const double bz, double& vx,
                                  double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (!m_eVelE.empty()) {
    // Interpolation in user table.
    return Medium::ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }
  // Calculate the mobility
  double mu = -m_eMobility;
  const double b = sqrt(bx * bx + by * by + bz * bz);
  if (b < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    // Hall mobility
    const double muH = m_eHallFactor * mu;
    const double mu2 = muH * muH
    const double eb = bx * ex + by * ey + bz * ez;
    const double f = muH / (1. + mu2 * b * b);
    // Compute the drift velocity using the Langevin equation.
    vx = f * (ex + muH * (ey * bz - ez * by) + mu2 * bx * eb);
    vy = f * (ey + muH * (ez * bx - ex * bz) + mu2 * by * eb);
    vz = f * (ez + muH * (ex * by - ey * bx) + mu2 * bz * eb);
  }
  return true;
}

bool MediumGaAs::ElectronTownsend(const double ex, const double ey,
                                  const double ez, const double bx,
                                  const double by, const double bz,
                                  double& alpha) {
  alpha = 0.;
  if (!m_eAlp.empty()) {
    // Interpolation in user table.
    return Medium::ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  return false;
}

bool MediumGaAs::ElectronAttachment(const double ex, const double ey,
                                    const double ez, const double bx,
                                    const double by, const double bz,
                                    double& eta) {
  eta = 0.;
  if (!m_eAtt.empty()) {
    // Interpolation in user table.
    return Medium::ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return true;
}

bool MediumGaAs::HoleVelocity(const double ex, const double ey, const double ez,
                              const double bx, const double by, const double bz,
                              double& vx, double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (!m_hVelE.empty()) {
    // Interpolation in user table.
    return Medium::HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }
  // Calculate the mobility
  double mu = m_hMobility;
  const double b = sqrt(bx * bx + by * by + bz * bz);
  if (b < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    // Hall mobility
    const double muH = m_hHallFactor * mu;
    const double mu2 = muH * muH;
    const double eb = bx * ex + by * ey + bz * ez;
    const double f = mu / (1. + mu2 * b * b);
    // Compute the drift velocity using the Langevin equation.
    vx = f * (ex + muH * (ey * bz - ez * by) + mu2 * bx * eb);
    vy = f * (ey + muH * (ez * bx - ex * bz) + mu2 * by * eb);
    vz = f * (ez + muH * (ex * by - ey * bx) + mu2 * bz * eb);
  }
  return true;
}

bool MediumGaAs::HoleTownsend(const double ex, const double ey, const double ez,
                              const double bx, const double by, const double bz,
                              double& alpha) {
  alpha = 0.;
  if (!m_hAlp.empty()) {
    // Interpolation in user table.
    return Medium::HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  return false;
}

bool MediumGaAs::HoleAttachment(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& eta) {
  eta = 0.;
  if (!m_hAtt.empty()) {
    // Interpolation in user table.
    return Medium::HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return true;
}

void MediumGaAs::SetLowFieldMobility(const double mue, const double muh) {
  if (mue <= 0. || muh <= 0.) {
    std::cerr << m_className << "::SetLowFieldMobility:\n"
              << "    Mobility must be greater than zero.\n";
    return;
  }

  m_eMobility = mue;
  m_hMobility = muh;
  m_hasUserMobility = true;
  m_isChanged = true;
}
}
