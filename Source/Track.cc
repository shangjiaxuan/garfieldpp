#include <iostream>
#include <algorithm>
#include <cctype>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/Track.hh"
#include "Garfield/ViewDrift.hh"

namespace Garfield {

Track::Track() : m_mass(MuonMass) { SetBetaGamma(3.); }

void Track::SetParticle(const std::string& part) {
  std::string id = part;
  std::transform(id.begin(), id.end(), id.begin(), 
                 [](unsigned char c) -> unsigned char { 
                   return std::toupper(c);
                 });
  m_isElectron = false;
  if (id == "ELECTRON" || id == "E-") {
    m_q = -1;
    m_mass = ElectronMass;
    m_spin = 1;
    m_isElectron = true;
    m_particleName = "e-";
  } else if (id == "POSITRON" || id == "E+") {
    m_q = 1;
    m_mass = ElectronMass;
    m_spin = 1;
    m_particleName = "e+";
  } else if (id == "MUON" || id == "MU" || id == "MU-") {
    m_q = -1;
    m_mass = MuonMass;
    m_spin = 1;
    m_particleName = "mu-";
  } else if (id == "MU+") {
    m_q = 1;
    m_mass = MuonMass;
    m_spin = 1;
    m_particleName = "mu+";
  } else if (id == "PION" || id == "PI" || id == "PI-") {
    m_q = -1;
    m_mass = 139.57018e6;
    m_spin = 0;
    m_particleName = "pi-";
  } else if (id == "PI+") {
    m_q = 1;
    m_mass = 139.57018e6;
    m_spin = 0;
    m_particleName = "pi+";
  } else if (id == "KAON" || id == "K" || id == "K-") {
    m_q = -1;
    m_mass = 493.677e6;
    m_spin = 0;
    m_particleName = "K-";
  } else if (id == "K+") {
    m_q = 1;
    m_mass = 493.677e6;
    m_spin = 0;
    m_particleName = "K+";
  } else if (id == "PROTON" || id == "P") {
    m_q = 1;
    m_mass = ProtonMass;
    m_spin = 1;
    m_particleName = "p";
  } else if (id == "ANTI-PROTON" || id == "ANTIPROTON" || 
             id == "P-BAR" || id == "PBAR") {
    m_q = -1;
    m_mass = ProtonMass;
    m_spin = 1;
    m_particleName = "pbar";
  } else if (id == "DEUTERON" || id == "D") {
    m_q = 1;
    m_mass = 1875.612793e6;
    m_spin = 2;
    m_particleName = "d";
  } else if (id == "ALPHA") {
    m_q = 2;
    m_mass = 3.727379240e9;
    m_spin = 0;
    m_particleName = "alpha";
  } else {
    std::cerr << m_className << "::SetParticle:\n"
              << "    Particle " << part << " is not defined.\n";
  }
}

void Track::SetEnergy(const double e) {
  if (e <= m_mass) {
    std::cerr << m_className << "::SetEnergy:\n"
              << "    Particle energy must be greater than the mass.\n";
    return;
  }

  m_energy = e;
  const double gamma = m_energy / m_mass;
  m_beta2 = 1. - 1. / (gamma * gamma);
  m_isChanged = true;
}

void Track::SetBetaGamma(const double bg) {
  if (bg <= 0.) {
    std::cerr << m_className << "::SetBetaGamma:\n"
              << "    Particle speed must be greater than zero.\n";
    return;
  }

  const double bg2 = bg * bg;
  m_energy = m_mass * sqrt(1. + bg2);
  m_beta2 = bg2 / (1. + bg2);
  m_isChanged = true;
}

void Track::SetBeta(const double beta) {
  if (beta <= 0. || beta >= 1.) {
    std::cerr << m_className << "::SetBeta:\n"
              << "    Beta must be between zero and one.\n";
    return;
  }

  m_beta2 = beta * beta;
  m_energy = m_mass * sqrt(1. / (1. - m_beta2));
  m_isChanged = true;
}

void Track::SetGamma(const double gamma) {
  if (gamma <= 1.) {
    std::cerr << m_className << "::SetGamma:\n"
              << "    Gamma must be greater than one.\n";
    return;
  }

  m_energy = m_mass * gamma;
  m_beta2 = 1. - 1. / (gamma * gamma);
  m_isChanged = true;
}

void Track::SetMomentum(const double p) {
  if (p <= 0.) {
    std::cerr << m_className << "::SetMomentum:\n"
              << "    Particle momentum must be greater than zero.\n";
    return;
  }

  m_energy = sqrt(m_mass * m_mass + p * p);
  const double bg = p / m_mass;
  m_beta2 = bg * bg / (1. + bg * bg);
  m_isChanged = true;
}

void Track::SetKineticEnergy(const double ekin) {
  if (ekin <= 0.) {
    std::cerr << m_className << "::SetKineticEnergy:\n"
              << "    Kinetic energy must be greater than zero.\n";
    return;
  }

  m_energy = m_mass + ekin;
  const double gamma = 1. + ekin / m_mass;
  m_beta2 = 1. - 1. / (gamma * gamma);
  m_isChanged = true;
}

void Track::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }

  m_sensor = s;
}

void Track::EnablePlotting(ViewDrift* view) {
  if (!view) {
    std::cerr << m_className << "::EnablePlotting: Null pointer.\n";
    return;
  }

  m_viewer = view;
  m_usePlotting = true;
}

void Track::DisablePlotting() {
  m_usePlotting = false;
  m_viewer = nullptr;
}

void Track::PlotNewTrack(const double x0, const double y0, const double z0) {
  if (!m_usePlotting || !m_viewer) return;

  m_viewer->NewChargedParticleTrack(1, m_plotId, x0, y0, z0);
}

void Track::PlotCluster(const double x0, const double y0, const double z0) {
  if (m_plotId < 0 || !m_usePlotting || !m_viewer) {
    std::cerr << m_className << "::PlotCluster: No track set. Program bug!\n";
    return;
  }
  m_viewer->AddTrackPoint(m_plotId, x0, y0, z0);
}
}
