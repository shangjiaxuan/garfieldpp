#include <iostream>
#include <cstdio>
#include <cmath>
#include <numeric>

#include "Garfield/DriftLineRKF.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"

namespace {

std::string PrintVec(const std::array<double, 3>& x) {

  return "(" + std::to_string(x[0]) + ", " + std::to_string(x[1]) + ", " +
         std::to_string(x[2]) + ")";
}

double Mag(const std::array<double, 3>& x) {
 
  return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

double Mag(const double x, const double y, const double z) {

  return sqrt(x * x + y * y + z * z);
}

double Mag(const double x, const double y) {

  return sqrt(x * x + y * y);
}

double Mag2(const double x, const double y) {

  return x * x + y * y;
}

}

namespace Garfield {

using Vec = std::array<double, 3>;

DriftLineRKF::DriftLineRKF() {
  m_t.reserve(1000);
  m_x.reserve(1000);
}

void DriftLineRKF::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }
  m_sensor = s;
}

void DriftLineRKF::SetIntegrationAccuracy(const double eps) {
  if (eps > 0.) {
    m_accuracy = eps;
  } else {
    std::cerr << m_className << "::SetIntegrationAccuracy:\n"
              << "    Accuracy must be greater than zero.\n";
  }
}

void DriftLineRKF::SetMaximumStepSize(const double ms) {
  if (ms > 0.) {
    m_maxStepSize = ms;
    m_useStepSizeLimit = true;
  } else {
    std::cerr << m_className << "::SetMaximumStepSize:\n"
              << "    Step size must be greater than zero.\n";
  }
}

void DriftLineRKF::EnablePlotting(ViewDrift* view) {
  if (!view) {
    std::cerr << m_className << "::EnablePlotting: Null pointer.\n";
    return;
  }
  m_view = view;
}

void DriftLineRKF::DisablePlotting() { m_view = nullptr; }

void DriftLineRKF::SetGainFluctuationsFixed(const double gain) {

  if (gain > 1.) {
    std::cout << m_className << "::SetGainFluctuationsFixed: "
              << "Avalanche size set to " << gain << ".\n";
  } else {
    std::cout << m_className << "::SetGainFluctuationsFixed:\n    "
              << "Avalanche size will be given by "
              << "the integrated Townsend coefficient.\n";
  }
  m_gain = gain; 
  m_gainFluctuations = GainFluctuations::None;
}

void DriftLineRKF::SetGainFluctuationsPolya(const double theta, 
                                            const double mean) {

  if (theta < 0.) {
    std::cerr << m_className << "::SetGainFluctuationsPolya: "
              << "Shape parameter must be >= 0.\n";
    return;
  }  
  if (mean > 1.) {
    std::cout << m_className << "::SetGainFluctuationsPolya: "
              << "Mean avalanche size set to " << mean << ".\n";
  } else {
    std::cout << m_className << "::SetGainFluctuationsPolya:\n    "
              << "Mean avalanche size will be given by "
              << "the integrated Townsend coefficient.\n";
  }
  m_gain = mean;
  m_theta = theta;
  m_gainFluctuations = GainFluctuations::Polya;
} 

bool DriftLineRKF::DriftElectron(const double x0, const double y0,
                                 const double z0, const double t0) {
  m_particle = Particle::Electron;
  if (!DriftLine(x0, y0, z0, t0, Particle::Electron, m_t, m_x, m_status)) {
    return false;
  }
  std::vector<double> ne(m_t.size(), 1.);
  std::vector<double> ni(m_t.size(), 0.);
  double scale = 1.;
  if (m_doAvalanche) Avalanche(Particle::Electron, m_x, ne, ni, scale);
  if (m_doSignal) {
    ComputeSignal(Particle::Electron, scale * m_scaleE, m_t, m_x, ne);
  }
  if (m_doAvalanche && m_doIonTail) AddIonTail(m_t, m_x, ni, scale);
  return true;
}

bool DriftLineRKF::AddIonTail(const std::vector<double>& te,
                              const std::vector<std::array<double, 3> >& xe,
                              const std::vector<double>& ni,
                              const double scale) {
  // SIGETR, SIGIOR
  const unsigned int nPoints = te.size();
  if (nPoints < 2) return false;
  if (ni.size() != nPoints) return false;
  // Loop over the electron track.
  for (unsigned int i = 1; i < nPoints; ++i) {
    // Skip points where there are no ions yet.
    if (scale * ni[i] < 1.) continue;
    // Skip also points with a negligible contribution.
    // if (scale * ni[i] < threshold * m_nI) continue;
    // Compute the ion drift line.
    const auto& x0 = xe[i];
    std::vector<double> ti;
    std::vector<std::array<double, 3> > xi;
    int stat = 0;
    if (!DriftLine(x0[0], x0[1], x0[2], te[i], Particle::Ion, ti, xi, stat)) {
      std::cerr << m_className << "::AddIonTail:\n"
                << "    Unable to obtain an ion tail; tail not added.\n";
      return false;
    }
    if (m_debug) {
      std::cout << m_className << "::AddIonTail: Origin = " << PrintVec(x0)
                << ", n = " << ti.size() << ", status = " << stat << "\n";
    }
    // Compute the contribution of the drift line to the signal.
    ComputeSignal(Particle::Ion, scale * m_scaleI * ni[i], ti, xi, {});
  }
  return true;
}

bool DriftLineRKF::DriftPositron(const double x0, const double y0,
                                 const double z0, const double t0) {

  m_particle = Particle::Positron;
  if (!DriftLine(x0, y0, z0, t0, Particle::Positron, m_t, m_x, m_status)) {
    return false;
  }
  if (m_doSignal) {
    ComputeSignal(Particle::Positron, m_scaleE, m_t, m_x, {});
  }
  return true;
}

bool DriftLineRKF::DriftHole(const double x0, const double y0, const double z0,
                             const double t0) {
  m_particle = Particle::Hole;
  if (!DriftLine(x0, y0, z0, t0, Particle::Hole, m_t, m_x, m_status)) {
    return false;
  }
  if (m_doSignal) ComputeSignal(Particle::Hole, m_scaleH, m_t, m_x, {});
  return true;
}

bool DriftLineRKF::DriftIon(const double x0, const double y0, const double z0,
                            const double t0) {
  m_particle = Particle::Ion;
  if (!DriftLine(x0, y0, z0, t0, Particle::Ion, m_t, m_x, m_status)) {
    return false;
  }
  if (m_doSignal) ComputeSignal(Particle::Ion, m_scaleI, m_t, m_x, {});
  return true;
}

bool DriftLineRKF::DriftNegativeIon(const double x0, const double y0,
                                    const double z0, const double t0) {

  m_particle = Particle::NegativeIon;
  if (!DriftLine(x0, y0, z0, t0, Particle::NegativeIon, m_t, m_x, m_status)) {
    return false;
  }
  if (m_doSignal) {
    ComputeSignal(Particle::NegativeIon, m_scaleI, m_t, m_x, {});
  }
  return true;
}

bool DriftLineRKF::DriftLine(const double xi, const double yi, const double zi,
                             const double ti, const Particle particle,
                             std::vector<double>& ts,
                             std::vector<Vec>& xs, int& flag) {

  // -----------------------------------------------------------------------
  //    DLCALC - Subroutine doing the actual drift line calculations. 
  //             The calculations are based on a Runge-Kutta-Fehlberg method
  //             which has the advantage of controlling the stepsize and the
  //             error while needing only relatively few calls to EFIELD.
  //             Full details are given in the reference quoted below.
  //    VARIABLES : H          : Current stepsize (it is in fact a delta t).
  //                HPREV      : Stores the previous value of H (comparison)
  //                INITCH     : Used for checking initial stepsize (1 = ok)
  //                Other variables such as F0, F1, F2, F3, PHII, PHIII,
  //                CI. ,CII. , BETA.. etc   are explained in the reference.
  //    REFERENCE : Stoer + Bulirsch, Einfuhrung in die Numerische
  //                Mathematic II, chapter 7, page 122, 1978, HTB, Springer.
  // -----------------------------------------------------------------------


  // Reset the drift line.
  ts.clear();
  xs.clear();
  // Reset the status flag.
  flag = StatusAlive;

  // Check if the sensor is defined.
  if (!m_sensor) {
    std::cerr << m_className << "::DriftLine: Sensor is not defined.\n";
    flag = StatusCalculationAbandoned;
    return false;
  }

  // Get the sensor's bounding box.
  double xmin = 0., xmax = 0.;
  double ymin = 0., ymax = 0.;
  double zmin = 0., zmax = 0.;
  bool bbox = m_sensor->GetArea(xmin, ymin, zmin, xmax, ymax, zmax);

  // Initialise the current position.
  Vec x0 = {xi, yi, zi};
  // Make sure the initial position is at a valid location.
  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;
  if (GetField(x0, ex, ey, ez, bx, by, bz, medium) != 0) { 
    std::cerr << m_className << "::DriftLine:\n"
              << "    No valid field at initial position.\n";
    flag = StatusLeftDriftMedium;
    return false;
  }

  // Set the numerical constants for the RKF integration.
  constexpr double c10 = 214. / 891.;
  constexpr double c11 = 1. / 33.;
  constexpr double c12 = 650. / 891.;
  constexpr double c20 = 533. / 2106.;
  constexpr double c22 = 800. / 1053.;
  constexpr double c23 = -1. / 78.;

  constexpr double b10 = 1. / 4.;
  constexpr double b20 = -189. / 800.;
  constexpr double b21 = 729. / 800.;
  constexpr double b30 = 214. / 891.;
  constexpr double b31 = 1. / 33.;
  constexpr double b32 = 650. / 891.;

  // Set the charge of the drifting particle.
  const double charge = particle == Particle::Electron ? -1 : 1;

  // Initialise the particle velocity.
  Vec v0 = {0., 0., 0.};
  if (!GetVelocity(ex, ey, ez, bx, by, bz, medium, particle, v0)) {
    std::cerr << m_className << "::DriftLine:\n"
              << "    Cannot retrieve drift velocity.\n";
    return false;
  }

  const double speed0 = Mag(v0);
  if (speed0 < Small) {
    std::cerr << m_className << "::DriftLine:\n"
              << "    Zero velocity at initial position.\n";
    return false;
  }

  // Initialise time step and previous time step.
  double h = m_accuracy / speed0;
  double hprev = h;
  double t0 = ti;

  // Set the initial point.
  ts.push_back(t0);
  xs.push_back(x0);

  int initCycle = 3;
  bool ok = true;
  while (ok) {
    // Get the velocity at the first probe point.
    Vec x1 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x1[i] += h * b10 * v0[i];
    }
    Vec v1;
    int stat = 0;
    if (!GetVelocity(x1, particle, v1, stat)) {
      flag = StatusCalculationAbandoned;
      break;
    } else if (stat != 0) {
      if (m_debug) {
        std::cout << m_className << "::DriftLine: Point 1 outside.\n";
      }
      if (!Terminate(x0, x1, particle, ts, xs)) {
        flag = StatusCalculationAbandoned;
      } else {
        flag = stat;
      }
      break;
    }
    // Get the velocity at the second probe point.
    Vec x2 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x2[i] += h * (b20 * v0[i] + b21 * v1[i]);
    }
    Vec v2;
    if (!GetVelocity(x2, particle, v2, stat)) {
      flag = StatusCalculationAbandoned;
      break;
    } else if (stat != 0) {
      if (m_debug) {
        std::cout << m_className << "::DriftLine: Point 2 outside.\n";
      }
      if (!Terminate(x0, x2, particle, ts, xs)) {
        flag = StatusCalculationAbandoned;
      } else {
        flag = stat;
      }
      break;
    }
    // Get the velocity at the third probe point.
    Vec x3 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x3[i] += h * (b30 * v0[i] + b31 * v1[i] + b32 * v2[i]);
    }
    Vec v3;
    if (!GetVelocity(x3, particle, v3, stat)) {
      flag = StatusCalculationAbandoned;
      break;
    } else if (stat != 0) {
      if (m_debug) {
        std::cout << m_className << "::DriftLine: Point 3 outside.\n";
      }
      if (!Terminate(x0, x3, particle, ts, xs)) {
        flag = StatusCalculationAbandoned;
      } else {
        flag = stat;
      }
      break;
    }
    // Check if we crossed a wire.
    double xw = 0., yw = 0., zw = 0., rw = 0.;
    if (m_sensor->IsWireCrossed(x0[0], x0[1], x0[2], 
                                x1[0], x1[1], x1[2], xw, yw, zw, true, rw) ||
        m_sensor->IsWireCrossed(x0[0], x0[1], x0[2], 
                                x2[0], x2[1], x2[2], xw, yw, zw, true, rw) ||
        m_sensor->IsWireCrossed(x0[0], x0[1], x0[2], 
                                x3[0], x3[1], x3[2], xw, yw, zw, true, rw)) {
      if (m_debug) {
        std::cout << m_className << "::DriftLine: Crossed wire.\n";
      }
      if (DriftToWire(xw, yw, rw, particle, ts, xs, stat)) {
        flag = stat;
      } else if (h > Small) {
        h *= 0.5;
        continue;
      } else {
        std::cerr << m_className << "::DriftLine: Step size too small. Stop.\n";
        flag = StatusCalculationAbandoned;
      }
      break;
    }
    // Check if we are inside the trap radius of a wire.
    if (particle != Particle::Ion) {
      if (m_sensor->IsInTrapRadius(charge, x1[0], x1[1], x1[2], xw, yw, rw)) {
        if (!DriftToWire(xw, yw, rw, particle, ts, xs, flag)) {
          flag = StatusCalculationAbandoned;
        }
        break;
      }
      if (m_sensor->IsInTrapRadius(charge, x2[0], x2[1], x2[2], xw, yw, rw)) {
        if (!DriftToWire(xw, yw, rw, particle, ts, xs, flag)) {
          flag = StatusCalculationAbandoned;
        }
        break;
      }
      if (m_sensor->IsInTrapRadius(charge, x3[0], x3[1], x3[2], xw, yw, rw)) {
        if (!DriftToWire(xw, yw, rw, particle, ts, xs, flag)) {
          flag = StatusCalculationAbandoned;
        }
        break;
      }
    } 
    // Calculate the correction terms.
    Vec phi1 = {0., 0., 0.};
    Vec phi2 = {0., 0., 0.};
    for (unsigned int i = 0; i < 3; ++i) {
      phi1[i] = c10 * v0[i] + c11 * v1[i] + c12 * v2[i];
      phi2[i] = c20 * v0[i] + c22 * v2[i] + c23 * v3[i];
    }
    // Check if the step length is valid.
    const double phi1mag = Mag(phi1);
    if (phi1mag < Small) {
      std::cerr << m_className << "::DriftLine: Step has zero length. Stop.\n";
      flag = StatusCalculationAbandoned;
      break;
    } else if (m_useStepSizeLimit && h * phi1mag > m_maxStepSize) {
      if (m_debug) {
        std::cout << m_className << "::DriftLine: Step is considered too long. "
                  << "H is reduced.\n";
      }
      h = 0.5 * m_maxStepSize / phi1mag;
      continue;
    } else if (bbox) {
      // Don't allow h to become too large in view of the time resolution.
      if (h * fabs(phi1[0]) > 0.1 * fabs(xmax - xmin) ||
          h * fabs(phi1[1]) > 0.1 * fabs(ymax - ymin)) {
        h *= 0.5;
        if (m_debug) {
          std::cout << m_className << "::DriftLine: Step is considered too long. "
                    << "H is divided by two.\n";
        }
        continue;
      }
    } else if (m_rejectKinks && xs.size() > 1) {
      const unsigned int np = xs.size();
      const auto& x = xs[np - 1];
      const auto& xp = xs[np - 2];
      if (phi1[0] * (x[0] - xp[0]) + phi1[1] * (x[1] - xp[1]) +
          phi1[2] * (x[2] - xp[2]) < 0.) {
        std::cerr << m_className << "::DriftLine: Bending angle > 90 degree.\n";
        flag = StatusSharpKink;
        break;
      }
    }
    if (m_debug) std::cout << m_className << "::DriftLine: Step size ok.\n";
    // Update the position and time.
    for (unsigned int i = 0; i < 3; ++i) x0[i] += h * phi1[i];
    t0 += h;
    // Check the new position.
    m_sensor->ElectricField(x0[0], x0[1], x0[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      // The new position is not inside a valid drift medium.
      // Terminate the drift line.
      if (m_debug) {
        std::cout << m_className << "::DriftLine: Point outside. Terminate.\n";
      }
      if (!Terminate(xs.back(), x0, particle, ts, xs)) {
        flag = StatusCalculationAbandoned;
      }
      break;
    }
    // Add the new point to the drift line.
    ts.push_back(t0);
    xs.push_back(x0);
    // Adjust the step size according to the accuracy of the two estimates.
    hprev = h;
    const double dphi = fabs(phi1[0] - phi2[0]) + fabs(phi1[1] - phi2[1]) +
                        fabs(phi1[2] - phi2[2]);
    if (dphi > 0) {
      h = sqrt(h * m_accuracy / dphi);
      if (m_debug) {
        std::cout << m_className << "::DriftLine: Adapting H to " << h << ".\n";
      }
    } else {
      h *= 2;
      if (m_debug) {
        std::cout << m_className << "::DriftLine: H increased by factor two.\n";
      }
    }
    // Make sure that H is different from zero; this should always be ok.
    if (h < Small) {
      std::cerr << m_className << "::DriftLine: Step size is zero. Stop.\n";
      flag = StatusCalculationAbandoned;
      break;
    }
    // Check the initial step size.
    if (initCycle > 0 && h < 0.2 * hprev) {
      if (m_debug) {
        std::cout << m_className << "::DriftLine: Reinitialise step size.\n";
      }
      --initCycle;
      t0 = ti;
      x0 = {xi, yi, zi};
      ts = {t0};
      xs = {x0};
      continue;
    }
    initCycle = 0;
    // Don't allow H to grow too quickly
    if (h > 10 * hprev) {
      h = 10 * hprev;
      if (m_debug) {
        std::cout << m_className << "::DriftLine: H restricted to 10 times "
                  << "the previous value.\n";
      }
    }
    // Stop in case H tends to become too small.
    if (h * (fabs(phi1[0]) + fabs(phi1[1]) + fabs(phi1[2])) < m_accuracy) {
      std::cerr << m_className << "::DriftLine: Step size has become smaller "
                << "than int. accuracy. Stop.\n";
      flag = StatusCalculationAbandoned;
      break;
    }
    // Update the velocity.
    v0 = v3;
  }
  if (m_view) {
    // If requested, add the drift line to a plot.
    int id = 0;
    const unsigned int nPoints = m_x.size();
    if (particle == Particle::Ion || particle == Particle::NegativeIon) {
      m_view->NewIonDriftLine(nPoints, id, xi, yi, zi);
    } else if (particle == Particle::Electron ||
               particle == Particle::Positron) {
      m_view->NewElectronDriftLine(nPoints, id, xi, yi, zi);
    } else if (particle == Particle::Hole) {
      m_view->NewHoleDriftLine(nPoints, id, xi, yi, zi);
    }
    for (unsigned int i = 0; i < nPoints; ++i) {
      const auto& x = m_x[i];
      m_view->SetDriftLinePoint(id, i, x[0], x[1], x[2]);
    }
  }
  if (flag == StatusCalculationAbandoned) return false;
  return true;
}

bool DriftLineRKF::Avalanche(const Particle particle,
                             const std::vector<Vec>& xs,
                             std::vector<double>& ne,
                             std::vector<double>& ni, double& scale) {

  // SIGETR
  const unsigned int nPoints = xs.size();
  if (nPoints < 2) return true;
  // Locations and weights for 6-point Gaussian integration
  constexpr double tg[6] = {-0.932469514203152028, -0.661209386466264514,
                            -0.238619186083196909, 0.238619186083196909,
                            0.661209386466264514,  0.932469514203152028};
  constexpr double wg[6] = {0.171324492379170345, 0.360761573048138608,
                            0.467913934572691047, 0.467913934572691047,
                            0.360761573048138608, 0.171324492379170345};

  ne.assign(nPoints, 1.);
  ni.assign(nPoints, 0.);
  bool start = false;
  bool overflow = false;
  // Loop over the drift line.
  for (unsigned int i = 1; i < nPoints; ++i) {
    const auto& xp = xs[i - 1];
    const auto& x = xs[i];
    const Vec dx = {x[0] - xp[0], x[1] - xp[1], x[2] - xp[2]};
    // Calculate the integrated Townsend and attachment coefficients.
    double alpsum = 0.;
    double etasum = 0.;
    for (unsigned int j = 0; j < 6; ++j) {
      const double f = 0.5 * (1. + tg[j]);
      Vec xj = xp;
      for (unsigned int k = 0; k < 3; ++k) xj[k] += f * dx[k];
      double ex = 0., ey = 0., ez = 0.;
      double bx = 0., by = 0., bz = 0.;
      Medium* medium = nullptr;
      if (GetField(xj, ex, ey, ez, bx, by, bz, medium) != 0) {
        std::cerr << m_className << "::Avalanche:\n    Invalid field at "
                  << "drift line point " << i  << ", segment " << j << ".\n";
        continue;
      }
      double alp = 0.;
      if (!GetAlpha(ex, ey, ez, bx, by, bz, medium, particle, alp)) {
        std::cerr << m_className << "::Avalanche:\n    Cannot retrieve alpha at "
                  << "drift line point " << i  << ", segment " << j << ".\n";
        continue;
      }
      double eta = 0.;
      if (!GetEta(ex, ey, ez, bx, by, bz, medium, particle, eta)) {
        std::cerr << m_className << "::Avalanche:\n    Cannot retrieve eta at "
                  << "drift line point " << i  << ", segment " << j << ".\n";
        continue;
      }
      alpsum += wg[j] * alp;
      etasum += wg[j] * eta;
    }
    alpsum *= 0.5;
    etasum *= 0.5;
    if (alpsum > 1.e-6 && !start) {
      if (m_debug) {
        std::cout << m_className << "::Avalanche: Avalanche starts at step " 
                  << i << ".\n";
      }
      start = true;
    }
    const double d = Mag(dx);
    // Update the numbers of electrons and ions.
    constexpr double expmax = 30.;
    const double logp = log(std::max(1., ne[i - 1]));
    if (logp + d * (alpsum - etasum) > expmax) {
      overflow = true;
      ne[i] = exp(expmax);
    } else { 
      ne[i] = ne[i - 1] * exp(d * (alpsum - etasum));
    }
    if (logp + d * alpsum > expmax) {
      overflow = true;
      ni[i] = exp(expmax);
    } else {
      ni[i] = ne[i - 1] * (exp(d * alpsum) - 1); 
    }
  } 
  if (overflow) {
    std::cerr << m_className << "::Avalanche:\n    "
              << "Warning: Integrating the Townsend coefficients "
              << "would lead to exponential overflow.\n    "
              << "Avalanche truncated.\n";
  }
  const double qe = ne.back();
  const double qi = std::accumulate(ni.begin(), ni.end(), 0.);
  scale = 1.;
  if (qi > 1. && 
      !(m_gainFluctuations == GainFluctuations::None && m_gain < 1.)) {
    const double gain = m_gain > 1. ? m_gain : GetGain();
    double q1 = gain;
    if (m_gainFluctuations == GainFluctuations::Polya) {
      for (unsigned int i = 0; i < 100; ++i) {
        q1 = gain * RndmPolya(m_theta);
        if (q1 >= 1.) break;
      }
      q1 = std::max(q1, 1.);
    }
    q1 *= GetLoss();
    scale = (q1 + 1.) / (qi + 1.); 
  }
  if (m_debug) {
    std::cout << m_className << "::Avalanche:\n    "
              << "Final number of electrons: " << qe << "\n    "
              << "Number of ions:            " << qi << "\n    "
              << "Charge scaling factor:     " << scale << "\n    "
              << "Avalanche development:\n Step      Electrons     Ions\n";
    for (unsigned int i = 0; i < nPoints; ++i) {
      std::printf("%6d %15.7f %15.7f\n", i, scale * ne[i], scale * ni[i]);
    }
  }
  m_nE = scale * qe;
  m_nI = scale * qi;
  return true;
}

double DriftLineRKF::GetArrivalTimeSpread(const double eps) {

  // -----------------------------------------------------------------------
  //    DLCDF1 - Routine returning the integrated diffusion coefficient of
  //             the current drift line. The routine uses an adaptive
  //             Simpson integration.
  // -----------------------------------------------------------------------

  const unsigned int nPoints = m_x.size();
  // Return straight away if there is only one point.
  if (nPoints < 2) return 0.;
  const Particle particle = m_particle;

  // First get a rough estimate of the result.
  double crude = 0.;
  double varPrev = 0.;
  for (unsigned int i = 0; i < nPoints; ++i) {
    // Get the drift velocity and diffusion coefficients at this step.
    double ex = 0., ey = 0., ez = 0.;
    double bx = 0., by = 0., bz = 0.;
    Medium* medium = nullptr;
    if (GetField(m_x[i], ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::GetArrivalTimeSpread:\n"
                << "    Invalid drift line point " << i << ".\n";
      continue;
    }
    Vec v;
    if (!GetVelocity(ex, ey, ez, bx, by, bz, medium, particle, v)) {
      std::cerr << m_className << "::GetArrivalTimeSpread:\n"
              << "    Cannot retrieve drift velocity at point " << i << ".\n";
      continue;
    }
    const double speed = Mag(v);
    if (speed < Small) {
      std::cerr << m_className << "::GetArrivalTimeSpread:\n"
                << "    Zero drift velocity at point " << i << ".\n";
      continue;
    }
    double dl = 0., dt = 0.;
    if (!GetDiffusion(ex, ey, ez, bx, by, bz, medium, particle, dl, dt)) {
      std::cerr << m_className << "::GetArrivalTimeSpread:\n"
              << "    Cannot retrieve diffusion at point " << i << ".\n";
      continue;
    }
    const double sigma = dl / speed;
    const double var = sigma * sigma;
    if (i > 0) {
      const auto& x = m_x[i];
      const auto& xPrev = m_x[i - 1];
      const double d = Mag(x[0] - xPrev[0], x[1] - xPrev[1], x[2] - xPrev[2]);
      crude += 0.5 * d * (var + varPrev);
    }
    varPrev = var;
  }
  crude = sqrt(crude);

  const double tol = eps * crude; 
  double sum = 0.;
  for (unsigned int i = 0; i < nPoints - 1; ++i) {
    sum += IntegrateDiffusion(m_x[i], m_x[i + 1], particle, tol);
  }
  return sqrt(sum);
}

double DriftLineRKF::GetGain(const double eps) {

  // -----------------------------------------------------------------------
  //    DLCTW1 - Routine returning the multiplication factor for the current
  //             drift line. The routine uses an adaptive Simpson style
  //             integration.
  // -----------------------------------------------------------------------

  const unsigned int nPoints = m_x.size();
  // Return straight away if there is only one point.
  if (nPoints < 2) return 1.;
  const Particle particle = m_particle;
  if (particle == Particle::Ion) return 1.;
  if (m_status == StatusCalculationAbandoned) return 1.;

  // First get a rough estimate of the result.
  double crude = 0.;
  double alphaPrev = 0.;
  for (unsigned int i = 0; i < nPoints; ++i) {
    // Get the Townsend coefficient at this step.
    double ex = 0., ey = 0., ez = 0.;
    double bx = 0., by = 0., bz = 0.;
    Medium* medium = nullptr;
    if (GetField(m_x[i], ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::GetGain:\n"
                << "    Invalid drift line point " << i << ".\n";
      continue;
    }
    double alpha = 0.;
    if (!GetAlpha(ex, ey, ez, bx, by, bz, medium, particle, alpha)) {
      std::cerr << m_className << "::GetGain:\n"
                << "    Cannot retrieve alpha at point " << i << ".\n";
      continue;
    }
    if (i > 0) {
      const auto& x = m_x[i];
      const auto& xPrev = m_x[i - 1];
      const double d = Mag(x[0] - xPrev[0], x[1] - xPrev[1], x[2] - xPrev[2]);
      crude += 0.5 * d * (alpha + alphaPrev);
    }
    alphaPrev = alpha;
  }
  // Stop if the rough estimate is negligibly small.
  if (crude < Small) return 1.;

  // Calculate the integration tolerance based on the rough estimate.
  const double tol = eps * crude;
  double sum = 0.;
  for (unsigned int i = 0; i < nPoints - 1; ++i) {
    sum += IntegrateAlpha(m_x[i], m_x[i + 1], particle, tol);
  }
  return exp(sum);
}

double DriftLineRKF::GetLoss(const double eps) {

  // -----------------------------------------------------------------------
  //    DLCAT1 - Routine returning the attachment losses for the current
  //             drift line. The routine uses an adaptive Simpson style
  //             integration.
  // -----------------------------------------------------------------------

  const unsigned int nPoints = m_x.size();
  // Return straight away if there is only one point.
  if (nPoints < 2) return 1.;
  const Particle particle = m_particle;
  if (particle == Particle::Ion) return 1.;
  if (m_status == StatusCalculationAbandoned) return 1.;

  // First get a rough estimate of the result.
  double crude = 0.;
  double etaPrev = 0.;
  for (unsigned int i = 0; i < nPoints; ++i) {
    // Get the Townsend coefficient at this step.
    double ex = 0., ey = 0., ez = 0.;
    double bx = 0., by = 0., bz = 0.;
    Medium* medium = nullptr;
    if (GetField(m_x[i], ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::GetLoss:\n"
                << "    Invalid drift line point " << i << ".\n";
      continue;
    }
    double eta = 0.;
    if (!GetEta(ex, ey, ez, bx, by, bz, medium, particle, eta)) {
      std::cerr << m_className << "::GetLoss:\n"
                << "    Cannot retrieve eta at point " << i << ".\n";
      continue;
    }
    if (i > 0) {
      const auto& x = m_x[i];
      const auto& xPrev = m_x[i - 1];
      const double d = Mag(x[0] - xPrev[0], x[1] - xPrev[1], x[2] - xPrev[2]);
      crude += 0.5 * d * (eta + etaPrev);
    }
    etaPrev = eta;
  }

  // Calculate the integration tolerance based on the rough estimate.
  const double tol = eps * crude;
  double sum = 0.;
  for (unsigned int i = 0; i < nPoints - 1; ++i) {
    sum += IntegrateEta(m_x[i], m_x[i + 1], particle, tol);
  }
  return exp(-sum);
}

int DriftLineRKF::GetField(const std::array<double, 3>& x,
                           double& ex, double& ey, double& ez,
                           double& bx, double& by, double& bz,
                           Medium*& medium) const {
  int status = 0;
  m_sensor->MagneticField(x[0], x[1], x[2], bx, by, bz, status);
  m_sensor->ElectricField(x[0], x[1], x[2], ex, ey, ez, medium, status);
  return status;
}

bool DriftLineRKF::GetVelocity(const std::array<double, 3>& x,
                               const Particle particle,
                               std::array<double, 3>& v, int& status) const {
  v.fill(0.);
  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;
  // Stop if we are outside a valid drift medium.
  status = GetField(x, ex, ey, ez, bx, by, bz, medium);
  if (status != 0) return true;
  if (particle == Particle::Electron) {
    return medium->ElectronVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Ion) {
    return medium->IonVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Hole) {
    return medium->HoleVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Positron) {
    const bool ok = medium->ElectronVelocity(ex, ey, ez, bx, by, bz, 
                                             v[0], v[1], v[2]);
    for (unsigned int i = 0; i < 3; ++i) v[i] *= -1;
    return ok;
  } else if (particle == Particle::NegativeIon) {
    const bool ok = medium->IonVelocity(ex, ey, ez, bx, by, bz, 
                                        v[0], v[1], v[2]);
    for (unsigned int i = 0; i < 3; ++i) v[i] *= -1;
    return ok;
  } 
  std::cerr << m_className << "::GetVelocity:\n"
            << "    Cannot retrieve drift velocity at " << PrintVec(x) << ".\n";
  return false;
}

bool DriftLineRKF::GetVelocity(const double ex, const double ey,
                               const double ez, const double bx,
                               const double by, const double bz, 
                               Medium* medium, const Particle particle,
                               std::array<double, 3>& v) const {
  v.fill(0.);
  if (particle == Particle::Electron) {
    return medium->ElectronVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Ion) {
    return medium->IonVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Hole) {
    return medium->HoleVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Positron) {
    const bool ok = medium->ElectronVelocity(ex, ey, ez, bx, by, bz, 
                                             v[0], v[1], v[2]);
    for (unsigned int i = 0; i < 3; ++i) v[i] *= -1;
    return ok; 
  } else if (particle == Particle::NegativeIon) {
    const bool ok = medium->IonVelocity(ex, ey, ez, bx, by, bz, 
                                        v[0], v[1], v[2]);
    for (unsigned int i = 0; i < 3; ++i) v[i] *= -1;
    return ok; 
  }
  return false;
}

bool DriftLineRKF::GetDiffusion(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, 
                                Medium* medium, const Particle particle,
                                double& dl, double& dt) const {
  if (particle == Particle::Electron || particle == Particle::Positron) {
    return medium->ElectronDiffusion(ex, ey, ez, bx, by, bz, dl, dt);
  } else if (particle == Particle::Ion || particle == Particle::NegativeIon) {
    return medium->IonDiffusion(ex, ey, ez, bx, by, bz, dl, dt);
  } else if (particle == Particle::Hole) {
    return medium->HoleDiffusion(ex, ey, ez, bx, by, bz, dl, dt);
  }
  return false;
}

bool DriftLineRKF::GetAlpha(const double ex, const double ey, const double ez,
                            const double bx, const double by, const double bz,
                            Medium* medium, const Particle particle,
                            double& alpha) const {
  if (particle == Particle::Electron || particle == Particle::Positron) {
    return medium->ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  } else if (particle == Particle::Hole) {
    return medium->HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  return false;
}

bool DriftLineRKF::GetEta(const double ex, const double ey, const double ez,
                          const double bx, const double by, const double bz,
                          Medium* medium, const Particle particle,
                          double& eta) const {
  if (particle == Particle::Electron) {
    return medium->ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
  } else if (particle == Particle::Hole) {
    return medium->HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return false;
}

bool DriftLineRKF::Terminate(const std::array<double, 3>& xx0,
                             const std::array<double, 3>& xx1,
                             const Particle particle,
                             std::vector<double>& ts,
                             std::vector<Vec>& xs) {

  // -----------------------------------------------------------------------
  //    DLCFMP - Terminates drift line calculation by making a last step
  //             to the boundary of the mesh or the drift medium.
  //    VARIABLES : XX0: Last point in drift medium.
  //                XX1: Estimated step, outside drift medium.
  // -----------------------------------------------------------------------

  // Check the validity of the initial point.
  Vec vv0 = {0., 0., 0.};
  int status = 0;
  if (!GetVelocity(xx0, particle, vv0, status)) {
    std::cerr << m_className << "::Terminate:\n"
              << "    Cannot retrieve initial drift velocity.\n";
    return false;
  } else if (status != 0) {
    std::cerr << m_className << "::Terminate:\n"
              << "    No valid field at initial point. Program bug!\n";
    return false;
  }
  double speed = Mag(vv0);
  if (speed < Small) {
    std::cerr << m_className << "::Terminate:\n"
              << "    Zero velocity at initial position.\n";
    return false;
  }

  // Final point just inside the medium.
  Vec x0 = xx0;
  // Final point just outside the medium.
  Vec x1 = xx1;
  // Perform some bisections.
  const unsigned int nBisections = 20;
  for (unsigned int i = 0; i < nBisections; ++i) {
    // Quit bisection when interval becomes too small.
    bool small = true;
    for (unsigned int j = 0; j < 3; ++j) {
      if (fabs(x1[j] - x0[j]) > 1.e-6 * (fabs(x0[j]) + fabs(x1[j]))) {
        small = false;
        break;
      }
    } 
    if (small) {
      if (m_debug) {
        std::cout << m_className << "::Terminate:\n"
                  << "    Bisection ended at cycle " << i << ".\n";
      }
      break; 
    }
    // Calculate the mid point and evaluate the field.
    Vec xm;
    for (unsigned int j = 0; j < 3; ++j) xm[j] = 0.5 * (x0[j] + x1[j]);
    double ex = 0., ey = 0., ez = 0.;
    Medium* medium = nullptr;
    m_sensor->ElectricField(xm[0], xm[1], xm[2], ex, ey, ez, medium, status);
    if (status == 0) {
      x0 = xm;
    } else {
      x1 = xm;
    }
  }

  // Compute drift velocity at the end of the step.
  Vec v0;
  if (!GetVelocity(x0, particle, v0, status) || status != 0) {
    std::cerr << m_className << "::Terminate:\n"
              << "    Warning: Unable to compute mean velocity at last step.\n";
  } else {
    speed = 0.5 * (speed + Mag(v0));
  }
  // Calculate the time step.
  const double dt = Mag(x0[0] - xx0[0], x0[1] - xx0[1], x0[2] - xx0[2]) / speed;
  // Add the last point, just inside the drift area.
  ts.push_back(ts.back() + dt);
  xs.push_back(x0);
  m_status = StatusLeftDriftMedium;
  return true;
}

bool DriftLineRKF::DriftToWire(const double xw, const double yw,
                               const double rw, const Particle particle,
                               std::vector<double>& ts, 
                               std::vector<Vec>& xs, int& stat) {

  // -----------------------------------------------------------------------
  //   DLCWIR - Terminates drift line calculation by making some last steps
  //            towards the surface of the wire on which it is supposed to
  //            end. The precision is controlled in order to obtain a
  //            good estimate of the total remaining drift-time.
  // -----------------------------------------------------------------------

  // Get the starting point.
  Vec x0 = xs.back();
  double t0 = ts.back() - ts.front();
  if (m_debug) {
    std::cout << m_className << "::DriftToWire:\n    Drifting from ("
              << x0[0] << ", " << x0[1] << ") to wire at ("
              << xw << ", " << yw << ") with radius " << rw << " cm.\n";
  }

  // Get the initial drift velocity.
  Vec v0;
  int status = 0;
  if (!GetVelocity(x0, particle, v0, status)) {
    std::cerr << m_className << "::DriftToWire:\n"
              << "    Cannot retrieve initial drift velocity.\n";
    return false;
  } else if (status != 0) {
    std::cerr << m_className << "::DriftToWire:\n"
              << "    No valid field at initial point. Program bug!\n";
    return false;
  }

  // Estimate the time needed to reach the wire
  // assuming a straight-line trajectory and constant velocity.
  double dt = (Mag(xw - x0[0], yw - x0[1]) - rw) / Mag(v0[0], v0[1]);
  if (m_debug) {
    std::cout << m_className << "::DriftToWire: "
              << "Estimated time needed to reach the wire: " << dt << " ns.\n";
  }

  constexpr unsigned int nMaxSplit = 10;
  unsigned int nSplit = 0;
  // Move towards the wire.
  bool onwire = false;
  const double r2 = rw * rw;
  while (!onwire && dt > 1.e-6 * t0) {
    // Calculate the estimated end point.
    Vec x1 = x0;
    for (unsigned int j = 0; j < 3; ++j) x1[j] += dt * v0[j];
    // Make sure we are not moving away from the wire.
    const double xinp0 = (x1[0] - x0[0]) * (xw - x0[0]) + 
                         (x1[1] - x0[1]) * (yw - x0[1]);
    if (xinp0 < 0.) {
      if (m_debug) {
        std::cerr << m_className << "::DriftToWire:\n"
                  << "    Particle moves away from the wire. Quit.\n";
      }
      return false;
    }
    // Check if the end point is inside the wire or the wire was crossed.
    const double xinp1 = (x0[0] - x1[0]) * (xw - x1[0]) + 
                         (x0[1] - x1[1]) * (yw - x1[1]);
    if (xinp1 < 0.) {
      if (Mag2(xw - x1[0], yw - x1[1]) <= r2) {
        onwire = true;
        if (m_debug) std::cout << m_className << "::DriftToWire: Inside.\n";
      }
    } else {
      if (m_debug) std::cout << m_className << "::DriftToWire: Wire crossed.\n";
      onwire = true;
    }
    if (onwire) {
      const double dw0 = Mag(xw - x0[0], yw - x0[1]);
      x1[0] = xw - (rw + BoundaryDistance) * (xw - x0[0]) / dw0;
      x1[1] = yw - (rw + BoundaryDistance) * (yw - x0[1]) / dw0;
      dt = Mag(x1[0] - x0[0], x1[1] - x0[1]) / Mag(v0[0], v0[1]);
      x1[2] = x0[2] + dt * v0[2];
    }
    // Calculate the drift velocity at the end point.
    Vec v1;
    if (!GetVelocity(x1, particle, v1, status)) {
      std::cerr << m_className << "::DriftToWire:\n"
                << "    Cannot retrieve drift velocity at end point. Quit.\n";
      return false;
    } else if (status != 0) {
      std::cerr << m_className << "::DriftToWire:\n"
                << "    End point is not in a valid drift medium. Quit.\n";
      return false;
    }
    // Get a point halfway between for an accuracy check.
    Vec xm;
    for (unsigned int j = 0; j < 3; ++j) xm[j] = 0.5 * (x0[j] + x1[j]);
    // Calculate the drift velocity at the mid point.
    Vec vm;
    if (!GetVelocity(xm, particle, vm, status)) {
      std::cerr << m_className << "::DriftToWire:\n"
                << "    Cannot retrieve drift velocity at mid point. Quit.\n";
      return false;
    } else if (status != 0) {
      std::cerr << m_className << "::DriftToWire:\n"
                << "    Mid point is not in a valid drift medium. Quit.\n";
      return false;
    }
    // Make sure the velocities are non-zero.
    const double speed0 = Mag(v0[0], v0[1]);
    const double speed1 = Mag(v1[0], v1[1]);
    const double speedm = Mag(vm[0], vm[1]);
    if (speed0 < Small || speed1 < Small || speedm < Small) {
      std::cerr << m_className << "::DriftToWire: Zero velocity. Stop.\n";
      return false;
    }
    const double p0 = 1. / speed0;
    const double p1 = 1. / speed1;
    const double pm = 1. / speedm;
    // Compare first and second order estimates.
    const double d = Mag(x0[0] - x1[0], x0[1] - x1[1]);
    const double tol = 1.e-4 * (1. + fabs(t0));
    if (d * fabs(p0 - 2 * pm + p1) / 3. > tol && nSplit < nMaxSplit) {
      // Accuracy was not good enough so halve the step time.
      if (m_debug) {
        std::cout << m_className << "::DriftToWire: Reducing step size.\n";
      }
      dt *= 0.5;
      onwire = false;
      ++nSplit;
      continue;
    }
    const double t1 = t0 + d * (p0 + 4 * pm + p1) / 6.;
    // Add a new point to the drift line.
    ts.push_back(ts[0] + t1);
    xs.push_back(x1);
    // Proceed to the next step.
    x0 = x1;
    t0 = t1;
    v0 = v1;
  }
  // Get the wire index (status code inside the wire).
  double ex = 0., ey = 0., ez = 0.;
  Medium* medium = nullptr;
  m_sensor->ElectricField(xw, yw, 0., ex, ey, ez, medium, stat);
  return true;
}

void DriftLineRKF::PrintDriftLine() const {

  std::cout << m_className << "::PrintDriftLine:\n";
  if (m_x.empty()) {
    std::cout << "    No drift line present.\n";
    return;
  }
  if (m_particle == Particle::Electron) {
    std::cout << "    Particle: electron\n";
  } else if (m_particle == Particle::Ion) {
    std::cout << "    Particle: ion\n";
  } else if (m_particle == Particle::Hole) {
    std::cout << "    Particle: hole\n";
  } else {
    std::cout << "    Particle: unknown\n";
  }
  std::cout << "    Status: " << m_status << "\n"
            << "  Step       time [ns]        "
            << "x [cm]          y [cm]          z [cm]\n";
  const unsigned int nPoints = m_x.size();
  for (unsigned int i = 0; i < nPoints; ++i) {
    std::printf("%6d %15.7f %15.7f %15.7f %15.7f\n", 
                i, m_t[i], m_x[i][0], m_x[i][1], m_x[i][2]);
  }
 
}

void DriftLineRKF::GetEndPoint(double& x, double& y, double& z, double& t,
                               int& stat) const {
  if (m_x.empty()) {
    x = y = z = t = 0.;
    stat = m_status;
    return;
  }
  const auto& p = m_x.back();
  x = p[0];
  y = p[1];
  z = p[2];
  t = m_t.back();
  stat = m_status;
}

void DriftLineRKF::GetDriftLinePoint(const unsigned int i, double& x, double& y,
                                     double& z, double& t) const {
  if (i >= m_x.size()) {
    std::cerr << m_className << "::GetDriftLinePoint: Index out of range.\n";
    return;
  }

  const auto& p = m_x[i];
  x = p[0];
  y = p[1];
  z = p[2];
  t = m_t[i];
}

double DriftLineRKF::IntegrateDiffusion(const std::array<double, 3>& xi,
                                        const std::array<double, 3>& xe,
                                        const Particle particle, 
                                        const double tol) {

  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;

  // Make sure the starting point is valid.
  Vec x0 = xi;
  if (GetField(x0, ex, ey, ez, bx, by, bz, medium) != 0) {
    std::cerr << m_className << "::IntegrateDiffusion: Invalid starting point "
              << PrintVec(xi) << ".\n";
    return 0.;
  }

  // Determine the drift velocity at the starting point.
  Vec v0;
  if (!GetVelocity(ex, ey, ez, bx, by, bz, medium, particle, v0)) {
    std::cerr << m_className << "::IntegrateDiffusion:\n"
              << "    Cannot retrieve drift velocity at initial point.\n";
    return 0.;
  }
  double speed0 = Mag(v0);
  if (speed0 < Small) {
    std::cerr << m_className << "::IntegrateDiffusion:\n"
              << "    Zero velocity at starting point.\n";
    return 0.;
  }
  // Determine the diffusion coefficients at the initial point.
  double dl0 = 0., dt0 = 0.;
  if (!GetDiffusion(ex, ey, ez, bx, by, bz, medium, particle, dl0, dt0)) {
    std::cerr << m_className << "::IntegrateDiffusion:\n"
              << "    Cannot retrieve diffusion at initial point.\n";
    return 0.;
  }
  const double sigma0 = dl0 / speed0;
  double var0 = sigma0 * sigma0;

  // Make sure the end point is valid.
  Vec x1 = xe;
  if (GetField(x1, ex, ey, ez, bx, by, bz, medium) != 0) {
    std::cerr << m_className << "::IntegrateDiffusion: Invalid end point "
              << PrintVec(xe) << ".\n";
    return 0.;
  }
  // Determine the drift velocity at the end point.
  Vec v1;
  if (!GetVelocity(ex, ey, ez, bx, by, bz, medium, particle, v1)) {
    std::cerr << m_className << "::IntegrateDiffusion:\n"
              << "    Cannot retrieve drift velocity at end point.\n";
    return 0.;
  }
  double speed1 = Mag(v1);
  if (speed1 < Small) {
    std::cerr << m_className << "::IntegrateDiffusion:\n"
              << "    Zero velocity at end point.\n";
    return 0.;
  }
  // Determine the diffusion coefficients at the end point.
  double dl1 = 0., dt1 = 0.;
  if (!GetDiffusion(ex, ey, ez, bx, by, bz, medium, particle, dl1, dt1)) {
    std::cerr << m_className << "::IntegrateDiffusion:\n"
              << "    Cannot retrieve diffusion at initial point.\n";
    return 0.;
  }
  double sigma1 = dl1 / speed1;
  double var1 = sigma1 * sigma1;

  double integral = 0.;
  bool done = false;
  while (!done) {
    // Check if we are close to the end point.
    if (Mag(xe[0] - x0[0], xe[1] - x0[1], xe[2] - x0[2]) < 1.e-6) {
      done = true;
      break;
    }
    const double d = Mag(x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]);
    if (d < 1.e-6) {
      // Step length has become very small.
      if (m_debug) {
        std::cout << m_className << "::IntegrateDiffusion: Small step.\n";
      }
      integral += var0 * d;
      // Proceed with the next step.
      x0 = x1;
      x1 = xe;
      continue;
    }
    // Calculate drift velocity and diffusion at the end point of the step.
    if (GetField(x1, ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::IntegrateDiffusion: Invalid end point.\n";
      break;
    }
    if (!GetVelocity(ex, ey, ez, bx, by, bz, medium, particle, v1)) {
      std::cerr << m_className << "::IntegrateDiffusion:\n"
                << "    Cannot retrieve drift velocity at end point.\n";
      break;
    }
    speed1 = Mag(v1);
    if (speed1 < Small) {
      std::cerr << m_className << "::IntegrateDiffusion:\n"
                << "    Zero drift velocity at end point.\n";
      break;
    } 
    if (!GetDiffusion(ex, ey, ez, bx, by, bz, medium, particle, dl1, dt1)) {
      std::cerr << m_className << "::IntegrateDiffusion:\n"
                << "    Cannot retrieve diffusion at end point.\n";
      break;
    }
    // Determine drift velocity and diffusion at the mid point of the step.
    Vec xm;
    for (unsigned int i = 0; i < 3; ++i) xm[i] = 0.5 * (x0[i] + x1[i]);
    if (GetField(xm, ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::IntegrateDiffusion: Invalid mid point.\n";
      break;
    }
    Vec vm;
    if (!GetVelocity(ex, ey, ez, bx, by, bz, medium, particle, vm)) {
      std::cerr << m_className << "::IntegrateDiffusion:\n"
                << "    Cannot retrieve drift velocity at mid point.\n";
      break;
    }
    const double speedm = Mag(vm);
    if (speedm < Small) {
      std::cerr << m_className << "::IntegrateDiffusion:\n"
                << "    Zero drift velocity at mid point.\n";
      break;
    } 
    double dlm = 0., dtm = 0.;
    if (!GetDiffusion(ex, ey, ez, bx, by, bz, medium, particle, dlm, dtm)) {
      std::cerr << m_className << "::IntegrateDiffusion:\n"
                << "    Cannot retrieve diffusion at mid point.\n";
      break;
    }
    sigma1 = dl1 / speed1;
    var1 = sigma1 * sigma1;
    const double sigmam = dlm / speedm;
    const double varm = sigmam * sigmam;
    // Compare first and second order estimates 
    // (integrals calculated using trapezoidal and Simpson's rule).
    if (fabs(var0 - 2 * varm + var1) * sqrt(d * 2 / (var0 + var1)) / 6. < tol) {
      // Accuracy is good enough.
      integral += d * (var0 + 4 * varm + var1) / 6.;
      // Proceed to the next step.
      x0 = x1;
      x1 = xe;
      var0 = var1;
    } else {
      // Accuracy is not good enough, so halve the step.
      x1 = xm;
      var1 = varm;
    }
  }
  return integral;
}

double DriftLineRKF::IntegrateAlpha(const std::array<double, 3>& xi, 
                                    const std::array<double, 3>& xe,
                                    const Particle particle, const double tol) {

  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;

  // Make sure the starting point is valid.
  Vec x0 = xi;
  if (GetField(x0, ex, ey, ez, bx, by, bz, medium) != 0) {
    std::cerr << m_className << "::IntegrateAlpha: Invalid starting point "
              << PrintVec(xi) << ".\n";
    return 0.;
  }
  // Determine the Townsend coefficient at the initial point.
  double alpha0 = 0.;
  if (!GetAlpha(ex, ey, ez, bx, by, bz, medium, particle, alpha0)) {
    std::cerr << m_className << "::IntegrateAlpha:\n"
              << "    Cannot retrieve Townsend coefficient at initial point.\n";
    return 0.;
  }
  // Make sure the end point is valid.
  Vec x1 = xe;
  if (GetField(x1, ex, ey, ez, bx, by, bz, medium) != 0) {
    std::cerr << m_className << "::IntegrateAlpha: Invalid end point "
              << PrintVec(xe) << ".\n";
    return 0.;
  }
  // Determine the Townsend coefficient at the end point.
  double alpha1 = 0.;
  if (!GetAlpha(ex, ey, ez, bx, by, bz, medium, particle, alpha1)) {
    std::cerr << m_className << "::IntegrateAlpha:\n"
              << "    Cannot retrieve Townsend coefficient at end point.\n";
    return 0.;
  }
  double integral = 0.;
  bool done = false;
  while (!done) {
    // Check if we are close to the end point.
    if (Mag(xe[0] - x0[0], xe[1] - x0[1], xe[2] - x0[2]) < 1.e-6) {
      done = true;
      break;
    }
    const double d = Mag(x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]);
    if (d < 1.e-6) {
      // Step length has become very small.
      if (m_debug) {
        std::cout << m_className << "::IntegrateAlpha: Small step.\n";
      }
      integral += alpha0 * d;
      // Proceed with the next step.
      x0 = x1;
      x1 = xe;
      continue;
    }
    // Calculate the Townsend coefficient at the end point of the step.
    if (GetField(x1, ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::IntegrateAlpha: Invalid end point.\n";
      break;
    }
    if (!GetAlpha(ex, ey, ez, bx, by, bz, medium, particle, alpha1)) {
      std::cerr << m_className << "::IntegrateAlpha:\n"
                << "    Cannot retrieve Townsend coefficient at end point.\n";
      break;
    }
    // Calculate the Townsend coefficient at the mid point of the step.
    Vec xm;
    for (unsigned int i = 0; i < 3; ++i) xm[i] = 0.5 * (x0[i] + x1[i]);    
    if (GetField(xm, ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::IntegrateAlpha: Invalid mid point.\n";
      break;
    }
    double alpham = 0.;
    if (!GetAlpha(ex, ey, ez, bx, by, bz, medium, particle, alpham)) {
      std::cerr << m_className << "::IntegrateAlpha:\n"
                << "    Cannot retrieve Townsend coefficient at mid point.\n";
      break;
    }
    // Compare first and second order estimates.
    if (d * fabs(alpha0 - 2 * alpham + alpha1) / 3. < tol) {
      // Accuracy is good enough.
      integral += d * (alpha0 + 4 * alpham + alpha1) / 6.;
      // Proceed to the next step.
      x0 = x1;
      x1 = xe;
      alpha0 = alpha1;
    } else {
      // Accuracy is not good enough, so halve the step.
      x1 = xm;
      alpha1 = alpham;
    }
  }
  return integral;
}

double DriftLineRKF::IntegrateEta(const std::array<double, 3>& xi, 
                                  const std::array<double, 3>& xe,
                                  const Particle particle, const double tol) {

  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;

  // Make sure the starting point is valid.
  Vec x0 = xi;
  if (GetField(x0, ex, ey, ez, bx, by, bz, medium) != 0) {
    std::cerr << m_className << "::IntegrateEta: Invalid starting point "
              << PrintVec(xi) << ".\n";
    return 0.;
  }
  // Determine the attachment coefficient at the initial point.
  double eta0 = 0.;
  if (!GetEta(ex, ey, ez, bx, by, bz, medium, particle, eta0)) {
    std::cerr << m_className << "::IntegrateEta:\n"
              << "    Cannot retrieve att. coefficient at initial point.\n";
    return 0.;
  }
  // Make sure the end point is valid.
  Vec x1 = xe;
  if (GetField(x1, ex, ey, ez, bx, by, bz, medium) != 0) {
    std::cerr << m_className << "::IntegrateEta: Invalid end point "
              << PrintVec(xe) << ".\n";
    return 0.;
  }
  // Determine the attachment coefficient at the end point.
  double eta1 = 0.;
  if (!GetEta(ex, ey, ez, bx, by, bz, medium, particle, eta1)) {
    std::cerr << m_className << "::IntegrateEta:\n"
              << "    Cannot retrieve att. coefficient at end point.\n";
    return 0.;
  }
  double integral = 0.;
  bool done = false;
  while (!done) {
    // Check if we are close to the end point.
    if (Mag(xe[0] - x0[0], xe[1] - x0[1], xe[2] - x0[2]) < 1.e-6) {
      done = true;
      break;
    }
    const double d = Mag(x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]);
    if (d < 1.e-6) {
      // Step length has become very small.
      if (m_debug) {
        std::cout << m_className << "::IntegrateEta: Small step.\n";
      }
      integral += eta0 * d;
      // Proceed with the next step.
      x0 = x1;
      x1 = xe;
      continue;
    }
    // Calculate the attachment coefficient at the end point of the step.
    if (GetField(x1, ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::IntegrateEta: Invalid end point.\n";
      break;
    }
    if (!GetEta(ex, ey, ez, bx, by, bz, medium, particle, eta1)) {
      std::cerr << m_className << "::IntegrateEta:\n"
                << "    Cannot retrieve att. coefficient at end point.\n";
      break;
    }
    // Calculate the attachment coefficient at the mid point of the step.
    Vec xm;
    for (unsigned int i = 0; i < 3; ++i) xm[i] = 0.5 * (x0[i] + x1[i]);    
    if (GetField(xm, ex, ey, ez, bx, by, bz, medium) != 0) {
      std::cerr << m_className << "::IntegrateEta: Invalid mid point.\n";
      break;
    }
    double etam = 0.;
    if (!GetEta(ex, ey, ez, bx, by, bz, medium, particle, etam)) {
      std::cerr << m_className << "::IntegrateEta:\n"
                << "    Cannot retrieve att. coefficient at mid point.\n";
      break;
    }
    // Compare first and second order estimates.
    if (d * fabs(eta0 - 2 * etam + eta1) / 3. < tol) {
      // Accuracy is good enough.
      integral += d * (eta0 + 4 * etam + eta1) / 6.;
      // Proceed to the next step.
      x0 = x1;
      x1 = xe;
      eta0 = eta1;
    } else {
      // Accuracy is not good enough, so halve the step.
      x1 = xm;
      eta1 = etam;
    }
  }
  return integral;
}

void DriftLineRKF::ComputeSignal(const Particle particle, const double scale,
                                 const std::vector<double>& ts,
                                 const std::vector<Vec>& xs,
                                 const std::vector<double>& ne) const {

  const unsigned int nPoints = ts.size();
  if (nPoints < 2) return;
  const double q = particle == Particle::Electron ? -1 * scale : scale;

  // Get the drift velocity at each point.
  std::vector<std::array<double, 3> > vs;
  for (const auto& x : xs) {
    int stat = 0;
    Vec v;
    if (!GetVelocity(x, particle, v, stat)) {
      std::cerr << m_className << "::ComputeSignal:\n"
                << "    Cannot retrieve velocity at " << PrintVec(x) << "\n";
    }
    vs.push_back(std::move(v));
  }
  m_sensor->AddSignal(q, ts, xs, vs, ne, m_navg);
}

bool DriftLineRKF::FieldLine(const double xi, const double yi, const double zi,
                             std::vector<std::array<float, 3> >& xl, 
                             const bool electron) {

  xl.clear();

  // Check if the sensor is defined.
  if (!m_sensor) {
    std::cerr << m_className << "::FieldLine: Sensor is not defined.\n";
    return false;
  }

  // Get the sensor's bounding box.
  double xmin = 0., xmax = 0.;
  double ymin = 0., ymax = 0.;
  double zmin = 0., zmax = 0.;
  bool bbox = m_sensor->GetArea(xmin, ymin, zmin, xmax, ymax, zmax);

  // Make sure the initial position is at a valid location.
  double ex = 0., ey = 0., ez = 0.;
  Medium* medium = nullptr;
  int stat = 0;
  m_sensor->ElectricField(xi, yi, zi, ex, ey, ez, medium, stat);
  if (!medium || stat != 0) {
    std::cerr << m_className << "::FieldLine:\n"
              << "    No valid field at initial position.\n";
    return false;
  }
  Vec x0 = {xi, yi, zi};
  Vec f0 = {ex, ey, ez};
  if (electron) for (auto& f : f0) f *= -1; 

  // Set the numerical constants for the RKF integration.
  constexpr double c10 = 214. / 891.;
  constexpr double c11 = 1. / 33.;
  constexpr double c12 = 650. / 891.;
  constexpr double c20 = 533. / 2106.;
  constexpr double c22 = 800. / 1053.;
  constexpr double c23 = -1. / 78.;

  constexpr double b10 = 1. / 4.;
  constexpr double b20 = -189. / 800.;
  constexpr double b21 = 729. / 800.;
  constexpr double b30 = 214. / 891.;
  constexpr double b31 = 1. / 33.;
  constexpr double b32 = 650. / 891.;

  const double fmag0 = Mag(f0);
  if (fmag0 < Small) {
    std::cerr << m_className << "::FieldLine:\n"
              << "    Zero field at initial position.\n";
    return false;
  }

  // Initialise time step and previous time step.
  double h = m_accuracy / fmag0;
  double hprev = h;

  // Set the initial point.
  xl.push_back({float(x0[0]), float(x0[1]), float(x0[2])});

  int initCycle = 3;
  bool ok = true;
  while (ok) {
    // Get the field at the first probe point.
    Vec x1 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x1[i] += h * b10 * f0[i];
    }
    m_sensor->ElectricField(x1[0], x1[1], x1[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      if (m_debug) {
        std::cout << m_className << "::FieldLine: Point 1 outside.\n";
      }
      Terminate(x0, x1, xl);
      return true;
    }
    Vec f1 = {ex, ey, ez};
    if (electron) for (auto& f : f1) f *= -1;
    // Get the field at the second probe point.
    Vec x2 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x2[i] += h * (b20 * f0[i] + b21 * f1[i]);
    }
    m_sensor->ElectricField(x2[0], x2[1], x2[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      if (m_debug) {
        std::cout << m_className << "::FieldLine: Point 2 outside.\n";
      }
      Terminate(x0, x2, xl);
      return true;
    }
    Vec f2 = {ex, ey, ez};
    if (electron) for (auto& f : f2) f *= -1;
    // Get the field at the third probe point.
    Vec x3 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x3[i] += h * (b30 * f0[i] + b31 * f1[i] + b32 * f2[i]);
    }
    m_sensor->ElectricField(x3[0], x3[1], x3[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      if (m_debug) {
        std::cout << m_className << "::FieldLine: Point 3 outside.\n";
      }
      Terminate(x0, x3, xl);
      return true;
    }
    Vec f3 = {ex, ey, ez};
    if (electron) for (auto& f : f3) f *= -1;
    // Check if we crossed a wire.
    double xw = 0., yw = 0., zw = 0., rw = 0.;
    if (m_sensor->IsWireCrossed(x0[0], x0[1], x0[2], 
                                x1[0], x1[1], x1[2], xw, yw, zw, false, rw) ||
        m_sensor->IsWireCrossed(x0[0], x0[1], x0[2], 
                                x2[0], x2[1], x2[2], xw, yw, zw, false, rw) ||
        m_sensor->IsWireCrossed(x0[0], x0[1], x0[2], 
                                x3[0], x3[1], x3[2], xw, yw, zw, false, rw)) {
      // TODO!
      xl.push_back({float(xw), float(yw), float(zw)});
      return true;
      // Drift to wire.
      if (h > Small) {
        h *= 0.5;
        continue;
      } else {
        std::cerr << m_className << "::FieldLine: Step size too small. Stop.\n";
        return false;
      }
    }
    // Calculate the correction terms.
    Vec phi1 = {0., 0., 0.};
    Vec phi2 = {0., 0., 0.};
    for (unsigned int i = 0; i < 3; ++i) {
      phi1[i] = c10 * f0[i] + c11 * f1[i] + c12 * f2[i];
      phi2[i] = c20 * f0[i] + c22 * f2[i] + c23 * f3[i];
    }
    // Check if the step length is valid.
    const double phi1mag = Mag(phi1);
    if (phi1mag < Small) {
      std::cerr << m_className << "::FieldLine: Step has zero length. Stop.\n";
      break;
    } else if (m_useStepSizeLimit && h * phi1mag > m_maxStepSize) {
      if (m_debug) {
        std::cout << m_className << "::FieldLine: Step is considered too long. "
                  << "H is reduced.\n";
      }
      h = 0.5 * m_maxStepSize / phi1mag;
      continue;
    } else if (bbox) {
      // Don't allow h to become too large.
      if (h * fabs(phi1[0]) > 0.1 * fabs(xmax - xmin) ||
          h * fabs(phi1[1]) > 0.1 * fabs(ymax - ymin)) {
        h *= 0.5;
        if (m_debug) {
          std::cout << m_className << "::FieldLine: Step is considered too long. "
                    << "H is divided by two.\n";
        }
        continue;
      }
    }
    if (m_debug) std::cout << m_className << "::FieldLine: Step size ok.\n";
    // Update the position.
    for (unsigned int i = 0; i < 3; ++i) x0[i] += h * phi1[i];
    // Check the new position.
    m_sensor->ElectricField(x0[0], x0[1], x0[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      // The new position is not inside a valid drift medium.
      // Terminate the drift line.
      if (m_debug) {
        std::cout << m_className << "::FieldLine: Point outside. Terminate.\n";
      }
      std::array<double, 3> xp = {xl.back()[0], xl.back()[1], xl.back()[2]};
      Terminate(xp, x0, xl);
      return true;
    }
    // Add the new point to the drift line.
    xl.push_back({float(x0[0]), float(x0[1]), float(x0[2])});
    // Adjust the step size according to the accuracy of the two estimates.
    hprev = h;
    const double dphi = fabs(phi1[0] - phi2[0]) + fabs(phi1[1] - phi2[1]) +
                        fabs(phi1[2] - phi2[2]);
    if (dphi > 0) {
      h = sqrt(h * m_accuracy / dphi);
      if (m_debug) {
        std::cout << m_className << "::FieldLine: Adapting H to " << h << ".\n";
      }
    } else {
      h *= 2;
      if (m_debug) {
        std::cout << m_className << "::FieldLine: H increased by factor two.\n";
      }
    }
    // Make sure that H is different from zero; this should always be ok.
    if (h < Small) {
      std::cerr << m_className << "::FieldLine: Step size is zero. Stop.\n";
      return false;
    }
    // Check the initial step size.
    if (initCycle > 0 && h < 0.2 * hprev) {
      if (m_debug) {
        std::cout << m_className << "::FieldLine: Reinitialise step size.\n";
      }
      --initCycle;
      x0 = {xi, yi, zi};
      xl.clear();
      xl.push_back({float(xi), float(yi), float(zi)});
      continue;
    }
    initCycle = 0;
    // Don't allow H to grow too quickly
    if (h > 10 * hprev) {
      h = 10 * hprev;
      if (m_debug) {
        std::cout << m_className << "::FieldLine: H restricted to 10 times "
                  << "the previous value.\n";
      }
    }
    // Stop in case H tends to become too small.
    if (h * (fabs(phi1[0]) + fabs(phi1[1]) + fabs(phi1[2])) < m_accuracy) {
      std::cerr << m_className << "::FieldLine: Step size has become smaller "
                << "than int. accuracy. Stop.\n";
      return false;
    }
    // Update the field.
    f0 = f3;
  }
  return true;
}

void DriftLineRKF::Terminate(const std::array<double, 3>& xx0,
                             const std::array<double, 3>& xx1,
                             std::vector<std::array<float, 3> >& xs) {

  // Final point just inside the medium.
  Vec x0 = xx0;
  // Final point just outside the medium.
  Vec x1 = xx1;
  // Perform some bisections.
  const unsigned int nBisections = 20;
  for (unsigned int i = 0; i < nBisections; ++i) {
    // Quit bisection when interval becomes too small.
    bool small = true;
    for (unsigned int j = 0; j < 3; ++j) {
      if (fabs(x1[j] - x0[j]) > 1.e-6 * (fabs(x0[j]) + fabs(x1[j]))) {
        small = false;
        break;
      }
    } 
    if (small) {
      if (m_debug) {
        std::cout << m_className << "::Terminate:\n"
                  << "    Bisection ended at cycle " << i << ".\n";
      }
      break; 
    }
    // Calculate the mid point and evaluate the field.
    Vec xm;
    for (unsigned int j = 0; j < 3; ++j) xm[j] = 0.5 * (x0[j] + x1[j]);
    double ex = 0., ey = 0., ez = 0.;
    Medium* medium = nullptr;
    int status = 0;
    m_sensor->ElectricField(xm[0], xm[1], xm[2], ex, ey, ez, medium, status);
    if (status == 0) {
      x0 = xm;
    } else {
      x1 = xm;
    }
  }

  xs.push_back({float(x0[0]), float(x0[1]), float(x0[2])});
}
}
