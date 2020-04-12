#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <algorithm>

#include <TAxis.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>

#include "Garfield/ComponentBase.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"

namespace {

void SampleRange(const double xmin, const double ymin, const double xmax,
                 const double ymax, TF2* f, double& zmin, double& zmax) {
  constexpr unsigned int n = 1000;
  const double dx = xmax - xmin;
  const double dy = ymax - ymin;
  zmin = std::numeric_limits<double>::max();
  zmax = -zmin;
  for (unsigned int i = 0; i < n; ++i) {
    const double z = f->Eval(xmin + Garfield::RndmUniform() * dx,
                             ymin + Garfield::RndmUniform() * dy);
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }
}

void SampleRange(TF1* f, double& ymin, double& ymax) {
  constexpr unsigned int n = 1000;
  ymin = std::numeric_limits<double>::max();
  ymax = -ymin;
  for (unsigned int i = 0; i < n; ++i) {
    const double y = f->Eval(Garfield::RndmUniform());
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
  }
}

double Interpolate(const std::array<double, 1000>& y,
                   const std::array<double, 1000>& x, const double xx) {

  const double tol = 1.e-6 * fabs(x.back() - x.front());
  if (xx < x[0]) return y[0];
  const auto it1 = std::upper_bound(x.cbegin(), x.cend(), xx);
  if (it1 == x.cend()) return y.back();
  const auto it0 = std::prev(it1);
  const double dx = (*it1 - *it0);
  if (dx < tol) return y[it0 - x.cbegin()];
  const double f = (xx - *it0) / dx;
  return y[it0 - x.cbegin()] * (1. - f) + f * y[it1 - x.cbegin()];
}

}

namespace Garfield {

ViewField::ViewField() : ViewBase("ViewField") { }

void ViewField::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }
  m_sensor = s;
  m_component = nullptr;
}

void ViewField::SetComponent(ComponentBase* c) {
  if (!c) {
    std::cerr << m_className << "::SetComponent: Null pointer.\n";
    return;
  }
  m_component = c;
  m_sensor = nullptr;
}

void ViewField::SetVoltageRange(const double vmin, const double vmax) {
  m_vmin = std::min(vmin, vmax);
  m_vmax = std::max(vmin, vmax);
  m_useAutoRange = false;
}

void ViewField::SetElectricFieldRange(const double emin, const double emax) {
  m_emin = std::min(emin, emax);
  m_emax = std::max(emin, emax);
  m_useAutoRange = false;
}

void ViewField::SetWeightingFieldRange(const double wmin, const double wmax) {
  m_wmin = std::min(wmin, wmax);
  m_wmax = std::max(wmin, wmax);
  m_useAutoRange = false;
}

void ViewField::SetNumberOfContours(const unsigned int n) {
  if (n > 0) m_nContours = n;
}

void ViewField::SetNumberOfSamples1d(const unsigned int n) {
  m_nSamples1d = std::max(4u, n);
}

void ViewField::SetNumberOfSamples2d(const unsigned int nx,
                                     const unsigned int ny) {
  m_nSamples2dX = std::max(4u, nx);
  m_nSamples2dY = std::max(4u, ny);
}

void ViewField::PlotContour(const std::string& option) {
  Draw2d(option, true, false, "", "CONT4Z");
}

void ViewField::Plot(const std::string& option, const std::string& drawopt) {
  Draw2d(option, false, false, "", drawopt);
}

void ViewField::PlotProfile(const double x0, const double y0, const double z0,
                            const double x1, const double y1, const double z1,
                            const std::string& option) {
  DrawProfile(x0, y0, z0, x1, y1, z1, option, false, "");
}

void ViewField::PlotWeightingField(const std::string& label,
                                   const std::string& option,
                                   const std::string& drawopt) {
  Draw2d(option, false, true, label, drawopt);
}

void ViewField::PlotContourWeightingField(const std::string& label,
                                          const std::string& option) {
  Draw2d(option, true, true, label, "CONT4Z");
}

void ViewField::PlotProfileWeightingField(const std::string& label,
                                          const double x0, const double y0,
                                          const double z0, const double x1,
                                          const double y1, const double z1,
                                          const std::string& option) {
  DrawProfile(x0, y0, z0, x1, y1, z1, option, true, label);
}


ViewField::Parameter ViewField::GetPar(const std::string& option,
                                       std::string& title) const {

  std::string opt;
  std::transform(option.begin(), option.end(), 
                 std::back_inserter(opt), toupper);
  if (opt == "V" || opt == "P" || opt == "PHI" || 
      opt.find("VOLT") != std::string::npos ||
      opt.find("POT") != std::string::npos) {
    title = "potential";
    return Parameter::Potential;
  } else if (opt == "E" || opt == "FIELD" || opt == "NORM" ||
             opt.find("MAG") != std::string::npos) {
    title = "field";
    return Parameter::Magnitude;
  } else if (opt.find("X") != std::string::npos) {
    title = "field (x-component)";
    return Parameter::Ex;
  } else if (opt.find("Y") != std::string::npos) {
    title = "field (y-component)";
    return Parameter::Ey;
  } else if (opt.find("Z") != std::string::npos) {
    title = "field (z-component)";
    return Parameter::Ez;
  }
  std::cerr << m_className << "::GetPar: Unknown option (" << option << ").\n";
  title = "potential";
  return Parameter::Potential;
}

void ViewField::Draw2d(const std::string& option, const bool contour, 
                       const bool wfield, const std::string& electrode,
                       const std::string& drawopt) {
  if (!m_sensor && !m_component) {
    std::cerr << m_className << "::Draw2d:\n"
              << "    Neither sensor nor component are defined.\n";
    return;
  }

  // Determine the x-y range.
  if (!SetPlotLimits()) return;

  // Determine the quantity to be plotted.
  std::string title;
  const Parameter par = GetPar(option, title);

  auto eval = [this, par, wfield, electrode](double* u, double* /*p*/) {
    // Transform to global coordinates.
    const double x = m_proj[0][0] * u[0] + m_proj[1][0] * u[1] + m_proj[2][0];
    const double y = m_proj[0][1] * u[0] + m_proj[1][1] * u[1] + m_proj[2][1];
    const double z = m_proj[0][2] * u[0] + m_proj[1][2] * u[1] + m_proj[2][2];
    return wfield ? Wfield(x, y, z, par, electrode) : Field(x, y, z, par);
  };
  const std::string fname = FindUnusedFunctionName("f2D");
  TF2 f2(fname.c_str(), eval, m_xMinPlot, m_xMaxPlot, m_yMinPlot, m_yMaxPlot, 0);

  // Set the x-y range.
  f2.SetRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot);

  // Set the z-range.
  double zmin = m_vmin;
  double zmax = m_vmax;
  if (wfield) {
    if (contour) {
      title = "Contours of the weighting " + title;
    } else {
      title = "Weighting " + title;
    }
    if (m_useAutoRange) {
      SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, &f2, 
                  zmin, zmax);
    } else if (par == Parameter::Potential) {
      zmin = 0.;
      zmax = 1.;
    } else {
      zmin = m_wmin;
      zmax = m_wmax;
    }
  } else {
    if (contour) {
      title = "Contours of the electric " + title;
    } else {
      title = "Electric " + title;
    }
    if (par == Parameter::Potential) {
      if (m_useAutoRange) {
        if (m_component) {
          if (!m_component->GetVoltageRange(zmin, zmax)) {
            SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, 
                        &f2, zmin, zmax);
          }
        } else if (m_sensor) {
          if (!m_sensor->GetVoltageRange(zmin, zmax)) {
            SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, 
                        &f2, zmin, zmax);
          }
        }
      } else {
        zmin = m_vmin;
        zmax = m_vmax;
      }
    } else {
      if (m_useAutoRange) {
        SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, 
                    &f2, zmin, zmax);
      } else {
        zmin = m_emin;
        zmax = m_emax;
      }
    }
  }
  f2.SetMinimum(zmin);
  f2.SetMaximum(zmax);

  // Set the contours if requested.
  if (contour) {
    std::vector<double> level(m_nContours, 0.);
    if (m_nContours > 1) {
      const double step = (zmax - zmin) / (m_nContours - 1.);
      for (unsigned int i = 0; i < m_nContours; ++i) {
        level[i] = zmin + i * step;
      }
    } else {
      level[0] = 0.5 * (zmax + zmin);
    }
    if (m_debug) {
      std::cout << m_className << "::Draw2d:\n"
                << "    Number of contours: " << m_nContours << "\n";
      for (unsigned int i = 0; i < m_nContours; ++i) {
        std::cout << "        Level " << i << " = " << level[i] << "\n";
      }
    }
    f2.SetContour(m_nContours, level.data());
  }

  // Set the resolution.
  f2.SetNpx(m_nSamples2dX);
  f2.SetNpy(m_nSamples2dY);

  // Set the labels.
  std::string labels = ";" + LabelX() + ";" + LabelY();
  f2.SetTitle(labels.c_str());

  auto canvas = GetCanvas();
  canvas->cd();
  canvas->SetTitle(title.c_str());
  f2.DrawCopy(drawopt.c_str());
  gPad->SetRightMargin(0.15);
  gPad->Update();
}

void ViewField::DrawProfile(const double x0, const double y0, const double z0,
                            const double x1, const double y1, const double z1,
                            const std::string& option, 
                            const bool wfield, const std::string& electrode) {
  if (!m_sensor && !m_component) {
    std::cerr << m_className << "::DrawProfile:\n"
              << "    Neither sensor nor component are defined.\n";
    return;
  }

  // Check the distance between the two points.
  const double dx = x1 - x0;
  const double dy = y1 - y0;
  const double dz = z1 - z0;
  if (dx * dx + dy * dy + dz * dz <= 0.) {
    std::cerr << m_className << "::DrawProfile:\n"
              << "    Start and end points coincide.\n";
    return;
  }

  // Determine the quantity to be plotted.
  std::string title;
  const Parameter par = GetPar(option, title);

  auto eval = [this, par, wfield, electrode, 
               x0, y0, z0, dx, dy, dz](double* u, double* /*p*/) {
    // Get the position.
    const double t = u[0];
    const double x = x0 + t * dx;
    const double y = y0 + t * dy;
    const double z = z0 + t * dz;
    return wfield ? Wfield(x, y, z, par, electrode) : Field(x, y, z, par);
  };

  const std::string fname = FindUnusedFunctionName("fProfile");
  TF1 f1(fname.c_str(), eval, 0., 1., 0);

  double fmin = m_vmin;
  double fmax = m_vmax;
  if (wfield) {
    title = "weighting " + title;
    if (par == Parameter::Potential) {
      fmin = 0.;
      fmax = 1.;
    } else {
      if (m_useAutoRange) {
        SampleRange(&f1, fmin, fmax);
      } else {
        fmin = m_wmin;
        fmax = m_wmax;
      }
    }
  } else {
    title = "electric " + title;
    if (par == Parameter::Potential) {
      if (m_useAutoRange) {
        if (m_component) {
          if (!m_component->GetVoltageRange(fmin, fmax)) {
            SampleRange(&f1, fmin, fmax);
          }
        } else if (m_sensor) {
          if (!m_sensor->GetVoltageRange(fmin, fmax)) {
            SampleRange(&f1, fmin, fmax);
          }
        }
      } else {
        fmin = m_vmin;
        fmax = m_vmax;
      }
    } else {
      if (m_useAutoRange) {
        SampleRange(&f1, fmin, fmax);
      } else {
        fmin = m_emin;
        fmax = m_emax;
      }
    }
  }
  f1.SetMinimum(fmin);
  f1.SetMaximum(fmax);

  std::string labels = ";normalised distance;";
  if (par == Parameter::Potential) {
    labels += "#phi";
    if (wfield) {
      labels += "_w";
    } else {
      labels += " [V]";
    }
  } else {
    labels += "#it{E}";
    if (wfield) {
      labels += "_{w";
      if (par != Parameter::Magnitude) labels += ",";
    } else if (par != Parameter::Magnitude) {
      labels += "_{";
    }
    if (par == Parameter::Ex) {
      labels += "x";
    } else if (par == Parameter::Ey) {
      labels += "y";
    } else if (par == Parameter::Ez) {
      labels += "z";
    }
    if (wfield || par != Parameter::Magnitude) labels += "}";
    if (wfield) {
      labels += " [1/cm]";
    } else {
      labels += " [V/cm]";
    }
  }
  f1.SetTitle(labels.c_str());
  f1.SetNpx(m_nSamples1d);

  auto canvas = GetCanvas();
  canvas->cd();
  title = "Profile plot of the " + title;
  canvas->SetTitle(title.c_str());
  f1.DrawCopy();
  gPad->Update();
}

bool ViewField::SetPlotLimits() {

  if (m_userPlotLimits) return true;
  double xmin = 0., ymin = 0., xmax = 0., ymax = 0.;
  if (m_userBox) {
    if (PlotLimitsFromUserBox(xmin, ymin, xmax, ymax)) {
      m_xMinPlot = xmin;
      m_xMaxPlot = xmax;
      m_yMinPlot = ymin;
      m_yMaxPlot = ymax;
      return true;
    } 
  }
  // Try to get the area/bounding box from the sensor/component.
  bool ok = false;

  if (m_sensor) {
    ok = PlotLimits(m_sensor, xmin, ymin, xmax, ymax);
  } else {
    ok = PlotLimits(m_component, xmin, ymin, xmax, ymax);
  } 
  if (ok) {
    m_xMinPlot = xmin;
    m_xMaxPlot = xmax;
    m_yMinPlot = ymin;
    m_yMaxPlot = ymax;
  }
  return ok;
}

double ViewField::Field(const double x, const double y, const double z,
                        const Parameter par) const {

  // Compute the field.
  double ex = 0., ey = 0., ez = 0., volt = 0.;
  int status = 0;
  Medium* medium = nullptr;
  if (!m_sensor) {
    m_component->ElectricField(x, y, z, ex, ey, ez, volt, medium, status);
  } else {
    m_sensor->ElectricField(x, y, z, ex, ey, ez, volt, medium, status);
  }
  if (m_useStatus && status != 0) return m_vBkg;
  switch (par) {
    case Parameter::Potential:
      return volt;
      break;
    case Parameter::Magnitude:
      return sqrt(ex * ex + ey * ey + ez * ez);
      break;
    case Parameter::Ex:
      return ex;
      break;
    case Parameter::Ey:
      return ey;
      break;
    case Parameter::Ez:
      return ez;
      break;
    default:
      break;
  }
  return volt;
}

double ViewField::Wfield(const double x, const double y, const double z,
                         const Parameter par, 
                         const std::string& electrode) const {

  if (par == Parameter::Potential) {
    if (m_sensor) {
      return m_sensor->WeightingPotential(x, y, z, electrode);
    } else {
      return m_component->WeightingPotential(x, y, z, electrode);
    }
  }

  double ex = 0., ey = 0., ez = 0.;
  if (!m_sensor) {
    m_component->WeightingField(x, y, z, ex, ey, ez, electrode);
  } else {
    m_sensor->WeightingField(x, y, z, ex, ey, ez, electrode);
  }

  switch (par) {
    case Parameter::Magnitude:
      return sqrt(ex * ex + ey * ey + ez * ez);
      break;
    case Parameter::Ex:
      return ex;
      break;
    case Parameter::Ey:
      return ey;
      break;
    case Parameter::Ez:
      return ez;
      break;
    default:
      break;
  }
  return 0.;
} 

bool ViewField::EqualFluxIntervals(
    const double x0, const double y0, const double z0,
    const double x1, const double y1, const double z1,
    std::vector<double>& xf, std::vector<double>& yf,
    std::vector<double>& zf,
    const unsigned int nPoints) const {

  if (nPoints < 2) {
    std::cerr << m_className << "::EqualFluxIntervals:\n"
              << "    Number of flux lines must be > 1.\n";
    return false;
  }

  // Set integration intervals.
  constexpr unsigned int nV = 5;
  // Compute the inplane vector normal to the track.
  const double xp = (y1 - y0) * m_plane[2] - (z1 - z0) * m_plane[1];
  const double yp = (z1 - z0) * m_plane[0] - (x1 - x0) * m_plane[2];
  const double zp = (x1 - x0) * m_plane[1] - (y1 - y0) * m_plane[0];
  // Compute the total flux, accepting positive and negative parts.
  double q = 0.;
  if (m_component) {
    q = m_component->IntegrateFluxLine(x0, y0, z0, x1, y1, z1, 
                                       xp, yp, zp, 20 * nV, 0);
  } else {
    q = m_sensor->IntegrateFluxLine(x0, y0, z0, x1, y1, z1, 
                                    xp, yp, zp, 20 * nV, 0);
  }
  const int isign = q > 0 ? +1 : -1;
  if (m_debug) {
    std::cout << m_className << "::EqualFluxIntervals:\n";
    std::printf("    Total flux: %15.e8\n", q);
  }
  // Compute the 1-sided flux in a number of steps.
  double fsum = 0.;
  unsigned int nOtherSign = 0;
  double s0 = -1.;
  double s1 = -1.;
  constexpr size_t nSteps = 1000;
  std::array<double, nSteps> sTab;
  std::array<double, nSteps> fTab;
  constexpr double ds = 1. / nSteps;
  const double dx = (x1 - x0) * ds;
  const double dy = (y1 - y0) * ds;
  const double dz = (z1 - z0) * ds;
  for (size_t i = 0; i < nSteps; ++i) {
    const double x = x0 + i * dx;
    const double y = y0 + i * dy;
    const double z = z0 + i * dz;
    if (m_component) {
      q = m_component->IntegrateFluxLine(x, y, z, x + dx, y + dy, z + dz, 
                                         xp, yp, zp, nV, isign);
    } else {
      q = m_sensor->IntegrateFluxLine(x, y, z, x + dx, y + dy, z + dz, 
                                      xp, yp, zp, nV, isign);
    }
    sTab[i] = (i + 1) * ds;
    if (q > 0) {
      fsum += q;
      if (s0 < -0.5) s0 = i * ds;
      s1 = (i + 1) * ds;
    }
    if (q < 0) ++nOtherSign;
    fTab[i] = fsum;
  }
  if (m_debug) {
    std::printf("    Used flux: %15.8e V. Start: %10.3f End: %10.3f\n",
                fsum, s0, s1);
  }
  // Make sure that the sum is positive.
  if (fsum <= 0) {
    std::cerr << m_className << "::EqualFluxIntervals:\n"
              << "    1-Sided flux integral is not > 0.\n";
    return false;
  } else if (s0 < -0.5 || s1 < -0.5 || s1 <= s0) {
    std::cerr << m_className << "::EqualFluxIntervals:\n"
              << "    No flux interval without sign change found.\n";
    return false;
  } else if (nOtherSign > 0) {
    std::cerr << m_className << "::EqualFluxIntervals:\n"
              << " The flux changes sign over the line.\n";
  }
  // Normalise the flux.
  const double scale = (nPoints - 1) / fsum;
  for (size_t i = 0; i < nSteps; ++i) fTab[i] *= scale;

  // Compute new cluster position.
  for (size_t i = 0; i < nPoints; ++i) {
    double s = std::min(s1, std::max(s0, Interpolate(sTab, fTab, i)));
    xf.push_back(x0 + s * (x1 - x0));
    yf.push_back(y0 + s * (y1 - y0));
    zf.push_back(z0 + s * (z1 - z0));
  }
  return true;
}

bool ViewField::FixedFluxIntervals(
    const double x0, const double y0, const double z0,
    const double x1, const double y1, const double z1,
    std::vector<double>& xf, std::vector<double>& yf,
    std::vector<double>& zf, const double interval) const {

  if (interval <= 0.) {
    std::cerr << m_className << "::FixedFluxIntervals:\n"
              << "    Flux interval must be > 0.\n";
    return false;
  }

  // Set integration intervals.
  constexpr unsigned int nV = 5;
  // Compute the inplane vector normal to the track.
  const double xp = (y1 - y0) * m_plane[2] - (z1 - z0) * m_plane[1];
  const double yp = (z1 - z0) * m_plane[0] - (x1 - x0) * m_plane[2];
  const double zp = (x1 - x0) * m_plane[1] - (y1 - y0) * m_plane[0];
  // Compute the total flux, accepting positive and negative parts.
  double q = 0.;
  if (m_component) {
    q = m_component->IntegrateFluxLine(x0, y0, z0, x1, y1, z1, 
                                       xp, yp, zp, 20 * nV, 0);
  } else {
    q = m_sensor->IntegrateFluxLine(x0, y0, z0, x1, y1, z1, 
                                    xp, yp, zp, 20 * nV, 0);
  }
  const int isign = q > 0 ? +1 : -1;
  if (m_debug) {
    std::cout << m_className << "::FixedFluxIntervals:\n";
    std::printf("    Total flux: %15.8e V\n", q);
  }
  // Compute the 1-sided flux in a number of steps.
  double fsum = 0.;
  unsigned int nOtherSign = 0;
  double s0 = -1.;
  constexpr size_t nSteps = 1000;
  std::array<double, nSteps> sTab;
  std::array<double, nSteps> fTab;
  constexpr double ds = 1. / nSteps;
  const double dx = (x1 - x0) * ds;
  const double dy = (y1 - y0) * ds;
  const double dz = (z1 - z0) * ds;
  for (size_t i = 0; i < nSteps; ++i) {
    const double x = x0 + i * dx;
    const double y = y0 + i * dy;
    const double z = z0 + i * dz;
    if (m_component) {
      q = m_component->IntegrateFluxLine(x, y, z, x + dx, y + dy, z + dz, 
                                         xp, yp, zp, nV, isign);
    } else {
      q = m_sensor->IntegrateFluxLine(x, y, z, x + dx, y + dy, z + dz, 
                                      xp, yp, zp, nV, isign);
    }
    sTab[i] = (i + 1) * ds;
    if (q > 0) {
      fsum += q;
      if (s0 < -0.5) s0 = i * ds;
    }
    if (q < 0) ++nOtherSign;
    fTab[i] = fsum;
  }
  // Make sure that the sum is positive.
  if (m_debug) {
    std::printf("    Used flux: %15.8e V. Start offset: %10.3f\n", fsum, s0);
  }
  if (fsum <= 0) {
    std::cerr << m_className << "::FixedFluxIntervals:\n"
              << "    1-Sided flux integral is not > 0.\n";
    return false;
  } else if (s0 < -0.5) {
    std::cerr << m_className << "::FixedFluxIntervals:\n"
              << "    No flux interval without sign change found.\n";
    return false;
  } else if (nOtherSign > 0) {
    std::cerr << m_className << "::FixedFluxIntervals:\n"
              << "    Warning: The flux changes sign over the line.\n";
  }

  double f = 0.;
  while (f < fsum) {
    const double s = Interpolate(sTab, fTab, f);
    f += interval;
    xf.push_back(x0 + s * (x1 - x0));
    yf.push_back(y0 + s * (y1 - y0));
    zf.push_back(z0 + s * (z1 - z0));
  }
  return true;
}
}
