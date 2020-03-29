#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <limits>
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
}

namespace Garfield {

ViewField::ViewField() : ViewBase("ViewField") {
  SetDefaultProjection();
}

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

void ViewField::SetArea(const double xmin, const double ymin, const double xmax,
                        const double ymax) {
  // Check range, assign if non-null.
  if (xmin == xmax || ymin == ymax) {
    std::cerr << m_className << "::SetArea: Null area range is not permitted.\n"
              << "      " << xmin << " < x < " << xmax << "\n"
              << "      " << ymin << " < y < " << ymax << "\n";
    return;
  }
  m_xMinPlot = std::min(xmin, xmax);
  m_yMinPlot = std::min(ymin, ymax);
  m_xMaxPlot = std::max(xmin, xmax);
  m_yMaxPlot = std::max(ymin, ymax);
  m_userPlotLimits = true;
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


void ViewField::SetDefaultProjection() {
  // Default projection: x-y at z=0
  m_proj[0][0] = 1;
  m_proj[1][0] = 0;
  m_proj[2][0] = 0;
  m_proj[0][1] = 0;
  m_proj[1][1] = 1;
  m_proj[2][1] = 0;
  m_proj[0][2] = 0;
  m_proj[1][2] = 0;
  m_proj[2][2] = 0;

  // Plane description
  m_plane[0] = 0;
  m_plane[1] = 0;
  m_plane[2] = 1;
  m_plane[3] = 0;

  // Prepare axis labels.
  Labels();
}

void ViewField::Labels() {
  // Initialisation of the x-axis label
  strcpy(m_xLabel, "\0");
  char buf[100];

  const double tol = 1.e-4;
  // x portion
  if (fabs(m_proj[0][0] - 1) < tol) {
    strcat(m_xLabel, "x");
  } else if (fabs(m_proj[0][0] + 1) < tol) {
    strcat(m_xLabel, "-x");
  } else if (m_proj[0][0] > tol) {
    sprintf(buf, "%g x", m_proj[0][0]);
    strcat(m_xLabel, buf);
  } else if (m_proj[0][0] < -tol) {
    sprintf(buf, "%g x", m_proj[0][0]);
    strcat(m_xLabel, buf);
  }

  // y portion
  if (strlen(m_xLabel) > 0) {
    if (m_proj[0][1] < -tol) {
      strcat(m_xLabel, " - ");
    } else if (m_proj[0][1] > tol) {
      strcat(m_xLabel, " + ");
    }
    if (fabs(m_proj[0][1] - 1) < tol || fabs(m_proj[0][1] + 1) < tol) {
      strcat(m_xLabel, "y");
    } else if (fabs(m_proj[0][1]) > tol) {
      sprintf(buf, "%g y", fabs(m_proj[0][1]));
      strcat(m_xLabel, buf);
    }
  } else {
    if (fabs(m_proj[0][1] - 1) < tol) {
      strcat(m_xLabel, "y");
    } else if (fabs(m_proj[0][1] + 1) < tol) {
      strcat(m_xLabel, "-y");
    } else if (m_proj[0][1] > tol) {
      sprintf(buf, "%g y", m_proj[0][1]);
      strcat(m_xLabel, buf);
    } else if (m_proj[0][1] < -tol) {
      sprintf(buf, "%g y", m_proj[0][1]);
      strcat(m_xLabel, buf);
    }
  }

  // z portion
  if (strlen(m_xLabel) > 0) {
    if (m_proj[0][2] < -tol) {
      strcat(m_xLabel, " - ");
    } else if (m_proj[0][2] > tol) {
      strcat(m_xLabel, " + ");
    }
    if (fabs(m_proj[0][2] - 1) < tol || fabs(m_proj[0][2] + 1) < tol) {
      strcat(m_xLabel, "z");
    } else if (fabs(m_proj[0][2]) > tol) {
      sprintf(buf, "%g z", fabs(m_proj[0][2]));
      strcat(m_xLabel, buf);
    }
  } else {
    if (fabs(m_proj[0][2] - 1) < tol) {
      strcat(m_xLabel, "z");
    } else if (fabs(m_proj[0][2] + 1) < tol) {
      strcat(m_xLabel, "-z");
    } else if (m_proj[0][2] > tol) {
      sprintf(buf, "%g z", m_proj[0][2]);
      strcat(m_xLabel, buf);
    } else if (m_proj[0][2] < -tol) {
      sprintf(buf, "%g z", m_proj[0][2]);
      strcat(m_xLabel, buf);
    }
  }

  // Unit
  strcat(m_xLabel, " [cm]");

  // Initialisation of the y-axis label
  strcpy(m_yLabel, "\0");

  // x portion
  if (fabs(m_proj[1][0] - 1) < tol) {
    strcat(m_yLabel, "x");
  } else if (fabs(m_proj[1][0] + 1) < tol) {
    strcat(m_yLabel, "-x");
  } else if (m_proj[1][0] > tol) {
    sprintf(buf, "%g x", m_proj[1][0]);
    strcat(m_yLabel, buf);
  } else if (m_proj[1][0] < -tol) {
    sprintf(buf, "%g x", m_proj[1][0]);
    strcat(m_yLabel, buf);
  }

  // y portion
  if (strlen(m_yLabel) > 0) {
    if (m_proj[1][1] < -tol) {
      strcat(m_yLabel, " - ");
    } else if (m_proj[1][1] > tol) {
      strcat(m_yLabel, " + ");
    }
    if (fabs(m_proj[1][1] - 1) < tol || fabs(m_proj[1][1] + 1) < tol) {
      strcat(m_yLabel, "y");
    } else if (fabs(m_proj[1][1]) > tol) {
      sprintf(buf, "%g y", fabs(m_proj[1][1]));
      strcat(m_yLabel, buf);
    }
  } else {
    if (fabs(m_proj[1][1] - 1) < tol) {
      strcat(m_yLabel, "y");
    } else if (fabs(m_proj[1][1] + 1) < tol) {
      strcat(m_yLabel, "-y");
    } else if (m_proj[1][1] > tol) {
      sprintf(buf, "%g y", m_proj[1][1]);
      strcat(m_yLabel, buf);
    } else if (m_proj[1][1] < -tol) {
      sprintf(buf, "%g y", m_proj[1][1]);
      strcat(m_yLabel, buf);
    }
  }

  // z portion
  if (strlen(m_yLabel) > 0) {
    if (m_proj[1][2] < -tol) {
      strcat(m_yLabel, " - ");
    } else if (m_proj[1][2] > tol) {
      strcat(m_yLabel, " + ");
    }
    if (fabs(m_proj[1][2] - 1) < tol || fabs(m_proj[1][2] + 1) < tol) {
      strcat(m_yLabel, "z");
    } else if (fabs(m_proj[1][2]) > tol) {
      sprintf(buf, "%g z", fabs(m_proj[1][2]));
      strcat(m_yLabel, buf);
    }
  } else {
    if (fabs(m_proj[1][2] - 1) < tol) {
      strcat(m_yLabel, "z");
    } else if (fabs(m_proj[1][2] + 1) < tol) {
      strcat(m_yLabel, "-z");
    } else if (m_proj[1][2] > tol) {
      sprintf(buf, "%g z", m_proj[1][2]);
      strcat(m_yLabel, buf);
    } else if (m_proj[1][2] < -tol) {
      sprintf(buf, "%g z", m_proj[1][2]);
      strcat(m_yLabel, buf);
    }
  }

  // Unit
  strcat(m_yLabel, " [cm]");

  // Initialisation of the plane label
  strcpy(m_description, "\0");

  // x portion
  if (fabs(m_plane[0] - 1) < tol) {
    strcat(m_description, "x");
  } else if (fabs(m_plane[0] + 1) < tol) {
    strcat(m_description, "-x");
  } else if (m_plane[0] > tol) {
    sprintf(buf, "%g x", m_plane[0]);
    strcat(m_description, buf);
  } else if (m_plane[0] < -tol) {
    sprintf(buf, "%g x", m_plane[0]);
    strcat(m_description, buf);
  }

  // y portion
  if (strlen(m_description) > 0) {
    if (m_plane[1] < -tol) {
      strcat(m_description, " - ");
    } else if (m_plane[1] > tol) {
      strcat(m_description, " + ");
    }
    if (fabs(m_plane[1] - 1) < tol || fabs(m_plane[1] + 1) < tol) {
      strcat(m_description, "y");
    } else if (fabs(m_plane[1]) > tol) {
      sprintf(buf, "%g y", fabs(m_plane[1]));
      strcat(m_description, buf);
    }
  } else {
    if (fabs(m_plane[1] - 1) < tol) {
      strcat(m_description, "y");
    } else if (fabs(m_plane[1] + 1) < tol) {
      strcat(m_description, "-y");
    } else if (m_plane[1] > tol) {
      sprintf(buf, "%g y", m_plane[1]);
      strcat(m_description, buf);
    } else if (m_plane[1] < -tol) {
      sprintf(buf, "%g y", m_plane[1]);
      strcat(m_description, buf);
    }
  }

  // z portion
  if (strlen(m_description) > 0) {
    if (m_plane[2] < -tol) {
      strcat(m_description, " - ");
    } else if (m_plane[2] > tol) {
      strcat(m_description, " + ");
    }
    if (fabs(m_plane[2] - 1) < tol || fabs(m_plane[2] + 1) < tol) {
      strcat(m_description, "z");
    } else if (fabs(m_plane[2]) > tol) {
      sprintf(buf, "%g z", fabs(m_plane[2]));
      strcat(m_description, buf);
    }
  } else {
    if (fabs(m_plane[2] - 1) < tol) {
      strcat(m_description, "z");
    } else if (fabs(m_plane[2] + 1) < tol) {
      strcat(m_description, "-z");
    } else if (m_plane[2] > tol) {
      sprintf(buf, "%g z", m_plane[2]);
      strcat(m_description, buf);
    } else if (m_plane[2] < -tol) {
      sprintf(buf, "%g z", m_plane[2]);
      strcat(m_description, buf);
    }
  }

  // Constant
  sprintf(buf, " = %g", m_plane[3]);
  strcat(m_description, buf);

  if (m_debug) {
    std::cout << m_className << "::Labels:\n"
              << "    x label: |" << m_xLabel << "|\n"
              << "    y label: |" << m_yLabel << "|\n"
              << "    plane:   |" << m_description << "|\n";
  }
}

void ViewField::SetPlane(const double fx, const double fy, const double fz,
                         const double x0, const double y0, const double z0) {
  // Calculate two in-plane vectors for the normal vector
  const double fnorm = sqrt(fx * fx + fy * fy + fz * fz);
  if (fnorm > 0 && fx * fx + fz * fz > 0) {
    const double fxz = sqrt(fx * fx + fz * fz);
    m_proj[0][0] = fz / fxz;
    m_proj[0][1] = 0;
    m_proj[0][2] = -fx / fxz;
    m_proj[1][0] = -fx * fy / (fxz * fnorm);
    m_proj[1][1] = (fx * fx + fz * fz) / (fxz * fnorm);
    m_proj[1][2] = -fy * fz / (fxz * fnorm);
    m_proj[2][0] = x0;
    m_proj[2][1] = y0;
    m_proj[2][2] = z0;
  } else if (fnorm > 0 && fy * fy + fz * fz > 0) {
    const double fyz = sqrt(fy * fy + fz * fz);
    m_proj[0][0] = (fy * fy + fz * fz) / (fyz * fnorm);
    m_proj[0][1] = -fx * fz / (fyz * fnorm);
    m_proj[0][2] = -fy * fz / (fyz * fnorm);
    m_proj[1][0] = 0;
    m_proj[1][1] = fz / fyz;
    m_proj[1][2] = -fy / fyz;
    m_proj[2][0] = x0;
    m_proj[2][1] = y0;
    m_proj[2][2] = z0;
  } else {
    std::cout << m_className << "::SetPlane:\n"
              << "    Normal vector has zero norm. No new projection set.\n";
  }

  // Store the plane description
  m_plane[0] = fx;
  m_plane[1] = fy;
  m_plane[2] = fz;
  m_plane[3] = fx * x0 + fy * y0 + fz * z0;

  // Make labels to be placed along the axes
  Labels();
}

void ViewField::Rotate(const double theta) {
  // Rotate the axes
  double auxu[3], auxv[3];
  const double ctheta = cos(theta);
  const double stheta = sin(theta);
  for (int i = 0; i < 3; ++i) {
    auxu[i] = ctheta * m_proj[0][i] - stheta * m_proj[1][i];
    auxv[i] = stheta * m_proj[0][i] + ctheta * m_proj[1][i];
  }
  for (int i = 0; i < 3; ++i) {
    m_proj[0][i] = auxu[i];
    m_proj[1][i] = auxv[i];
  }

  // Make labels to be placed along the axes
  Labels();
}

ViewField::PlotType ViewField::GetPlotType(const std::string& option,
                                           std::string& title) const {

  std::string opt;
  std::transform(option.begin(), option.end(), 
                 std::back_inserter(opt), toupper);
  if (opt == "V" || opt == "P" || opt == "PHI" || 
      opt.find("VOLT") != std::string::npos ||
      opt.find("POT") != std::string::npos) {
    title = "potential";
    return PlotType::Potential;
  } else if (opt == "E" || opt == "FIELD" || opt == "NORM" ||
             opt.find("MAG") != std::string::npos) {
    title = "field";
    return PlotType::Magnitude;
  } else if (opt.find("X") != std::string::npos) {
    title = "field (x-component)";
    return PlotType::Ex;
  } else if (opt.find("Y") != std::string::npos) {
    title = "field (y-component)";
    return PlotType::Ey;
  } else if (opt.find("Z") != std::string::npos) {
    title = "field (z-component)";
    return PlotType::Ez;
  }
  std::cerr << m_className << "::GetPlotType:\n    Unknown option (" << option
            << ").\n";
  title = "potential";
  return PlotType::Potential;
}

void ViewField::Draw2d(const std::string& option, const bool contour, 
                       const bool wfield, const std::string& electrode,
                       const std::string& drawopt) {
  if (!m_sensor && !m_component) {
    std::cerr << m_className << "::Draw2d:\n"
              << "    Neither sensor nor component are defined.\n";
    return;
  }

  // Determine the x-y range (unless specified explicitly by the user).
  if (!Range()) return;

  // Determine the quantity to be plotted.
  std::string title;
  const PlotType plotType = GetPlotType(option, title);

  auto eval = [this, plotType, wfield, electrode](double* u, double* /*p*/) {
    // Transform to global coordinates.
    const double x = m_proj[0][0] * u[0] + m_proj[1][0] * u[1] + m_proj[2][0];
    const double y = m_proj[0][1] * u[0] + m_proj[1][1] * u[1] + m_proj[2][1];
    const double z = m_proj[0][2] * u[0] + m_proj[1][2] * u[1] + m_proj[2][2];
    return wfield ? Wfield(x, y, z, plotType, electrode) : Field(x, y, z, plotType);
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
    } else if (plotType == PlotType::Potential) {
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
    if (plotType == PlotType::Potential) {
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
  std::string labels = ";" + std::string(m_xLabel) + 
                       ";" + std::string(m_yLabel);
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
  const PlotType plotType = GetPlotType(option, title);

  auto eval = [this, plotType, wfield, electrode, 
               x0, y0, z0, dx, dy, dz](double* u, double* /*p*/) {
    // Get the position.
    const double t = u[0];
    const double x = x0 + t * dx;
    const double y = y0 + t * dy;
    const double z = z0 + t * dz;
    return wfield ? Wfield(x, y, z, plotType, electrode) : Field(x, y, z, plotType);
  };

  const std::string fname = FindUnusedFunctionName("fProfile");
  TF1 f1(fname.c_str(), eval, 0., 1., 0);

  double fmin = m_vmin;
  double fmax = m_vmax;
  if (wfield) {
    title = "weighting " + title;
    if (plotType == PlotType::Potential) {
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
    if (plotType == PlotType::Potential) {
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
  if (plotType == PlotType::Potential) {
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
      if (plotType != PlotType::Magnitude) labels += ",";
    } else if (plotType != PlotType::Magnitude) {
      labels += "_{";
    }
    if (plotType == PlotType::Ex) {
      labels += "x";
    } else if (plotType == PlotType::Ey) {
      labels += "y";
    } else if (plotType == PlotType::Ez) {
      labels += "z";
    }
    if (wfield || plotType != PlotType::Magnitude) labels += "}";
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

bool ViewField::Range() {

  if (m_userPlotLimits) return true;
  // Try to get the area/bounding box from the sensor/component.
  double bbmin[3];
  double bbmax[3];
  if (m_sensor) {
    if (!m_sensor->GetArea(bbmin[0], bbmin[1], bbmin[2], bbmax[0], bbmax[1],
                           bbmax[2])) {
      std::cerr << m_className << "::Range:\n"
                << "    Sensor area is not defined.\n"
                << "    Please set the plot range explicitly (SetArea).\n";
      return false;
    }
  } else {
    if (!m_component->GetBoundingBox(bbmin[0], bbmin[1], bbmin[2], bbmax[0],
                                     bbmax[1], bbmax[2])) {
      std::cerr << m_className << "::Range:\n"
                << "    Bounding box of the component is not defined.\n"
                << "    Please set the plot range explicitly (SetArea).\n";
      return false;
    }
  }
  const double tol = 1.e-4;
  double umin[2] = {-std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max()};
  double umax[2] = {std::numeric_limits<double>::max(),
                    std::numeric_limits<double>::max()};
  for (unsigned int i = 0; i < 3; ++i) {
    bbmin[i] -= m_proj[2][i];
    bbmax[i] -= m_proj[2][i];
    for (unsigned int j = 0; j < 2; ++j) {
      if (fabs(m_proj[j][i]) < tol) continue;
      const double t1 = bbmin[i] / m_proj[j][i];
      const double t2 = bbmax[i] / m_proj[j][i];
      const double tmin = std::min(t1, t2);
      const double tmax = std::max(t1, t2);
      if (tmin > umin[j] && tmin < umax[j]) umin[j] = tmin;
      if (tmax < umax[j] && tmax > umin[j]) umax[j] = tmax;
    }
  }
  m_xMinPlot = umin[0];
  m_xMaxPlot = umax[0];
  m_yMinPlot = umin[1];
  m_yMaxPlot = umax[1];
  std::cout << m_className << "::Range: Setting plot range to "
            << m_xMinPlot << " < x < " << m_xMaxPlot << ", " 
            << m_yMinPlot << " < y < " << m_yMaxPlot << ".\n";
  return true;
}

double ViewField::Field(const double x, const double y, const double z,
                        const PlotType plotType) const {

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
  switch (plotType) {
    case PlotType::Potential:
      return volt;
      break;
    case PlotType::Magnitude:
      return sqrt(ex * ex + ey * ey + ez * ez);
      break;
    case PlotType::Ex:
      return ex;
      break;
    case PlotType::Ey:
      return ey;
      break;
    case PlotType::Ez:
      return ez;
      break;
    default:
      break;
  }
  return volt;
}

double ViewField::Wfield(const double x, const double y, const double z,
                         const PlotType plotType, 
                         const std::string& electrode) const {

  if (plotType == PlotType::Potential) {
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

  switch (plotType) {
    case PlotType::Magnitude:
      return sqrt(ex * ex + ey * ey + ez * ez);
      break;
    case PlotType::Ex:
      return ex;
      break;
    case PlotType::Ey:
      return ey;
      break;
    case PlotType::Ez:
      return ez;
      break;
    default:
      break;
  }
  return 0.;
} 

}
