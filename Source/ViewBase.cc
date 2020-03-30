#include <iostream>
#include <cstring>
#include <cmath>

#include <TROOT.h>

#include "Garfield/Plotting.hh"
#include "Garfield/ViewBase.hh"

namespace {

bool Invert(std::array<std::array<double, 3>, 3>& a) {

  // Compute cofactors.
  const double c11 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
  const double c12 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
  const double c13 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
  const double c21 = a[2][1] * a[0][2] - a[2][2] * a[0][1];
  const double c22 = a[2][2] * a[0][0] - a[2][0] * a[0][2];
  const double c23 = a[2][0] * a[0][1] - a[2][1] * a[0][0];
  const double c31 = a[0][1] * a[1][2] - a[0][2] * a[1][1];
  const double c32 = a[0][2] * a[1][0] - a[0][0] * a[1][2];
  const double c33 = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  const double t1 = fabs(a[0][0]);
  const double t2 = fabs(a[1][0]);
  const double t3 = fabs(a[2][0]);
  double det = 0.;
  double pivot = 0.;
  if (t2 < t1 && t3 < t1) {
    pivot = a[0][0];
    det = c22 * c33 - c23 * c32;
  } else if (t1 < t2 && t3 < t2) {
    pivot = a[1][0];
    det = c13 * c32 - c12 * c33;
  } else {
    pivot = a[2][0];
    det = c23 * c12 - c22 * c13;
  }
  if (det == 0.) return false;
  const double s = pivot / det;
  a[0][0] = s * c11;
  a[0][1] = s * c21;
  a[0][2] = s * c31;
  a[1][0] = s * c12;
  a[1][1] = s * c22;
  a[1][2] = s * c32;
  a[2][0] = s * c13;
  a[2][1] = s * c23;
  a[2][2] = s * c33;
  return true;
}

}

namespace Garfield {

ViewBase::ViewBase(const std::string& name) :
    m_className(name) { 

  plottingEngine.SetDefaultStyle();
}

ViewBase::~ViewBase() {
  if (!m_hasExternalCanvas && m_canvas) delete m_canvas;
}

void ViewBase::SetCanvas(TCanvas* c) {
  if (!c) return;
  if (!m_hasExternalCanvas && m_canvas) {
    delete m_canvas;
    m_canvas = nullptr;
  }
  m_canvas = c;
  m_hasExternalCanvas = true;
}

TCanvas* ViewBase::GetCanvas() {
  if (!m_canvas) {
    m_canvas = new TCanvas();
    if (m_hasExternalCanvas) m_hasExternalCanvas = false;
  }
  return m_canvas;
}

void ViewBase::SetArea(const double xmin, const double ymin, 
                       const double xmax, const double ymax) {
  // Check range, assign if non-null.
  if (xmin == xmax || ymin == ymax) {
    std::cerr << m_className << "::SetArea: Null area is not permitted.\n"
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

void ViewBase::SetArea(const double xmin, const double ymin, const double zmin,
                       const double xmax, const double ymax,
                       const double zmax) {
  // Check range, assign if non-null
  if (xmin == xmax || ymin == ymax || zmin == zmax) {
    std::cerr << m_className << "::SetArea: Null area range not permitted.\n";
    return;
  }
  m_xMinBox = std::min(xmin, xmax);
  m_yMinBox = std::min(ymin, ymax);
  m_zMinBox = std::min(zmin, zmax);
  m_xMaxBox = std::max(xmin, xmax);
  m_yMaxBox = std::max(ymin, ymax);
  m_zMaxBox = std::max(zmin, zmax);
  m_userBox = true;
}

void ViewBase::SetPlane(const double fx, const double fy, const double fz,
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

  UpdateProjectionMatrix();
}

void ViewBase::Rotate(const double theta) {
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

  UpdateProjectionMatrix();
}


std::string ViewBase::FindUnusedFunctionName(const std::string& s) const {
  int idx = 0;
  std::string fname = s + "_0";
  while (gROOT->GetListOfFunctions()->FindObject(fname.c_str())) {
    ++idx;
    fname = s + "_" + std::to_string(idx);
  }
  return fname;
}

std::string ViewBase::FindUnusedHistogramName(const std::string& s) const {
  int idx = 0;
  std::string hname = s + "_0";
  while (gDirectory->GetList()->FindObject(hname.c_str())) {
    ++idx;
    hname = s + "_" + std::to_string(idx);
  }
  return hname;
}

void ViewBase::UpdateProjectionMatrix() {

  m_prmat[0][0] = m_proj[0][0];
  m_prmat[1][0] = m_proj[0][1];
  m_prmat[2][0] = m_proj[0][2];
  m_prmat[0][1] = m_proj[1][0];
  m_prmat[1][1] = m_proj[1][1];
  m_prmat[2][1] = m_proj[1][2];
  const double vnorm = sqrt(m_plane[0] * m_plane[0] +
                            m_plane[1] * m_plane[1] +
                            m_plane[2] * m_plane[2]);
  if (vnorm <= 0.) {
    std::cerr << m_className << "::UpdateProjectionMatrix:\n"
              << "    Zero norm vector.\n";
    // TODO! Reset to default.
    return;
  }
  m_prmat[0][2] = m_plane[0] / vnorm;
  m_prmat[1][2] = m_plane[1] / vnorm;
  m_prmat[2][2] = m_plane[2] / vnorm;
  if (!Invert(m_prmat)) {
    std::cerr << m_className << "::UpdateProjectionMatrix:\n"
              << "    Inversion failed.\n";
  }
}

void ViewBase::ToPlane(const double x, const double y, const double z,
                       double& xp, double& yp) const {
  xp = m_prmat[0][0] * x + m_prmat[0][1] * y + m_prmat[0][2] * z;
  yp = m_prmat[1][0] * x + m_prmat[1][1] * y + m_prmat[1][2] * z;
}

void ViewBase::ToPlane(const float x, const float y, const float z,
                       float& xp, float& yp) const {
  xp = m_prmat[0][0] * x + m_prmat[0][1] * y + m_prmat[0][2] * z;
  yp = m_prmat[1][0] * x + m_prmat[1][1] * y + m_prmat[1][2] * z;
}

std::string ViewBase::LabelX() {
  // Initialisation of the x-axis label.
  char xLabel[50];
  strcpy(xLabel, "\0");
  char buf[100];

  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_proj[0][0] - 1) < tol) {
    strcat(xLabel, "x");
  } else if (fabs(m_proj[0][0] + 1) < tol) {
    strcat(xLabel, "-x");
  } else if (m_proj[0][0] > tol) {
    sprintf(buf, "%g x", m_proj[0][0]);
    strcat(xLabel, buf);
  } else if (m_proj[0][0] < -tol) {
    sprintf(buf, "%g x", m_proj[0][0]);
    strcat(xLabel, buf);
  }

  // y portion
  if (strlen(xLabel) > 0) {
    if (m_proj[0][1] < -tol) {
      strcat(xLabel, " - ");
    } else if (m_proj[0][1] > tol) {
      strcat(xLabel, " + ");
    }
    if (fabs(m_proj[0][1] - 1) < tol || fabs(m_proj[0][1] + 1) < tol) {
      strcat(xLabel, "y");
    } else if (fabs(m_proj[0][1]) > tol) {
      sprintf(buf, "%g y", fabs(m_proj[0][1]));
      strcat(xLabel, buf);
    }
  } else {
    if (fabs(m_proj[0][1] - 1) < tol) {
      strcat(xLabel, "y");
    } else if (fabs(m_proj[0][1] + 1) < tol) {
      strcat(xLabel, "-y");
    } else if (m_proj[0][1] > tol) {
      sprintf(buf, "%g y", m_proj[0][1]);
      strcat(xLabel, buf);
    } else if (m_proj[0][1] < -tol) {
      sprintf(buf, "%g y", m_proj[0][1]);
      strcat(xLabel, buf);
    }
  }

  // z portion
  if (strlen(xLabel) > 0) {
    if (m_proj[0][2] < -tol) {
      strcat(xLabel, " - ");
    } else if (m_proj[0][2] > tol) {
      strcat(xLabel, " + ");
    }
    if (fabs(m_proj[0][2] - 1) < tol || fabs(m_proj[0][2] + 1) < tol) {
      strcat(xLabel, "z");
    } else if (fabs(m_proj[0][2]) > tol) {
      sprintf(buf, "%g z", fabs(m_proj[0][2]));
      strcat(xLabel, buf);
    }
  } else {
    if (fabs(m_proj[0][2] - 1) < tol) {
      strcat(xLabel, "z");
    } else if (fabs(m_proj[0][2] + 1) < tol) {
      strcat(xLabel, "-z");
    } else if (m_proj[0][2] > tol) {
      sprintf(buf, "%g z", m_proj[0][2]);
      strcat(xLabel, buf);
    } else if (m_proj[0][2] < -tol) {
      sprintf(buf, "%g z", m_proj[0][2]);
      strcat(xLabel, buf);
    }
  }

  // Unit
  strcat(xLabel, " [cm]");
  return std::string(xLabel);

}

std::string ViewBase::LabelY() {

  char yLabel[50];
  // Initialisation of the y-axis label
  strcpy(yLabel, "\0");
  char buf[100];

  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_proj[1][0] - 1) < tol) {
    strcat(yLabel, "x");
  } else if (fabs(m_proj[1][0] + 1) < tol) {
    strcat(yLabel, "-x");
  } else if (m_proj[1][0] > tol) {
    sprintf(buf, "%g x", m_proj[1][0]);
    strcat(yLabel, buf);
  } else if (m_proj[1][0] < -tol) {
    sprintf(buf, "%g x", m_proj[1][0]);
    strcat(yLabel, buf);
  }

  // y portion
  if (strlen(yLabel) > 0) {
    if (m_proj[1][1] < -tol) {
      strcat(yLabel, " - ");
    } else if (m_proj[1][1] > tol) {
      strcat(yLabel, " + ");
    }
    if (fabs(m_proj[1][1] - 1) < tol || fabs(m_proj[1][1] + 1) < tol) {
      strcat(yLabel, "y");
    } else if (fabs(m_proj[1][1]) > tol) {
      sprintf(buf, "%g y", fabs(m_proj[1][1]));
      strcat(yLabel, buf);
    }
  } else {
    if (fabs(m_proj[1][1] - 1) < tol) {
      strcat(yLabel, "y");
    } else if (fabs(m_proj[1][1] + 1) < tol) {
      strcat(yLabel, "-y");
    } else if (m_proj[1][1] > tol) {
      sprintf(buf, "%g y", m_proj[1][1]);
      strcat(yLabel, buf);
    } else if (m_proj[1][1] < -tol) {
      sprintf(buf, "%g y", m_proj[1][1]);
      strcat(yLabel, buf);
    }
  }

  // z portion
  if (strlen(yLabel) > 0) {
    if (m_proj[1][2] < -tol) {
      strcat(yLabel, " - ");
    } else if (m_proj[1][2] > tol) {
      strcat(yLabel, " + ");
    }
    if (fabs(m_proj[1][2] - 1) < tol || fabs(m_proj[1][2] + 1) < tol) {
      strcat(yLabel, "z");
    } else if (fabs(m_proj[1][2]) > tol) {
      sprintf(buf, "%g z", fabs(m_proj[1][2]));
      strcat(yLabel, buf);
    }
  } else {
    if (fabs(m_proj[1][2] - 1) < tol) {
      strcat(yLabel, "z");
    } else if (fabs(m_proj[1][2] + 1) < tol) {
      strcat(yLabel, "-z");
    } else if (m_proj[1][2] > tol) {
      sprintf(buf, "%g z", m_proj[1][2]);
      strcat(yLabel, buf);
    } else if (m_proj[1][2] < -tol) {
      sprintf(buf, "%g z", m_proj[1][2]);
      strcat(yLabel, buf);
    }
  }

  // Unit
  strcat(yLabel, " [cm]");
  return std::string(yLabel);
}

std::string ViewBase::PlaneDescription() {

  char description[50];
  // Initialisation of the plane label
  strcpy(description, "\0");
  char buf[100];

  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_plane[0] - 1) < tol) {
    strcat(description, "x");
  } else if (fabs(m_plane[0] + 1) < tol) {
    strcat(description, "-x");
  } else if (m_plane[0] > tol) {
    sprintf(buf, "%g x", m_plane[0]);
    strcat(description, buf);
  } else if (m_plane[0] < -tol) {
    sprintf(buf, "%g x", m_plane[0]);
    strcat(description, buf);
  }

  // y portion
  if (strlen(description) > 0) {
    if (m_plane[1] < -tol) {
      strcat(description, " - ");
    } else if (m_plane[1] > tol) {
      strcat(description, " + ");
    }
    if (fabs(m_plane[1] - 1) < tol || fabs(m_plane[1] + 1) < tol) {
      strcat(description, "y");
    } else if (fabs(m_plane[1]) > tol) {
      sprintf(buf, "%g y", fabs(m_plane[1]));
      strcat(description, buf);
    }
  } else {
    if (fabs(m_plane[1] - 1) < tol) {
      strcat(description, "y");
    } else if (fabs(m_plane[1] + 1) < tol) {
      strcat(description, "-y");
    } else if (m_plane[1] > tol) {
      sprintf(buf, "%g y", m_plane[1]);
      strcat(description, buf);
    } else if (m_plane[1] < -tol) {
      sprintf(buf, "%g y", m_plane[1]);
      strcat(description, buf);
    }
  }

  // z portion
  if (strlen(description) > 0) {
    if (m_plane[2] < -tol) {
      strcat(description, " - ");
    } else if (m_plane[2] > tol) {
      strcat(description, " + ");
    }
    if (fabs(m_plane[2] - 1) < tol || fabs(m_plane[2] + 1) < tol) {
      strcat(description, "z");
    } else if (fabs(m_plane[2]) > tol) {
      sprintf(buf, "%g z", fabs(m_plane[2]));
      strcat(description, buf);
    }
  } else {
    if (fabs(m_plane[2] - 1) < tol) {
      strcat(description, "z");
    } else if (fabs(m_plane[2] + 1) < tol) {
      strcat(description, "-z");
    } else if (m_plane[2] > tol) {
      sprintf(buf, "%g z", m_plane[2]);
      strcat(description, buf);
    } else if (m_plane[2] < -tol) {
      sprintf(buf, "%g z", m_plane[2]);
      strcat(description, buf);
    }
  }

  // Constant
  sprintf(buf, " = %g", m_plane[3]);
  strcat(description, buf);
  return std::string(description);
}

}
