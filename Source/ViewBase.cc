#include <iostream>
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

void ViewBase::SetPlaneXY() {
  m_proj = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 0}}};
  m_plane = {0, 0, 1, 0};
  m_prmat = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
}

void ViewBase::SetPlaneXZ() {
  m_proj = {{{1, 0, 0}, {0, 0, 1}, {0, 0, 0}}};
  m_plane = {0, 1, 0, 0};
  UpdateProjectionMatrix();
}

void ViewBase::SetPlaneYZ() {
  m_proj = {{{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};
  m_plane = {1, 0, 0, 0};
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
              << "    Zero norm vector; reset to default.\n";
    SetPlaneXY();
    return;
  }
  m_prmat[0][2] = m_plane[0] / vnorm;
  m_prmat[1][2] = m_plane[1] / vnorm;
  m_prmat[2][2] = m_plane[2] / vnorm;
  if (!Invert(m_prmat)) {
    std::cerr << m_className << "::UpdateProjectionMatrix:\n"
              << "    Inversion failed; reset to default.\n";
    SetPlaneXY();
  }
}

void ViewBase::ToPlane(const double x, const double y, const double z,
                       double& xp, double& yp) const {
  xp = m_prmat[0][0] * x + m_prmat[0][1] * y + m_prmat[0][2] * z;
  yp = m_prmat[1][0] * x + m_prmat[1][1] * y + m_prmat[1][2] * z;
}

std::string ViewBase::LabelX() {

  std::string xLabel = "";
  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_proj[0][0] - 1) < tol) {
    xLabel = "#it{x}";
  } else if (fabs(m_proj[0][0] + 1) < tol) {
    xLabel = "#minus#it{x}";
  } else if (fabs(m_proj[0][0]) > tol) {
    xLabel = std::to_string(m_proj[0][0]) + " #it{x}";
  }

  // y portion
  if (!xLabel.empty()) {
    if (m_proj[0][1] < -tol) {
      xLabel += " #minus ";
    } else if (m_proj[0][1] > tol) {
      xLabel += " #plus ";
    }
    if (fabs(m_proj[0][1] - 1) < tol || fabs(m_proj[0][1] + 1) < tol) {
      xLabel += "#it{y}";
    } else if (fabs(m_proj[0][1]) > tol) {
      xLabel += std::to_string(fabs(m_proj[0][1])) + " #it{y}";
    }
  } else {
    if (fabs(m_proj[0][1] - 1) < tol) {
      xLabel = "#it{y}";
    } else if (fabs(m_proj[0][1] + 1) < tol) {
      xLabel = "#minus#it{y}";
    } else if (fabs(m_proj[0][1]) > tol) {
      xLabel = std::to_string(m_proj[0][1]) + " #it{y}";
    }
  }

  // z portion
  if (!xLabel.empty()) {
    if (m_proj[0][2] < -tol) {
      xLabel += " #minus ";
    } else if (m_proj[0][2] > tol) {
      xLabel += " #plus ";
    }
    if (fabs(m_proj[0][2] - 1) < tol || fabs(m_proj[0][2] + 1) < tol) {
      xLabel += "#it{z}";
    } else if (fabs(m_proj[0][2]) > tol) {
      xLabel += std::to_string(fabs(m_proj[0][2])) + " #it{z}";
    }
  } else {
    if (fabs(m_proj[0][2] - 1) < tol) {
      xLabel = "#it{z}";
    } else if (fabs(m_proj[0][2] + 1) < tol) {
      xLabel = "#minus#it{z}";
    } else if (fabs(m_proj[0][2]) > tol) {
      xLabel = std::to_string(m_proj[0][2]) + " #it{z}";
    }
  }

  // Unit
  xLabel += " [cm]";
  return xLabel;

}

std::string ViewBase::LabelY() {

  std::string yLabel = "";
  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_proj[1][0] - 1) < tol) {
    yLabel = "#it{x}";
  } else if (fabs(m_proj[1][0] + 1) < tol) {
    yLabel = "#minus#it{x}";
  } else if (fabs(m_proj[1][0]) > tol) {
    yLabel = std::to_string(m_proj[1][0]) + " #it{x}";
  }

  // y portion
  if (!yLabel.empty()) {
    if (m_proj[1][1] < -tol) {
      yLabel += " #minus ";
    } else if (m_proj[1][1] > tol) {
      yLabel += " #plus ";
    }
    if (fabs(m_proj[1][1] - 1) < tol || fabs(m_proj[1][1] + 1) < tol) {
      yLabel += "#it{y}";
    } else if (fabs(m_proj[1][1]) > tol) {
      yLabel += std::to_string(fabs(m_proj[1][1])) + " #it{y}";
    }
  } else {
    if (fabs(m_proj[1][1] - 1) < tol) {
      yLabel = "#it{y}";
    } else if (fabs(m_proj[1][1] + 1) < tol) {
      yLabel = "#minus#it{y}";
    } else if (fabs(m_proj[1][1]) > tol) {
      yLabel = std::to_string(m_proj[1][1]) + " #it{y}";
    }
  }

  // z portion
  if (!yLabel.empty()) {
    if (m_proj[1][2] < -tol) {
      yLabel += " #minus ";
    } else if (m_proj[1][2] > tol) {
      yLabel += " #plus ";
    }
    if (fabs(m_proj[1][2] - 1) < tol || fabs(m_proj[1][2] + 1) < tol) {
      yLabel += "#it{z}";
    } else if (fabs(m_proj[1][2]) > tol) {
      yLabel += std::to_string(fabs(m_proj[1][2])) + " #it{z}";
    }
  } else {
    if (fabs(m_proj[1][2] - 1) < tol) {
      yLabel = "#it{z}";
    } else if (fabs(m_proj[1][2] + 1) < tol) {
      yLabel = "#minus#it{z}";
    } else if (fabs(m_proj[1][2]) > tol) {
      yLabel = std::to_string(m_proj[1][2]) + " #it{z}";
    }
  }

  // Unit
  yLabel += " [cm]";
  return yLabel;
}

std::string ViewBase::PlaneDescription() {

  std::string description;

  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_plane[0] - 1) < tol) {
    description = "x";
  } else if (fabs(m_plane[0] + 1) < tol) {
    description = "-x";
  } else if (fabs(m_plane[0]) > tol) {
    description = std::to_string(m_plane[0]) + " x";
  }

  // y portion
  if (!description.empty()) {
    if (m_plane[1] < -tol) {
      description += " - ";
    } else if (m_plane[1] > tol) {
      description += " + ";
    }
    if (fabs(m_plane[1] - 1) < tol || fabs(m_plane[1] + 1) < tol) {
      description += "y";
    } else if (fabs(m_plane[1]) > tol) {
      description += std::to_string(fabs(m_plane[1])) + " y";
    }
  } else {
    if (fabs(m_plane[1] - 1) < tol) {
      description = "y";
    } else if (fabs(m_plane[1] + 1) < tol) {
      description = "-y";
    } else if (fabs(m_plane[1]) > tol) {
      description = std::to_string(m_plane[1]) + " y";
    }
  }

  // z portion
  if (!description.empty()) {
    if (m_plane[2] < -tol) {
      description += " - ";
    } else if (m_plane[2] > tol) {
      description += " + ";
    }
    if (fabs(m_plane[2] - 1) < tol || fabs(m_plane[2] + 1) < tol) {
      description += "z";
    } else if (fabs(m_plane[2]) > tol) {
      description += std::to_string(fabs(m_plane[2])) + " z";
    }
  } else {
    if (fabs(m_plane[2] - 1) < tol) {
      description = "z";
    } else if (fabs(m_plane[2] + 1) < tol) {
      description = "-z";
    } else if (fabs(m_plane[2]) > tol) {
      description = std::to_string(m_plane[2]) + " z";
    }
  }

  // Constant
  description += " = " + std::to_string(m_plane[3]);
  return description;
}

}
