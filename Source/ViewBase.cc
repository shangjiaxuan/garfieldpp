#include <iostream>
#include <cstring>

#include <TROOT.h>

#include "Garfield/Plotting.hh"
#include "Garfield/ViewBase.hh"

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

  // Make labels to be placed along the axes
  // TODO!
  // Labels();
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

  // Make labels to be placed along the axes
  // Labels();
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
