#include <cmath>
#include <iostream>
#include <limits>
#include <set>

#include <TAxis.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TH1F.h>

#include "Garfield/Sensor.hh"
#include "Garfield/ComponentBase.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewIsochrons.hh"

namespace {

double Interpolate(const std::vector<double>& y,
                   const std::vector<double>& x, const double xx) {

  const double tol = 1.e-6 * fabs(x.back() - x.front());
  if (xx < x.front()) return y.front();
  const auto it1 = std::upper_bound(x.cbegin(), x.cend(), xx);
  if (it1 == x.cend()) return y.back();
  const auto it0 = std::prev(it1);
  const double dx = (*it1 - *it0);
  if (dx < tol) return y[it0 - x.cbegin()];
  const double f = (xx - *it0) / dx;
  return y[it0 - x.cbegin()] * (1. - f) + f * y[it1 - x.cbegin()];
}

bool OnLine(const double x1, const double y1, const double x2, const double y2,
            const double u, const double v) {
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v)});
  epsx = std::max(1.e-10, epsx);
  epsy = std::max(1.e-10, epsy);

  if ((fabs(x1 - u) <= epsx && fabs(y1 - v) <= epsy) ||
      (fabs(x2 - u) <= epsx && fabs(y2 - v) <= epsy)) {
    // Point to be examined coincides with start or end.
    return true;
  } else if (fabs(x1 - x2) <= epsx && fabs(y1 - y2) <= epsy) {
    // The line (x1, y1) to (x2, y2) is in fact a point.
    return false;
  }
  double xc = 0., yc = 0.;
  if (fabs(u - x1) + fabs(v - y1) < fabs(u - x2) + fabs(v - y2)) {
    // (u, v) is nearer to (x1, y1).
    const double dx = (x2 - x1);
    const double dy = (y2 - y1);
    const double xl = ((u - x1) * dx + (v - y1) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x1;
      yc = y1;
    } else if (xl > 1.) {
      xc = x2;
      yc = y2;
    } else {
      xc = x1 + xl * dx;
      yc = y1 + xl * dy;
    }
  } else {
    // (u, v) is nearer to (x2, y2).
    const double dx = (x1 - x2);
    const double dy = (y1 - y2);
    const double xl = ((u - x2) * dx + (v - y2) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x2;
      yc = y2;
    } else if (xl > 1.) {
      xc = x1;
      yc = y1;
    } else {
      xc = x2 + xl * dx;
      yc = y2 + xl * dy;
    }
  }
  // See whether the point is on the line.
  if (fabs(u - xc) < epsx && fabs(v - yc) < epsy) {
    return true;
  }
  return false;
}

bool Crossing(const double x1, const double y1, const double x2,
              const double y2, const double u1, const double v1,
              const double u2, const double v2) {

  // Check for a point of one line located on the other line.
  if (OnLine(x1, y1, x2, y2, u1, v1) || OnLine(x1, y1, x2, y2, u2, v2) ||
      OnLine(u1, v1, u2, v2, x1, y1) || OnLine(u1, v1, u2, v2, x2, y2)) {
    return true;
  }
  // Matrix to compute the crossing point.
  std::array<std::array<double, 2>, 2> a;
  a[0][0] = y2 - y1;
  a[0][1] = v2 - v1;
  a[1][0] = x1 - x2;
  a[1][1] = u1 - u2;
  const double det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u1), fabs(u2)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v1), fabs(v2)});
  epsx = std::max(epsx, 1.e-10);
  epsy = std::max(epsy, 1.e-10);
  if (fabs(det) < epsx * epsy) {
    // Parallel, non-touching.
    return false;
  }
  // Crossing, non-trivial lines: solve crossing equations.
  const double aux = a[1][1];
  a[1][1] = a[0][0] / det;
  a[0][0] = aux / det;
  a[1][0] = -a[1][0] / det;
  a[0][1] = -a[0][1] / det;
  // Compute crossing point.
  const double xc = a[0][0] * (x1 * y2 - x2 * y1) + a[1][0] * (u1 * v2 - u2 * v1);
  const double yc = a[0][1] * (x1 * y2 - x2 * y1) + a[1][1] * (u1 * v2 - u2 * v1);
  // See whether the crossing point is on both lines.
  if (OnLine(x1, y1, x2, y2, xc, yc) && OnLine(u1, v1, u2, v2, xc, yc)) {
    // Intersecting lines.
    return true;
  }
  // Crossing point not on both lines.
  return false;
}

}

namespace Garfield {

ViewIsochrons::ViewIsochrons() : ViewBase("ViewIsochrons") {
  SetDefaultProjection();
}

void ViewIsochrons::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }

  m_sensor = s;
  m_component = nullptr;
}

void ViewIsochrons::SetComponent(ComponentBase* c) {
  if (!c) {
    std::cerr << m_className << "::SetComponent: Null pointer.\n";
    return;
  }

  m_component = c;
  m_sensor = nullptr;
}

void ViewIsochrons::SetArea(const double xmin, const double ymin, 
                            const double xmax, const double ymax) {
  // Check range, assign if non-null.
  if (xmin == xmax || ymin == ymax) {
    std::cerr << m_className << "::SetArea: Null area range is not permitted.\n"
              << "      " << xmin << " < x < " << xmax << "\n"
              << "      " << ymin << " < y < " << ymax << "\n";
    return;
  }
  m_xmin = std::min(xmin, xmax);
  m_ymin = std::min(ymin, ymax);
  m_xmax = std::max(xmin, xmax);
  m_ymax = std::max(ymin, ymax);
  m_hasUserArea = true;
}

void ViewIsochrons::SetAspectRatioSwitch(const double ar) {

  if (ar < 0.) {
    std::cerr << m_className << "::SetAspectRatioSwitch: Value must be > 0.\n";
    return;
  }
  m_aspectRatio = ar;
}

void ViewIsochrons::SetLoopThreshold(const double thr) {

  if (thr < 0.) {
    std::cerr << m_className << "::SetLoopThreshold: Value must be > 0.\n";
    return;
  }
  m_loopThreshold = thr;
}

void ViewIsochrons::SetConnectionThreshold(const double thr) {

  if (thr < 0.) {
    std::cerr << m_className << "::SetConnectionThreshold:\n"
              << "    Value must be > 0.\n";
    return;
  }
  m_connectionThreshold = thr;
}

void ViewIsochrons::PlotIsochrons(const double tstep,
    const std::vector<std::array<double, 3> >& points,
    const bool colour, const bool markers) {

  if (!m_sensor && !m_component) {
    std::cerr << m_className << "::PlotIsochrons:\n"
              << "    Neither sensor nor component are defined.\n";
    return;
  }
  SetupCanvas();
  if (!Range()) return;
  m_canvas->DrawFrame(m_xmin, m_ymin, m_xmax, m_ymax, ";x;y");
  //-----------------------------------------------------------------------
  //   DRFEQT - The main routine (DRFEQT) accumulates equal drift time data
  //   DRFEQP   which is plotted as a set of contours in the entry DRFEQP.
  //-----------------------------------------------------------------------
  std::vector<std::vector<std::array<double, 3> > > driftLines;
  std::vector<std::array<double, 3> > startPoints;  
  std::vector<std::array<double, 3> > endPoints;  
  std::vector<int> IXYT;
  // Accumulate drift lines.
  ComputeDriftLines(tstep, points, driftLines, startPoints, endPoints, IXYT);
  const unsigned int nDriftLines = driftLines.size();
  if (nDriftLines < 2) {
    std::cerr << m_className << "::PlotIsochrons: Too few drift lines.\n";
    return;
  }
  // Keep track of the largest number of contours.
  std::size_t nContours = 0;
  for (const auto& driftLine : driftLines) {
    nContours = std::max(nContours, driftLine.size());
  }
  if (nContours == 0) {
    std::cerr << m_className << "::PlotIsochrons: No contour lines.\n";
    return;
  }

  std::set<int> allStats;
  for (const auto stat : IXYT) allStats.insert(stat);

  // DRFEQP
  if (m_debug) {
    std::cout << m_className << "::PlotIsochrons:\n"
              << "    Drawing " << nContours << " contours, "
              << nDriftLines << " drift lines.\n";
    std::printf("    Connection threshold:   %10.3f\n", m_connectionThreshold);
    std::printf("    Aspect ratio threshold: %10.3f\n", m_aspectRatio);
    std::printf("    Loop closing threshold: %10.3f\n", m_loopThreshold);
    if (m_sortContours) {
      std::cout << "    Sort contours.\n";
    } else {
      std::cout << "    Do not sort contours.\n";
    }
    if (m_checkCrossings) {
      std::cout << "    Check for crossings.\n";
    } else {
      std::cout << "    Do not check for crossings.\n";
    }
    if (markers) {
      std::cout << "    Mark isochron points.\n";
    } else {
      std::cout << "    Draw isochron lines.\n";
    }
  }

  // Loop over the equal time contours.
  TGraph graph;
  graph.SetLineColor(kGray + 2);
  graph.SetMarkerColor(kGray + 2);
  graph.SetLineStyle(9);
  const double colRange = double(gStyle->GetNumberOfColors()) / nContours;
  for (unsigned int ic = 0; ic < nContours; ++ic) {
    if (colour) {
      const auto col = gStyle->GetColorPalette(int((ic + 0.99) * colRange));
      graph.SetLineColor(col);
      graph.SetMarkerColor(col);
    }
    graph.SetMarkerStyle(m_markerStyle);
    for (int stat : allStats) {
      std::vector<std::pair<std::array<double, 4>, unsigned int> > contour;
      // Loop over the drift lines, picking up the points when OK.
      for (unsigned int k = 0; k < nDriftLines; ++k) {
        const auto& dl = driftLines[k]; 
        // Reject any undesirable combinations.
        if (IXYT[k] != stat || ic >= dl.size()) continue;
        // Add the point to the contour line.
        std::array<double, 4> point = {dl[ic][0], dl[ic][1], dl[ic][2], 0.};
        contour.push_back(std::make_pair(point, k)); 
      }
      // Skip the plot of this contour if there are no points.
      if (contour.empty()) continue;
      bool circle = false;
      // If requested, sort the points on the contour line.
      if (m_sortContours && !markers) SortContour(contour, circle);
      // Plot this contour.
      if (markers) {
        // Simply mark the contours if this was requested.
        std::vector<double> xp;
        std::vector<double> yp;
        std::vector<double> zp;
        for (const auto& point : contour) {
          const double x = point.first[0];
          const double y = point.first[1];
          const double z = point.first[2];
          xp.push_back(m_proj[0][0] * x + m_proj[1][0] * y + z * m_plane[0]);
          yp.push_back(m_proj[0][1] * x + m_proj[1][1] * y + z * m_plane[1]);
          zp.push_back(m_proj[0][2] * x + m_proj[1][2] * y + z * m_plane[2]);
        }
        graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Psame");
        continue;
      }
      // Regular plotting.
      const double tolx = (m_xmax - m_xmin) * m_connectionThreshold;
      const double toly = (m_ymax - m_ymin) * m_connectionThreshold;
      // Flag to keep track if the segment is interrupted by a drift line
      // drift line or if it is too long.
      bool gap = false;
      // Flag for the first segment.
      bool firstgap = false;
      unsigned int ibegin = 0;
      // Coordinates to be plotted.
      std::vector<double> xp;
      std::vector<double> yp;
      std::vector<double> zp;
      const unsigned int nP = contour.size();
      for (unsigned int i = 0; i < nP; ++i) { 
        gap = false;
        const auto x0 = contour[i].first[0];
        const auto y0 = contour[i].first[1];
        const auto z0 = contour[i].first[2];
        xp.push_back(m_proj[0][0] * x0 + m_proj[1][0] * y0 + z0 * m_plane[0]);
        yp.push_back(m_proj[0][1] * x0 + m_proj[1][1] * y0 + z0 * m_plane[1]);
        zp.push_back(m_proj[0][2] * x0 + m_proj[1][2] * y0 + z0 * m_plane[2]);
        if (i == nP - 1) break; 
        const auto x1 = contour[i + 1].first[0];
        const auto y1 = contour[i + 1].first[1];
        // Reject contour segments which are long compared with AREA.
        if (fabs(x1 - x0) > tolx || fabs(y1 - y0) > toly) gap = true;
        // Get the indices of the drift lines corresponding 
        // to these two points on the contour line.
        const auto i0 = contour[i].second;
        const auto i1 = contour[i + 1].second; 
        // Set the BREAK flag if it crosses some stored drift line segment.
        if (m_checkCrossings && !gap) {
          for (unsigned int k = 0; k < nDriftLines; ++k) {
            const auto& dl = driftLines[k];
            for (unsigned int jc = 0; jc < dl.size(); ++jc) {
              if ((i0 == k || i1 == k) && (jc == ic || jc + 1 == ic)) {
                continue;
              }
              const auto& p0 = dl[jc];
              const auto& p1 = jc == dl.size() - 1 ? endPoints[k] : dl[jc + 1];
              if (Crossing(p0[0], p0[1], p1[0], p1[1], x0, y0, x1, y1)) {
                gap = true;
                break;
              }
            }
            if (gap) break;
            if ((i0 == k || i1 == k) && ic == 0) continue;
            const auto& p0 = startPoints[k];
            if (Crossing(p0[0], p0[1], dl[0][0], dl[0][1], 
                         x0, y0, x1, y1)) {
              gap = true;
              break;
            }
          }
        }
        // If there has been a break, plot what we have already.
        if (gap) {
          if (xp.size() > 1) {
            // Plot line.
            graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Lsame");
          } else if (ibegin != 0 || !circle) {
            // Plot markers.
            graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Psame");
          } else if (ibegin == 0) {
            // TODO!
            firstgap = true;
          }
          xp.clear();
          yp.clear();
          zp.clear();
          ibegin = i + 1;
        }
      }
      // Plot the remainder; if there is a break, put a * if FRSTBR is on.
      if (!gap && xp.size() > 1) {
        graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Lsame");
       } else if ((firstgap || !circle) && ibegin == nP) {
        // TODO!
        graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Psame");
      }
    }
  }

  // gPad->SetRightMargin(0.15);
  gPad->Update();
  graph.SetLineStyle(1);
  if (m_particle == Particle::Electron) {
    graph.SetLineColor(kOrange - 3);
  } else {
    graph.SetLineColor(kRed + 1);
  }
  for (unsigned int i = 0; i < nDriftLines; ++i) {
    std::vector<double> xp;
    std::vector<double> yp;
    const double x0 = startPoints[i][0];
    const double y0 = startPoints[i][1];
    const double z0 = startPoints[i][2];
    xp.push_back(m_proj[0][0] * x0 + m_proj[1][0] * y0 + z0 * m_plane[0]);
    yp.push_back(m_proj[0][1] * x0 + m_proj[1][1] * y0 + z0 * m_plane[1]);
    for (const auto& point : driftLines[i]) {
      const double x = point[0];
      const double y = point[1];
      const double z = point[2];
      xp.push_back(m_proj[0][0] * x + m_proj[1][0] * y + z * m_plane[0]);
      yp.push_back(m_proj[0][1] * x + m_proj[1][1] * y + z * m_plane[1]);
    }
    const double x1 = endPoints[i][0];
    const double y1 = endPoints[i][1];
    const double z1 = endPoints[i][2];
    xp.push_back(m_proj[0][0] * x1 + m_proj[1][0] * y1 + z1 * m_plane[0]);
    yp.push_back(m_proj[0][1] * x1 + m_proj[1][1] * y1 + z1 * m_plane[1]);
    graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Lsame");
  }
}

void ViewIsochrons::ComputeDriftLines(const double tstep,
  const std::vector<std::array<double, 3> >& points,
  std::vector<std::vector<std::array<double, 3> > >& driftLines,
  std::vector<std::array<double, 3> >& startPoints,
  std::vector<std::array<double, 3> >& endPoints,
  std::vector<int>& IXYT) {

  DriftLineRKF drift;
  Sensor sensor;
  if (m_sensor) {
    drift.SetSensor(m_sensor);
  } else {
    sensor.AddComponent(m_component);
    drift.SetSensor(&sensor);
  }
  for (const auto& point : points) {
    if (m_particle == Particle::Electron) {
      if (m_positive) {
        drift.DriftPositron(point[0], point[1], point[2], 0.);
      } else {
        drift.DriftElectron(point[0], point[1], point[2], 0.);
      }
    } else {
      if (m_positive) {
        drift.DriftIon(point[0], point[1], point[2], 0.);
      } else {
        drift.DriftNegativeIon(point[0], point[1], point[2], 0.);
      }
    }
    const unsigned int nu = drift.GetNumberOfDriftLinePoints();
    // Check that the drift line has enough steps.
    if (nu < 3) continue;
    int status = 0;
    double xf = 0., yf = 0., zf = 0., tf = 0.;
    drift.GetEndPoint(xf, yf, zf, tf, status); 
    std::vector<double> xu(nu, 0.);
    std::vector<double> yu(nu, 0.);
    std::vector<double> zu(nu, 0.);
    std::vector<double> tu(nu, 0.);
    for (unsigned int i = 0; i < nu; ++i) {
      drift.GetDriftLinePoint(i, xu[i], yu[i], zu[i], tu[i]);
    }
    // Find the number of points to be stored.
    const unsigned int nSteps = static_cast<unsigned int>(tu.back() / tstep);
    if (nSteps == 0) continue;
    std::vector<std::array<double, 3> > tab;
    // Interpolate at regular time intervals.
    for (unsigned int i = 0; i < nSteps; ++i) {
      const double t = (i + 1) * tstep;
      // tab.push_back(PLACO3(Interpolate(xu, tu, t),
      //                      Interpolate(yu, tu, t),
      //                      Interpolate(zu, tu, t)));
      std::array<double, 3> step = {Interpolate(xu, tu, t),
                                    Interpolate(yu, tu, t),
                                    Interpolate(zu, tu, t)};
      tab.push_back(step);
    }
    driftLines.push_back(std::move(tab));
    std::array<double, 3> start = {xu[0], yu[0], zu[0]};
    std::array<double, 3> end = {xu[nu - 1], yu[nu - 1], zu[nu - 1]};
    startPoints.push_back(std::move(start));
    endPoints.push_back(std::move(end));
    // Store the drift line return code.
    IXYT.push_back(status);
  }
}

void ViewIsochrons::SetDefaultProjection() {
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

void ViewIsochrons::Labels() {
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

void ViewIsochrons::SetPlane(const double fx, const double fy, const double fz,
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

void ViewIsochrons::Rotate(const double theta) {
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

void ViewIsochrons::SetupCanvas() {
  if (!m_canvas) {
    m_canvas = new TCanvas();
    m_canvas->SetTitle("Isochrons");
    m_hasExternalCanvas = false;
  }
  m_canvas->cd();
}

bool ViewIsochrons::Range() {

  if (m_hasUserArea) return true;
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
  m_xmin = umin[0];
  m_xmax = umax[0];
  m_ymin = umin[1];
  m_ymax = umax[1];
  return true; 
}

void ViewIsochrons::SortContour(
    std::vector<std::pair<std::array<double, 4>, unsigned int> >& contour,
    bool& circle) {

  if (contour.size() < 2) return;
  // First compute the centre of gravity.
  double xcog = 0.;
  double ycog = 0.;
  for (const auto& point : contour) {
    xcog += point.first[0];
    ycog += point.first[1];
  }
  const double scale = 1. / contour.size();
  xcog *= scale;
  ycog *= scale;
  // Compute angles wrt to the centre of gravity and principal axes.
  double sxx = 0.;
  double sxy = 0.;
  double syy = 0.;
  for (const auto& point : contour) {
    const double dx = point.first[0] - xcog;
    const double dy = point.first[1] - ycog;
    sxx += dx * dx;
    sxy += dx * dy;
    syy += dy * dy;
  }
  const double theta = 0.5 * atan2(2 * sxy, sxx - syy);
  const double ct = cos(theta);
  const double st = sin(theta);
  // Evaluate dispersions around the principal axes.
  sxx = 0.;
  syy = 0.;
  for (const auto& point : contour) {
    const double dx = point.first[0] - xcog;
    const double dy = point.first[1] - ycog;
    sxx += fabs(+ct * dx + st * dy);
    syy += fabs(-st * dx + ct * dy);
  }
  // Decide whether this is more linear or more circular.
  if (fabs(sxx) > m_aspectRatio * fabs(syy) || 
      fabs(syy) > m_aspectRatio * fabs(sxx)) {
    circle = false;
  } else {
    circle = true;
  }
  // Set a sorting coordinate accordingly.
  for (auto& point : contour) {
    const double dx = point.first[0] - xcog;
    const double dy = point.first[1] - ycog;
    point.first[3] = circle ? atan2(dy, dx) : ct * dx + st * dy;
  }
  // Sort the points.
  std::sort(contour.begin(), contour.end(), 
            [](const std::pair<std::array<double, 4>, int>& p1, 
               const std::pair<std::array<double, 4>, int>& p2) {
                 return p1.first[3] < p2.first[3]; }
           );
  if (!circle) return;
  // For circles, perhaps add the first point to the end of the list.
  // Compute breakpoint, total distance and maximum distance.
  double dsum = 0.;
  double dmax = -1.;
  unsigned int imax = 0;
  const unsigned int nPoints = contour.size();
  for (unsigned int j = 0; j < nPoints; ++j) {
    const auto& p1 = contour[j];
    const auto& p0 = j > 0 ? contour[j - 1] : contour.back();
    const double dx = p1.first[0] - p0.first[0];
    const double dy = p1.first[1] - p0.first[1];
    const double d = sqrt(dx * dx + dy * dy);
    if (j > 0) dsum += d;
    if (dmax < d) {
      dmax = d;
      imax = j;
    }
  }
  if (dmax < m_loopThreshold * dsum) {
    // If a true loop, close it.
    contour.push_back(contour[0]);
  } else {
    circle = false;
    if (imax > 0) {
      // Shift the points to make a line.
      std::rotate(contour.begin(), contour.begin() + imax, contour.end());
    }
  }
}

std::array<double, 3> ViewIsochrons::PLACO3(const double x1, const double y1, const double z1) {

  //-----------------------------------------------------------------------
  //   PLACO3 - Determines plane coordinates.
  //   (Last changed on 29/ 9/98.)
  //-----------------------------------------------------------------------

  std::vector<double> b = {x1, y1, z1};
  // Solve the equation.
  // Numerics::CERNLIB::dfeqn(3, FPRMAT, IPRMAT, b);
  // Return the solution.
  const double zcut = m_plane[0] * x1 + m_plane[1] * y1 + m_plane[2] * z1;
  std::array<double, 3> result = {b[0], b[1], zcut};
  return result;
}


}
