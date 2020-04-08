#include <iostream>
#include <cmath>
#include <limits>

#include <TGraph.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TAxis.h>
#include <TAxis3D.h>
#include <TH1F.h>

#include "Garfield/Plotting.hh"
#include "Garfield/ViewDrift.hh"

namespace Garfield {

ViewDrift::ViewDrift() : ViewBase("ViewDrift") {
  m_driftLines.reserve(1000);
  m_tracks.reserve(100);
  m_exc.reserve(1000);
  m_ion.reserve(1000);
  m_ion.reserve(1000);
}

ViewDrift::~ViewDrift() {
  Clear();
}

void ViewDrift::Clear() {
  m_driftLines.clear();
  m_tracks.clear();

  m_exc.clear();
  m_ion.clear();
  m_ion.clear();
}

void ViewDrift::SetClusterMarkerSize(const double size) {
  if (size > 0.) {
    m_markerSizeCluster = size;
  } else {
    std::cerr << m_className << "::SetClusterMarkerSize: Size must be > 0.\n";
  }
}

void ViewDrift::SetCollisionMarkerSize(const double size) {
  if (size > 0.) {
    m_markerSizeCollision = size;
  } else {
    std::cerr << m_className << "::SetCollisionMarkerSize: Size must be > 0.\n";
  }
}

void ViewDrift::NewElectronDriftLine(const unsigned int np, int& id,
                                     const float x0, const float y0,
                                     const float z0) {
  // Create a new electron drift line and add it to the list.
  std::array<float, 3> p = {x0, y0, z0};
  std::vector<std::array<float, 3> > dl(std::max(1U, np), p);
  m_driftLines.push_back(std::make_pair(std::move(dl), Particle::Electron));
  // Return the index of this drift line.
  id = m_driftLines.size() - 1;
}

void ViewDrift::NewHoleDriftLine(const unsigned int np, int& id,
                                 const float x0, const float y0,
                                 const float z0) {

  std::array<float, 3> p = {x0, y0, z0};
  std::vector<std::array<float, 3> > dl(std::max(1U, np), p);
  m_driftLines.push_back(std::make_pair(std::move(dl), Particle::Hole));
  // Return the index of this drift line.
  id = m_driftLines.size() - 1;
}

void ViewDrift::NewIonDriftLine(const unsigned int np, int& id, const float x0,
                                const float y0, const float z0) {

  std::array<float, 3> p = {x0, y0, z0};
  std::vector<std::array<float, 3> > dl(std::max(1U, np), p);
  m_driftLines.push_back(std::make_pair(std::move(dl), Particle::Ion));
  // Return the index of this drift line.
  id = m_driftLines.size() - 1;
}

void ViewDrift::AddPhoton(const float x0, const float y0, const float z0, 
                          const float x1, const float y1, const float z1) {
  std::array<float, 3> p0 = {x0, y0, z0};
  std::array<float, 3> p1 = {x1, y1, z1};
  m_photons.push_back({p0, p1});
}

void ViewDrift::NewChargedParticleTrack(const unsigned int np, int& id,
                                        const float x0, const float y0,
                                        const float z0) {
  // Create a new track and add it to the list.
  std::vector<std::array<float, 3> > track(std::max(1U, np));
  track[0] = {x0, y0, z0};
  m_tracks.push_back(std::move(track));
  // Return the index of this track.
  id = m_tracks.size() - 1;
}

void ViewDrift::SetDriftLinePoint(const unsigned int iL, const unsigned int iP,
                                  const float x, const float y,
                                  const float z) {
  if (iL >= m_driftLines.size() || iP >= m_driftLines[iL].first.size()) {
    std::cerr << m_className << "::SetDriftLinePoint: Index out of range.\n";
    return;
  }
  m_driftLines[iL].first[iP] = {x, y, z};
}

void ViewDrift::AddDriftLinePoint(const unsigned int iL, const float x,
                                  const float y, const float z) {
  if (iL >= m_driftLines.size()) {
    std::cerr << m_className << "::AddDriftLinePoint: Index out of range.\n";
    return;
  }
  std::array<float, 3> p = {x, y, z};
  m_driftLines[iL].first.push_back(std::move(p));
}

void ViewDrift::SetTrackPoint(const unsigned int iL, const unsigned int iP,
                              const float x, const float y, const float z) {
  if (iL >= m_tracks.size() || iP >= m_tracks[iL].size()) {
    std::cerr << m_className << "::SetTrackPoint: Index out of range.\n";
    return;
  }
  m_tracks[iL][iP] = {x, y, z};
}

void ViewDrift::AddTrackPoint(const unsigned int iL, const float x,
                              const float y, const float z) {
  if (iL >= m_tracks.size()) {
    std::cerr << m_className << "::AddTrackPoint: Index out of range.\n";
    return;
  }
  std::array<float, 3> p = {x, y, z};
  m_tracks[iL].push_back(std::move(p));
}

void ViewDrift::AddExcitation(const float x, const float y, const float z) {
  std::array<float, 3> p = {x, y, z};
  m_exc.push_back(std::move(p));
}

void ViewDrift::AddIonisation(const float x, const float y, const float z) {
  std::array<float, 3> p = {x, y, z};
  m_ion.push_back(std::move(p));
}

void ViewDrift::AddAttachment(const float x, const float y, const float z) {
  std::array<float, 3> p = {x, y, z};
  m_ion.push_back(std::move(p));
}

void ViewDrift::Plot(const bool twod, const bool axis) {
  if (twod) {
    Plot2d(axis);
  } else {
    Plot3d(axis);
  }
}

void ViewDrift::Plot2d(const bool axis) {
  auto canvas = GetCanvas();
  canvas->cd();
  canvas->SetTitle("Drift lines");
  // Check if the canvas range has already been set.
  bool rangeSet = true;
  if (canvas->GetListOfPrimitives()->GetSize() == 0 && 
      canvas->GetX1() == 0 && canvas->GetX2() == 1 && 
      canvas->GetY1() == 0 && canvas->GetY2() == 1) {
    rangeSet = false;
  }
  if (axis || !rangeSet) {
    // Determine the plot limits.
    if (!SetPlotLimits()) {
      std::cerr << m_className << "::Plot2d:\n"
                << "     Could not determine the plot limits.\n";
      return;
    }
  }
  if (axis) {
    auto frame = canvas->DrawFrame(m_xMinPlot, m_yMinPlot, 
                                   m_xMaxPlot, m_yMaxPlot);
    frame->GetXaxis()->SetTitle(LabelX().c_str());
    frame->GetYaxis()->SetTitle(LabelY().c_str());
  } else if (!rangeSet) {
    const double bm = canvas->GetBottomMargin();
    const double lm = canvas->GetLeftMargin();
    const double rm = canvas->GetRightMargin();
    const double tm = canvas->GetTopMargin();
    const double dx = m_xMaxPlot - m_xMinPlot;
    const double dy = m_yMaxPlot - m_yMinPlot;
    canvas->Range(m_xMinPlot - dx * (lm / (1. - rm - lm)),
                  m_yMinPlot - dy * (bm / (1. - tm - lm)),
                  m_xMaxPlot + dx * (rm / (1. - rm - lm)),
                  m_yMaxPlot + dy * (tm / (1. - tm - lm)));
  } 

  TGraph gr;
  gr.SetLineWidth(1);
  for (const auto& driftLine : m_driftLines) {
    const unsigned int nP = driftLine.first.size();
    if (nP < 2) continue;
    if (driftLine.second == Particle::Electron) {
      gr.SetLineColor(m_colElectron);
    } else if (driftLine.second == Particle::Hole) {
      gr.SetLineColor(m_colHole);
    } else {
      gr.SetLineColor(m_colIon);
    }
    std::vector<float> xgr;
    std::vector<float> ygr;
    auto x0 = driftLine.first[0];
    bool in0 = InBox(x0);
    if (in0) {
      double xp = 0., yp = 0.;
      ToPlane(x0[0], x0[1], x0[2], xp, yp);
      xgr.push_back(xp);
      ygr.push_back(yp);
    }
    for (unsigned int j = 1; j < nP; ++j) {
      auto x1 = driftLine.first[j];
      bool in1 = InBox(x1);
      if (in1 != in0) {
        double xp = 0., yp = 0.;
        std::array<float, 3> xc;
        Clip(x0, x1, xc);
        ToPlane(xc[0], xc[1], xc[2], xp, yp);
        xgr.push_back(xp);
        ygr.push_back(yp);
      } 
      if (in1) {
        double xp = 0., yp = 0.;
        ToPlane(x1[0], x1[1], x1[2], xp, yp);
        xgr.push_back(xp);
        ygr.push_back(yp);
      } else if (!xgr.empty()) {
        gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Lsame");
        xgr.clear();
        ygr.clear();
      }
      x0 = x1;
      in0 = in1;
    }
    if (!xgr.empty()) {
      gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Lsame");
    }
  }
  gPad->Update();

  gr.SetLineWidth(2);
  gr.SetLineColor(m_colTrack);
  for (const auto& track : m_tracks) {
    const auto nP = track.size();
    if (nP < 2) continue;
    std::vector<float> xgr;
    std::vector<float> ygr;
    auto x0 = track[0];
    bool in0 = InBox(x0);
    if (in0) {
      double xp = 0., yp = 0.;
      ToPlane(x0[0], x0[1], x0[2], xp, yp);
      xgr.push_back(xp);
      ygr.push_back(yp);
    }
    for (unsigned int j = 1; j < nP; ++j) {
      auto x1 = track[j];
      bool in1 = InBox(x1);
      if (in1 != in0) {
        double xp = 0., yp = 0.;
        std::array<float, 3> xc;
        Clip(x0, x1, xc);
        ToPlane(xc[0], xc[1], xc[2], xp, yp);
        xgr.push_back(xp);
        ygr.push_back(yp);
      } 
      if (in1) {
        double xp = 0., yp = 0.;
        ToPlane(x1[0], x1[1], x1[2], xp, yp);
        xgr.push_back(xp);
        ygr.push_back(yp);
      } else if (!xgr.empty()) {
        gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Lsame");
        xgr.clear();
        ygr.clear();
      }
      x0 = x1;
      in0 = in1;
    }
    if (!xgr.empty()) {
      gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Lsame");
    }
  }

  gr.SetLineColor(m_colPhoton);
  gr.SetLineStyle(2);
  for (const auto& photon : m_photons) {
    double xp0 = 0., yp0 = 0.;
    double xp1 = 0., yp1 = 0.;
    ToPlane(photon[0][0], photon[0][1], photon[0][2], xp0, yp0);
    ToPlane(photon[1][0], photon[1][1], photon[1][2], xp1, yp1);
    std::vector<float> xgr = {float(xp0), float(xp1)};
    std::vector<float> ygr = {float(yp0), float(yp1)};
    gr.DrawGraph(2, xgr.data(), ygr.data(), "Lsame"); 
  }

  gr.SetMarkerSize(m_markerSizeCollision);
  gr.SetMarkerStyle(20);
  if (!m_exc.empty()) {
    gr.SetMarkerColor(m_colExcitation);
    std::vector<float> xgr;
    std::vector<float> ygr;
    for (const auto& p : m_exc) {
      if (!InBox(p)) continue;
      double xp = 0., yp = 0.;
      ToPlane(p[0], p[1], p[2], xp, yp);
      xgr.push_back(xp);
      ygr.push_back(yp); 
    }
    if (!xgr.empty()) {
      gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Psame");
    }
  } 
  if (!m_ion.empty()) {
    gr.SetMarkerColor(m_colIonisation);
    std::vector<float> xgr;
    std::vector<float> ygr;
    for (const auto& p : m_ion) {
      if (!InBox(p)) continue;
      double xp = 0., yp = 0.;
      ToPlane(p[0], p[1], p[2], xp, yp);
      xgr.push_back(xp);
      ygr.push_back(yp); 
    }
    if (!xgr.empty()) {
      gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Psame");
    }
  }
  if (!m_att.empty()) {
    gr.SetMarkerColor(m_colAttachment);
    std::vector<float> xgr;
    std::vector<float> ygr;
    for (const auto& p : m_att) {
      if (!InBox(p)) continue;
      double xp = 0., yp = 0.;
      ToPlane(p[0], p[1], p[2], xp, yp);
      xgr.push_back(xp);
      ygr.push_back(yp); 
    }
    if (!xgr.empty()) {
      gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Psame");
    }
  }
 
  gPad->Update();
}

void ViewDrift::Clip(const std::array<float, 3>& x0, 
                     const std::array<float, 3>& x1,
                     std::array<float, 3>& xc) const {

  xc.fill(0.);
  const bool in0 = InBox(x0);
  const bool in1 = InBox(x1);
  if (in0 == in1) return;
  xc = in0 ? x1 : x0;
  const std::array<float, 3> dx = {x1[0] - x0[0], x1[1] - x0[1], 
                                   x1[2] - x0[2]};
  // Adjust x.
  if (dx[0] != 0. && (xc[0] < m_xMinBox || xc[0] > m_xMaxBox)) {
    const double b = xc[0] < m_xMinBox ? m_xMinBox : m_xMaxBox;
    const double s = (b - xc[0]) / dx[0];
    xc[0] = b;
    xc[1] += dx[1] * s;
    xc[2] += dx[2] * s;
  }
  if (dx[1] != 0. && (xc[1] < m_yMinBox || xc[1] > m_yMaxBox)) {
    const double b = xc[1] < m_yMinBox ? m_yMinBox : m_yMaxBox;
    const double s = (b - xc[1]) / dx[1];
    xc[0] += dx[0] * s;
    xc[1] = b;
    xc[2] += dx[2] * s;
  }
  // Adjust z.
  if (dx[2] != 0. && (xc[2] < m_zMinBox || xc[2] > m_zMaxBox)) {
    const double b = xc[2] < m_zMinBox ? m_zMinBox : m_zMaxBox;
    const double s = (b - xc[2]) / dx[2];
    xc[0] += dx[0] * s;
    xc[1] += dx[1] * s;
    xc[2] = b;
  }
}

void ViewDrift::Plot3d(const bool axis) {
  if (m_debug) std::cout << m_className << "::Plot: Plotting in 3D.\n";
  auto canvas = GetCanvas();
  if (axis) {
    if (!canvas->GetView()) {
      m_view.reset(TView::CreateView(1, 0, 0));
      m_view->SetRange(m_xMinBox, m_yMinBox, m_zMinBox, m_xMaxBox, m_yMaxBox, m_zMaxBox);
      m_view->ShowAxis();
      m_view->Top();
      canvas->SetView(m_view.get());
    }
  }
  for (const auto& driftLine : m_driftLines) {
    TPolyLine3D* pl = new TPolyLine3D(driftLine.first.size());
    for (const auto& p : driftLine.first) {
      pl->SetNextPoint(p[0], p[1], p[2]);
    }
    if (driftLine.second == Particle::Electron) {
      pl->SetLineColor(m_colElectron);
    } else if (driftLine.second == Particle::Hole) {
      pl->SetLineColor(m_colHole);
    } else {
      pl->SetLineColor(m_colIon);
    }
    pl->SetLineWidth(1);
    pl->SetBit(kCanDelete);
    pl->Draw("same");
  }

  for (const auto& track : m_tracks) {
    const unsigned int nPoints = track.size();
    TPolyMarker3D* pm = new TPolyMarker3D(nPoints, 20);
    TPolyLine3D* pl = new TPolyLine3D(nPoints);
    for (const auto& p : track) {
      pm->SetNextPoint(p[0], p[1], p[2]);
      pl->SetNextPoint(p[0], p[1], p[2]);
    }
    pm->SetMarkerColor(m_colTrack);
    pm->SetMarkerSize(m_markerSizeCluster);
    pm->SetBit(kCanDelete);
    pm->Draw("same");
    pl->SetLineColor(m_colTrack);
    pl->SetLineWidth(1);
    pl->SetBit(kCanDelete);
    pl->Draw("same");
  }
  if (!m_exc.empty()) {
    const unsigned int nP = m_exc.size();
    TPolyMarker3D* pm = new TPolyMarker3D(nP, 20);
    for (unsigned int i = 0; i < nP; ++i) {
      pm->SetPoint(i, m_exc[i][0], m_exc[i][1], m_exc[i][2]);
    }
    pm->SetMarkerColor(m_colExcitation);
    pm->SetMarkerSize(m_markerSizeCollision);
    pm->SetBit(kCanDelete);
    pm->Draw("same");
  }

  if (!m_ion.empty()) {
    const unsigned int nP = m_ion.size();
    TPolyMarker3D* pm = new TPolyMarker3D(nP, 20);
    for (unsigned int i = 0; i < nP; ++i) {
      pm->SetPoint(i, m_ion[i][0], m_ion[i][1], m_ion[i][2]);
    }
    pm->SetMarkerColor(m_colIonisation);
    pm->SetMarkerSize(m_markerSizeCollision);
    pm->SetBit(kCanDelete);
    pm->Draw("same");
  }

  if (!m_att.empty()) {
    const unsigned int nP = m_att.size();
    TPolyMarker3D* pm = new TPolyMarker3D(nP, 20);
    for (unsigned int i = 0; i < nP; ++i) {
      pm->SetPoint(i, m_att[i][0], m_att[i][1], m_att[i][2]);
    }
    pm->SetMarkerColor(m_colAttachment);
    pm->SetMarkerSize(m_markerSizeCollision);
    pm->SetBit(kCanDelete);
    pm->Draw("same");
  }

  canvas->Update();
  if (axis) {
    TAxis3D* ax3d = TAxis3D::GetPadAxis();
    ax3d->SetLabelColor(kGray + 2);
    ax3d->SetAxisColor(kGray + 2);
    ax3d->SetXTitle("x");
    ax3d->SetYTitle("y");
    ax3d->SetZTitle("z");
    canvas->Update();
  }
}

bool ViewDrift::SetPlotLimits() {

  if (m_userPlotLimits) return true;
  double xmin = 0., ymin = 0., xmax = 0., ymax = 0.;
  if (m_userBox) {
    if (PlotLimitsFromUserBox(xmin, ymin, xmax, ymax)) {
      m_xMinPlot = xmin;
      m_yMinPlot = ymin;
      m_xMaxPlot = xmax;
      m_yMaxPlot = ymax;
      return true;
    }
  }

  // Try to determine the limits from the drift lines themselves.
  std::array<double, 3> bbmin;
  std::array<double, 3> bbmax;
  bbmin.fill(std::numeric_limits<double>::max());
  bbmax.fill(-std::numeric_limits<double>::max());
  for (const auto& driftLine : m_driftLines) {
    for (const auto& p : driftLine.first) {
      for (unsigned int i = 0; i < 3; ++i) {
        bbmin[i] = std::min(bbmin[i], double(p[i])); 
        bbmax[i] = std::max(bbmax[i], double(p[i]));
      }
    }
  }
  for (const auto& track : m_tracks) {
    for (const auto& p : track) {
      for (unsigned int i = 0; i < 3; ++i) {
        bbmin[i] = std::min(bbmin[i], double(p[i])); 
        bbmax[i] = std::max(bbmax[i], double(p[i]));
      }
    }
  }
  return PlotLimits(bbmin, bbmax, 
                    m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot);
}

}
