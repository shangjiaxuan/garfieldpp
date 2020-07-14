#include <iostream>

#include <TGraph.h>
#include <TAxis.h>

#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewSignal.hh"

namespace Garfield {

ViewSignal::ViewSignal() : ViewBase("ViewSignal") {}

void ViewSignal::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }
  m_sensor = s;
}

void ViewSignal::SetRangeX(const double xmin, const double xmax) {

  if (fabs(xmax - xmin) < Small) {
    std::cerr << m_className << "::SetRangeX: Invalid range.\n";
    return;
  }
  m_xmin = std::min(xmin, xmax);
  m_xmax = std::max(xmin, xmax);
  m_userRangeX = true;
}

void ViewSignal::SetRangeY(const double ymin, const double ymax) {

  if (fabs(ymax - ymin) < Small) {
    std::cerr << m_className << "::SetRangeY: Invalid range.\n";
    return;
  }
  m_ymin = std::min(ymin, ymax);
  m_ymax = std::max(ymin, ymax);
  m_userRangeY = true;
}

void ViewSignal::PlotSignal(const std::string& label, const bool total,
                            const bool electron, const bool ion,
                            const bool delayed) {
  if (!m_sensor) {
    std::cerr << m_className << "::PlotSignal: Sensor is not defined.\n";
    return;
  }

  auto canvas = GetCanvas();
  canvas->cd();
  canvas->SetTitle("Signal");  

  unsigned int nBins = 100;
  double t0 = 0., dt = 1.;
  m_sensor->GetTimeWindow(t0, dt, nBins);
  const double t1 = t0 + nBins * dt;

  std::string xlabel = "time [ns]"; 
  std::string ylabel = m_labelY;
  if (ylabel.empty()) {
    ylabel = m_sensor->IsIntegrated(label) ? "signal [fC]" : "signal [fC / ns]";
  }
  if (total) {
    const auto hname = FindUnusedHistogramName("hSignal_");
    m_hSignal.reset(new TH1D(hname.c_str(), "", nBins, t0, t1));
    m_hSignal->SetLineColor(m_colTotal);
    m_hSignal->GetXaxis()->SetTitle(xlabel.c_str());
    m_hSignal->GetYaxis()->SetTitle(ylabel.c_str());
    m_hSignal->SetStats(0);
    for (unsigned int i = 0; i < nBins; ++i) {
      const double sig = m_sensor->GetSignal(label, i);
      m_hSignal->SetBinContent(i + 1, sig);
    }
    m_hSignal->DrawCopy("");
    if (m_userRangeX) m_hSignal->SetAxisRange(m_xmin, m_xmax, "X");
    if (m_userRangeY) m_hSignal->SetAxisRange(m_ymin, m_ymax, "Y");

    // Get and plot threshold crossings.
    const auto nCrossings = m_sensor->GetNumberOfThresholdCrossings();
    if (nCrossings > 0) {
      TGraph gCrossings;
      gCrossings.SetMarkerStyle(20);
      gCrossings.SetMarkerColor(m_colTotal);
      std::vector<double> xp;
      std::vector<double> yp;
      double time = 0., level = 0.;
      bool rise = true;
      for (unsigned int i = 0; i < nCrossings; ++i) {
        if (m_sensor->GetThresholdCrossing(i, time, level, rise)) {
          xp.push_back(time);
          yp.push_back(level);
        }
      }
      gCrossings.DrawGraph(xp.size(), xp.data(), yp.data(), "psame");
    } 

    if (delayed) {
      const auto hnamed = FindUnusedHistogramName("hDelayedSignal_");
      m_hDelayedSignal.reset(new TH1D(hnamed.c_str(), "", nBins, t0, t1));
      m_hDelayedSignal->SetLineColor(m_colDelayed[0]);
      m_hDelayedSignal->SetLineStyle(2);
      m_hDelayedSignal->SetStats(0);
      for (unsigned int i = 0; i < nBins; ++i) {
        const double sig = m_sensor->GetDelayedElectronSignal(label, i) +
                           m_sensor->GetDelayedIonSignal(label, i);
        m_hDelayedSignal->SetBinContent(i + 1, sig);
      }
      m_hDelayedSignal->DrawCopy("same");
    }
    gPad->Update();
  }

  // Plot the electron and ion signals if requested.
  if (electron) {
    const auto hname = FindUnusedHistogramName("hSignalElectrons_");
    m_hSignalElectrons.reset(new TH1D(hname.c_str(), "", nBins, t0, t1));
    m_hSignalElectrons->SetLineColor(m_colElectrons);
    m_hSignalElectrons->GetXaxis()->SetTitle(xlabel.c_str());
    m_hSignalElectrons->GetYaxis()->SetTitle(ylabel.c_str());
    m_hSignalElectrons->SetStats(0);
    for (unsigned int i = 0; i < nBins; ++i) {
      const double sig = m_sensor->GetElectronSignal(label, i);
      m_hSignalElectrons->SetBinContent(i + 1, sig);
    }
    m_hSignalElectrons->DrawCopy("same");
    if (!total) {
      if (m_userRangeX) m_hSignalElectrons->SetAxisRange(m_xmin, m_xmax, "X");
      if (m_userRangeY) m_hSignalElectrons->SetAxisRange(m_ymin, m_ymax, "Y");
    }
    if (delayed) {
      const auto hnamed = FindUnusedHistogramName("hDelayedSignalElectrons_");
      m_hDelayedSignalElectrons.reset(new TH1D(hnamed.c_str(), "", nBins, t0, t1));
      m_hDelayedSignalElectrons->SetLineColor(m_colDelayed[1]);
      m_hDelayedSignalElectrons->SetLineStyle(2);
      m_hDelayedSignalElectrons->SetStats(0);
      for (unsigned int i = 0; i < nBins; ++i) {
        const double sig = m_sensor->GetDelayedElectronSignal(label, i);
        m_hDelayedSignalElectrons->SetBinContent(i + 1, sig);
      }
      m_hDelayedSignalElectrons->DrawCopy("same");
    }
    gPad->Update();
  }
  if (ion) {
    const auto hname = FindUnusedHistogramName("hSignalIons_");
    m_hSignalIons.reset(new TH1D(hname.c_str(), "", nBins, t0, t1));
    m_hSignalIons->SetLineColor(m_colIons);
    m_hSignalIons->GetXaxis()->SetTitle(xlabel.c_str());
    m_hSignalIons->GetYaxis()->SetTitle(ylabel.c_str());
    m_hSignalIons->SetStats(0);
    for (unsigned int i = 0; i < nBins; ++i) {
      const double sig = m_sensor->GetIonSignal(label, i);
      m_hSignalIons->SetBinContent(i + 1, sig);
    }
    m_hSignalIons->DrawCopy("same");
    if (!(total || electron)) {
      if (m_userRangeX) m_hSignalIons->SetAxisRange(m_xmin, m_xmax, "X");
      if (m_userRangeY) m_hSignalIons->SetAxisRange(m_ymin, m_ymax, "Y");
    }
    if (delayed) {
      const auto hnamed = FindUnusedHistogramName("hDelayedSignalIons_");
      m_hDelayedSignalIons.reset(new TH1D(hnamed.c_str(), "", nBins, t0, t1));
      m_hDelayedSignalIons->SetLineColor(m_colDelayed[2]);
      m_hDelayedSignalIons->SetLineStyle(2);
      m_hDelayedSignalIons->SetStats(0);
      for (unsigned int i = 0; i < nBins; ++i) {
        const double sig = m_sensor->GetDelayedIonSignal(label, i);
        m_hDelayedSignalIons->SetBinContent(i + 1, sig);
      }
      m_hDelayedSignalIons->DrawCopy("same");
    }
    gPad->Update();
  }
}

}
