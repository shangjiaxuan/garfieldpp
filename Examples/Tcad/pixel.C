#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

// Transfer function
double transfer(double t) {
  constexpr double tR =  5.6;
  constexpr double tI =  1.8;
  constexpr double tA = 47.;
  constexpr double c1 = tA / ((tA - tI) * (tA - tI) * (tA - tR));
  constexpr double c2 = 1. / ((tA - tI) * tI * (tI - tR));
  constexpr double c3 = tR / ((tA - tR) * (tI - tR) * (tI - tR));
  constexpr double c4 = (tI * tI - tA * tR) / 
                        ((tA - tI) * (tA - tI) * (tI - tR) * (tI - tR));
  const double f1 = -exp(-t / tA) * c1;
  const double f2 =  exp(-t / tI) * t * c2; 
  const double f3 =  exp(-t / tR) * c3;
  const double f4 =  exp(-t / tI) * c4; 
  // constexpr double g = 0.07 / 0.46938;
  return tA * tR * (f1 + f2 + f3 + f4); 
}

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Sensor thickness.
  const double d = 100.e-4;
  // Strip pitch.
  const double pitch = 55.e-4;
  const double width = 3 * pitch;

  MediumSilicon si;

  // Import a two-dimensional TCAD field map.
  ComponentTcad2d fm;
  // Load the mesh (.grd file) and electric field (.dat).
  fm.Initialise("pixel_des.grd", "pixel_des.dat");
  fm.SetRangeZ(-width, width);
  // Associate the silicon regions in the field map with a medium object. 
  fm.SetMedium("Silicon", &si);

  ComponentAnalyticField wfield;
  wfield.AddPlaneY(0,    0., "bot");
  wfield.AddPlaneY(d, -100., "top");
  wfield.AddStripOnPlaneY('z', d, 
                          0.5 * width - 0.5 * pitch,
                          0.5 * width + 0.5 * pitch, "strip");
  wfield.AddReadout("strip");

  ViewField vField;
  constexpr bool plotField = true;
  if (plotField) {
    vField.SetComponent(&fm);
    vField.PlotContour("v");
  }

  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&fm);
  sensor.AddElectrode(&wfield, "strip");

  const int nSignalBins = 2000;
  const double tStep = 0.01;
  sensor.SetTimeWindow(0., tStep, nSignalBins);
  sensor.SetTransferFunction(transfer);
  // Threshold.
  const double thr1 = -1000. * ElementaryCharge;  
  std::cout << "Threshold: " << thr1 << " fC\n";

  // Charged-particle track.
  TrackHeed track;
  track.SetSensor(&sensor);
  track.SetParticle("pi");
  track.SetMomentum(180.e9);

  // Electron/hole transport.
  AvalancheMC drift;
  drift.SetDistanceSteps(1.e-4);
  drift.SetSensor(&sensor);
  drift.EnableSignalCalculation();

  ViewSignal vSignal;
  vSignal.SetSensor(&sensor);
  constexpr bool plotSignal = true;

  ViewDrift vDrift;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    vDrift.SetArea(0., 0., width, d);
    track.EnablePlotting(&vDrift);
  }

  const unsigned int nEvents = 1;
  for (unsigned int j = 0; j < nEvents; ++j) {
    sensor.ClearSignal();
    const double x0 = 0.5 * width + (RndmUniform() - 0.5) * pitch;
    const double t0 = 0.1; 
    track.NewTrack(x0, 0., 0., t0, 0., 1., 0.);
    double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., dummy = 0.;
    int ne = 0, nh = 0;
    unsigned int nc = 0;
    unsigned int nesum = 0;
    while (track.GetCluster(xc, yc, zc, tc, ne, nh, ec, dummy)) {
      ++nc;
      nesum += ne;
      if (nc % 100 == 0) std::cout << "    Cluster " << nc << "\n";
      drift.DisablePlotting();
      if (plotDrift && RndmUniform() < 0.05) {
        drift.EnablePlotting(&vDrift);
      }
      double xe, ye, ze, te, ee, dxe, dye, dze;
      for (int i = 0; i < ne; ++i) {
        track.GetElectron(i, xe, ye, ze, te, ee, dxe, dye, dze);
        drift.DriftElectron(xe, ye, ze, te);
      }
      double xh, yh, zh, th;
      for (int i = 0; i < nh; ++i) {
        track.GetIon(i, xh, yh, zh, th);
        drift.DriftHole(xh, yh, zh, th);
      }
    }
    std::cout << nesum << " electrons, " 
              << nesum * ElementaryCharge << " fC.\n";
    // Convolute the signal with the transfer function.
    sensor.ConvoluteSignals();
    // Plot the signal.
    if (!plotSignal) continue;
    vSignal.PlotSignal("strip");
  }

  if (plotDrift) {
    const bool twod = true;
    vDrift.Plot(twod, true);
  }

  app.Run(true);

}
