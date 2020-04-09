#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TAxis.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  // Make a medium
  MediumMagboltz gas;
  gas.SetComposition("ar");

  // Detector geometry
  // Gap [cm]
  constexpr double width = 1.;
  SolidBox box(0.5 * width, 0, 0, 0.5 * width, 10., 10.);
  GeometrySimple geo;
  geo.AddSolid(&box, &gas);

  // Make a component
  ComponentConstant cmp;
  cmp.SetGeometry(&geo);
  cmp.SetElectricField(100., 0., 0.);

  // Make a sensor
  Sensor sensor;
  sensor.AddComponent(&cmp);

  // Track class
  TrackHeed track;
  track.SetSensor(&sensor);
  track.SetParticle("pi");

  std::vector<double> bg = {0.5, 0.8, 1., 2., 3., 4., 5., 8., 10.,
    12., 15., 20., 50., 100., 200., 500.}; 

  const unsigned int nPoints = bg.size();
  TGraph gStoppingPower(nPoints);
  TGraph gClusterDensity(nPoints);
  for (unsigned int i = 0; i < nPoints; ++i) {
    track.SetBetaGamma(bg[i]);
    // Initial position and direction 
    double x0 = 0., y0 = 0., z0 = 0., t0 = 0.;
    double dx0 = 1., dy0 = 0., dz0 = 0.; 
    track.NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
    const double dedx = track.GetStoppingPower();
    const double imfp = track.GetClusterDensity();
    gStoppingPower.SetPoint(i, bg[i], 1.e-3 * dedx);
    gClusterDensity.SetPoint(i, bg[i], imfp); 
  }

  TCanvas cStoppingPower("cStoppingPower", "", 600, 600);
  cStoppingPower.cd();
  gStoppingPower.SetMarkerStyle(20);
  gStoppingPower.Draw("ap");
  gStoppingPower.GetXaxis()->SetTitle("#beta#gamma");
  gStoppingPower.GetYaxis()->SetTitle("d#it{E}/d#it{x} [keV/cm]");
  cStoppingPower.SetLogx();
  cStoppingPower.SetLogy();
  gStoppingPower.GetYaxis()->SetMoreLogLabels();
  cStoppingPower.Update();

  TCanvas cClusterDensity("cClusterDensity", "", 600, 600);
  cClusterDensity.cd();
  gClusterDensity.SetMarkerStyle(20);
  gClusterDensity.Draw("ap");
  gClusterDensity.GetXaxis()->SetTitle("#beta#gamma");
  gClusterDensity.GetYaxis()->SetTitle("cluster density [cm^{-1}]");
  cClusterDensity.SetLogx();
  cClusterDensity.Update();
  app.Run(true);

}
