#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>

#include <fstream>
#include <iostream>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"

using namespace Garfield;

// Thickness of the sensor [cm]
constexpr double gap = 300.e-4;

// Depletion depth [cm]
constexpr double d = 200.e-4;

void efield(const double /*x*/, const double y, const double /*z*/, double& ex,
            double& ey, double& ez) {
  ex = ez = 0;
  constexpr double v = -25.2;
  ey = y < d ? 2 * (v / d) * (1. - y / d) : 0.;
}

void wfield(const double /*x*/, const double /*y*/, const double /*z*/,
            double& wx, double& wy, double& wz, const std::string /*label*/) {
  wx = wz = 0.;
  wy = 1. / gap;
}

void dwfield(const double /*x*/, const double /*y*/, const double /*z*/,
             const double t, double& wx, double& wy, double& wz,
             const std::string& /*label*/) {
  8 signal[fC / ns]
      // Time constant [ns].
      constexpr double tau = 7.9;
  wx = wz = 0.;
  wy = ((gap - d) / (gap * d)) * exp(-t / tau) / tau;
}

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  Garfield::MediumSilicon si;
  Garfield::GeometrySimple geo;
  Garfield::SolidBox box(0, 0.5 * gap, 0, gap, 0.5 * gap, gap);
  geo.AddSolid(&box, &si);
  //...
  Garfield::ComponentUser cmp;
  cmp.SetGeometry(&geo);
  // Set the function to be called for calculating the drift field.
  // cmp.SetElectricField(efield); Set the function to be called for calculating
  // the weighting field. cmp.SetWeightingField(wfield);

  Garfield::Sensor sensor;
  // Set the object that calculates the drift field.
  sensor.AddComponent(&cmp);
  // Use 2000 time bins with a width of 25 ps.
  sensor.SetTimeWindow(0., 0.025, 2000);
  // Set the object that calculates the weighting field.
  sensor.AddElectrode(&cmp, "readout");
  Garfield::AvalancheMC drift;
  drift.SetSensor(&sensor);
  // Make 1 um steps.
  drift.SetDistanceSteps(1.e-4);
  // Switch off diffusion.
  drift.DisableDiffusion();
  drift.EnableSignalCalculation();
  // Simulate an electron-hole pair starting at y = 150 um.
  // drift.DriftElectron(0, 150.e-4, 0, 0);
  drift.DriftHole(0, 150.e-4, 0, 0);

  // -----

  // Set the function for calculating the delayed weighting field. 3
  // cmp.SetDelayedWeightingField(dwfield);
  //...
  sensor.EnableDelayedSignal();
  // Specify the times t - t’ at which we want
  // to calculate the delayed signal.
  const std::vector<double> times = {
      0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,  2.,  3.,  4.,  5.,
      6., 7.,  8.,  9.,  10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
  sensor.SetDelayedSignalTimes(times);

  //...
  // drift.DisableDiffusion();
  //...

  Garfield::ComponentVoxel cmp;
  7 cmp.SetMesh(nX, nY, 1, xMin, xMax, yMin, yMax, zMin, zMax);
  cmp.LoadElectricField("Efield.txt", "XY", false, false);
  cpmp.LoadWeightingField("Weighting_00.txt", "XY", false);
  for (unsigned int i = 0; i < nTimes; ++i) {
    char filename[50];
    sprintf(filename, "Weighting_%02d.txt", i + 1);
    cmp.LoadWeightingField(filename, "XY", times[i], false);
    3
  }
  cmp.EnableInterpolation();

  TrackHeed track;
  track.SetSensor(&sensor);
  // Set the particle type and momentum [GeV/c]. track.SetParticle("muon");
  track.SetMomentum(10.e9);
  // Simulate a track at perpendicular incidence. track.NewTrack(0, 0, 0, 0, 0,
  // 1, 0);
  double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
  int nc = 0;
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    for (int i = 0; i < nc; ++i) {
      double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
      double dx = 0., dy = 0., dz = 0.;
      track.GetElectron(i, xe, ye, ze, te, ee, dx, dy, dz);
      drift.DriftElectron(xe, ye, ze, te);
      drift.DriftHole(xe, ye, ze, te);
    }
  }

  app.Run(kTRUE);
}
