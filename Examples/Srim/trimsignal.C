#include <iostream>

#include <TApplication.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackTrim.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {

  // Application
  TApplication app("app", &argc, argv);

  // Define the medium.
  MediumSilicon si;

  // Define the geometry.
  ComponentConstant cmp;
  cmp.SetElectricField(0., 5000., 0.);
  // Thickness of the silicon layer [cm]
  constexpr double d = 100.e-4;
  cmp.SetArea(-d, 0., -d, d, d, d);
  cmp.SetMedium(&si); 
  // Define the weighting field and weighting potential
  // (parallel-plate electrode at y = 0).
  cmp.SetWeightingField(0, 1. / d, 0., "readout");
  cmp.SetWeightingPotential(0, 0, 0, 1.);

  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "readout");
  // Set the time bins for the induced current.
  const unsigned int nTimeBins = 1000;
  const double tmin =  0.;
  const double tmax = 10.;
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);
 
  // Read the TRIM output file. 
  TrackTrim tr;
  const std::string filename = "EXYZ.txt";
  // Import the first 100 ions.
  if (!tr.ReadFile(filename, 100)) {
    std::cerr << "Reading TRIM EXYZ file failed.\n";
    return 1;
  }
  // Connect the track to a sensor.
  tr.SetSensor(&sensor);

  DriftLineRKF drift;
  drift.SetSensor(&sensor);
  drift.SetMaximumStepSize(10.e-4);

  // Plot the track and the drift lines.
  ViewDrift driftView;
  tr.EnablePlotting(&driftView);
  drift.EnablePlotting(&driftView);
  
  // Simulate an ion track.
  tr.NewTrack(0., 0., 0., 0., 0., 1., 0.);
  // Loop over the clusters.
  double xc, yc, zc, tc, ec, ekin;
  int ne = 0;
  while (tr.GetCluster(xc, yc, zc, tc, ne, ec, ekin)) {
    // Simulate electron and ion drift lines starting 
    // from the cluster position. 
    // Scale the induced current by the number of electron/ion pairs 
    // in the cluster.
    drift.SetElectronSignalScalingFactor(ne);
    drift.DriftElectron(xc, yc, zc, tc);
    drift.SetHoleSignalScalingFactor(ne);
    drift.DriftHole(xc, yc, zc, tc);
  }
  driftView.SetArea(-2.e-4, 0., 2.e-4, 100.e-4);
  driftView.Plot(true);

  ViewSignal signalView;
  signalView.SetSensor(&sensor);
  signalView.PlotSignal("readout", true, true, true);

  app.Run();
  return 0;
}

