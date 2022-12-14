#include <iostream>

#include <TApplication.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackSrim.hh"
#include "Garfield/Random.hh"
#include "Garfield/ViewDrift.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {

  // Application
  TApplication app("app", &argc, argv);

  // Define the medium.
  MediumMagboltz gas;
  gas.SetComposition("ar");
  // Set temperature [K] and pressure [Torr].
  gas.SetPressure(760.0);
  gas.SetTemperature(293.15);

  // Make a component (with constant electric field).
  ComponentConstant cmp;
  // Define the active area and medium.
  constexpr double w = 10.;
  cmp.SetArea(0., -w, -w, w, w, w);
  cmp.SetMedium(&gas);
  cmp.SetElectricField(1000., 0., 0.);
  cmp.SetMagneticField(0., 0., 4.);
  
  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);
  
  // Create a track class.
  TrackSrim tr;
  // Connect the track to a sensor.
  tr.SetSensor(&sensor);
  // Read SRIM output from file.
  const std::string file = "Alpha_in_Ar.txt";
  if (!tr.ReadFile(file)) {
    std::cerr << "Reading SRIM file failed.\n";
    return 0;
  }
  // Set the initial kinetic energy of the particle (in eV).
  tr.SetKineticEnergy(9.e6);

  tr.SetClustersMaximum(200);
  // tr.EnableTransverseStraggling(false);

  // Generate and plot a track.
  ViewDrift viewer;
  tr.EnablePlotting(&viewer);
  tr.EnableDebugging();
  tr.NewTrack(0., 0., 0., 0., 1., 0., 0.);
  viewer.SetArea(0., -w, w, w);
  viewer.Plot(true, true);

  app.Run();
  return 0;
}
