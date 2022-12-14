#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TCanvas.h>

#include "Garfield/ComponentCST.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Load the field map.
  ComponentCST fm;
  fm.Initialise("drift_std_200_320_2500.efm", "m");
  fm.EnableMirrorPeriodicityX();
  fm.EnableMirrorPeriodicityY();
  fm.PrintRange();

  // Dimensions of the GEM [cm]
  constexpr double pitch = 0.014;

  ViewField fieldView;
  constexpr bool plotField = true;
  if (plotField) {
    fieldView.SetComponent(&fm);
    fieldView.SetPlaneXZ();
    // Set the plot limits in the current viewing plane.
    fieldView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
    TCanvas* cf = new TCanvas("cf", "", 600, 600);
    cf->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf);
    fieldView.PlotContour();
  }

  // Setup the gas.
  MediumMagboltz gas;
  gas.SetComposition("ar", 95., "ch4", 5.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise();

  // Load the ion mobilities.
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
  // Associate the gas with the corresponding field map material. 
  fm.SetGas(&gas);
  fm.PrintMaterials();

  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(&fm);
  sensor.SetArea(-5 * pitch, -5 * pitch, -0.01,
                  5 * pitch,  5 * pitch,  0.025);
  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);

  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.SetDistanceSteps(2.e-4);

  ViewDrift driftView;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    aval.EnablePlotting(&driftView);
    drift.EnablePlotting(&driftView);
  }

  constexpr unsigned int nEvents = 1;
  for (unsigned int i = 0; i < nEvents; ++i) { 
    std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position. 
    const double x0 = -0.5 * pitch + RndmUniform() * pitch;
    const double y0 = 0.;
    const double z0 = 0.02; 
    const double t0 = 0.;
    const double e0 = 0.1;
    aval.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    const unsigned int np = aval.GetNumberOfElectronEndpoints();
    for (unsigned int j = 0; j < np; ++j) {
      double xe1, ye1, ze1, te1, e1;
      double xe2, ye2, ze2, te2, e2;
      int status;
      aval.GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                  xe2, ye2, ze2, te2, e2, status);
      drift.DriftIon(xe1, ye1, ze1, te1);
    }
  }

  if (plotDrift) {
    TCanvas* cd = new TCanvas();
    constexpr bool plotMesh = true;
    if (plotMesh) {
      ViewFEMesh* meshView = new ViewFEMesh();
      meshView->SetArea(-2 * pitch, -2 * pitch, -0.02, 
                         2 * pitch,  2 * pitch, 0.02);
      meshView->SetCanvas(cd);
      meshView->SetComponent(&fm);
      // x-z projection.
      meshView->SetPlane(0, -1, 0, 0, 0, 0);
      meshView->SetFillMesh(true);
      // Set the colors of the metal and the kapton.
      meshView->SetColor(1, kGray);
      meshView->SetColor(2, kYellow + 3);
      meshView->EnableAxes();
      meshView->SetViewDrift(&driftView);
      meshView->Plot();
    } else {
      driftView.SetPlaneXZ();
      driftView.SetArea(-2 * pitch, -0.02, 2 * pitch, 0.02);
      driftView.SetCanvas(cd);
      constexpr bool twod = true;
      driftView.Plot(twod);
    }
  }

  app.Run(true);

}
