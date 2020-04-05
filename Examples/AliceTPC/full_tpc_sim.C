#include <iostream>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  
  // Switch between IROC and OROC.
  constexpr bool iroc = false;

  // Distance between rows of wires [cm].
  constexpr double gap = iroc ? 0.2 : 0.3;

  // Periodicity (wire spacing) [cm].
  constexpr double period = 0.25;

  // Wire diameters [cm].
  // Sense wires.
  constexpr double ds = 0.0020;
  // Cathode wires.
  constexpr double dc = 0.0075;
  // Gate wires.
  constexpr double dg = 0.0075;

  // Voltage settings [V].
  // Sense wires.
  constexpr double vs = iroc ? 1460. : 1570.;
  // Gate wires.
  constexpr double vg = -70.;
  constexpr double deltav = 90.;
 
  // HV plane (drift field).
  constexpr double yHV = 249.7;
  constexpr double vHV = -100000;
 
  // Setup the gas.
  MediumMagboltz gas;
  // Set the temperature [K] and pressure [Torr].
  gas.SetTemperature(293.15);
  gas.SetPressure(750.);
  // Set the gas mixture.
  gas.SetComposition("ne", 85.72, "co2", 9.52, "n2", 4.76);
  // Read the electron transport coefficients from a .gas file.
  gas.LoadGasFile("Ne_90_CO2_10_N2_5_with_mg.gas");
  // Read the ion mobility table from file.
  const std::string garfpath = std::getenv("GARFIELD_HOME");
  gas.LoadIonMobility(garfpath + "/Data/IonMobility_Ne+_Ne.txt");

  // Setup the electric field, using separate components for the 
  // electrons (gating open) and the ions (gating switched on).
  ComponentAnalyticField cmpe;
  ComponentAnalyticField cmpi;
  cmpe.SetMedium(&gas);
  cmpi.SetMedium(&gas);
  cmpe.SetPeriodicityX(period);
  cmpi.SetPeriodicityX(period);

  // Add the sense (anode) wires.
  constexpr double xs = 0.;
  constexpr double ys = gap;
  cmpe.AddWire(xs, ys, ds, vs, "s");
  cmpi.AddWire(xs, ys, ds, vs, "s");
  // Add the cathode wires.
  constexpr double xc = 0.5 * period;
  constexpr double yc = 2 * gap;
  cmpe.AddWire(xc, yc, dc, 0., "c");
  cmpi.AddWire(xc, yc, dc, 0., "c");

  // Add the gate wires.
  constexpr double xg1 = 0.25 * period;
  constexpr double xg2 = 0.75 * period;
  constexpr double yg = 2. * gap + 0.3;
  cmpe.AddWire(xg1, yg, dg, vg, "g", 100., 50., 19.3, 1);
  cmpe.AddWire(xg2, yg, dg, vg, "g", 100., 50., 19.3, 1);
  cmpi.AddWire(xg1, yg, dg, vg + deltav, "g+");
  cmpi.AddWire(xg2, yg, dg, vg - deltav, "g-");

  // Add the planes.
  cmpe.AddPlaneY(0., 0., "pad_plane");
  cmpi.AddPlaneY(0., 0., "pad_plane");
  cmpe.AddPlaneY(yHV, vHV, "HV");
  cmpi.AddPlaneY(yHV, vHV, "HV");

  // Set the magnetic field [T].
  cmpe.SetMagneticField(0, 0.5, 0);
  cmpi.SetMagneticField(0, 0.5, 0);

  cmpi.AddReadout("pad_plane");

  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmpe);
  sensor.AddComponent(&cmpi);
  sensor.AddElectrode(&cmpi, "pad_plane");
  const double xmin = -3 * period; 
  const double xmax =  3 * period;
  sensor.SetArea(xmin, 0., -1., xmax, yHV, 1.);

  // Plot isopotential contours.
  ViewField fieldView;
  fieldView.SetComponent(&cmpi);
  fieldView.SetArea(xmin, 0., xmax, 5 * gap);
  fieldView.SetVoltageRange(-400., 1000.);
  fieldView.PlotContour();

  DriftLineRKF drift;
  drift.SetSensor(&sensor);
  // Polya parameter for gain distribution (for Ne)
  constexpr double theta = 0.4;
  // Average gain.
  constexpr double gain = 10;

  // Set up the charged particle track.
  TrackHeed track;
  track.SetSensor(&sensor);
  track.SetParticle("pi");
  // Set the momentum [eV / c].
  track.SetMomentum(1.e9);

  // Plot the drift lines if requested.
  ViewDrift driftView;
  constexpr bool plotDriftLines = true;
  if (plotDriftLines) drift.EnablePlotting(&driftView);

  // Simulate a track.
  const double xt = xmin;
  const double yt = 0.5 * yHV;
  track.NewTrack(xt, yt, 0., 0., 1., 0., 0.);

  // Retrieve the clusters. 
  double xcl = 0., ycl = 0., zcl = 0., tcl = 0., ecl = 0., extra = 0.;
  int ncl = 0;
  while (track.GetCluster(xcl, ycl, zcl, tcl, ncl, ecl, extra)) {
    // Retrieve the electrons of the cluster.
    for (int i = 0; i < ncl; ++i) {
      double x0, y0, z0, t0, e0, dx0, dy0, dz0;
      track.GetElectron(i, x0, y0, z0, t0, e0, dx0, dy0, dz0);
      // Smear the coordinates to account for the drift up to the ROC.
      constexpr double dT = 0.0198;
      const double sigma = dT * sqrt(std::max(y0, 1.2) - 1.2);
      x0 += RndmGaussian() * sigma;
      z0 += RndmGaussian() * sigma;
      sensor.EnableComponent(0, true);
      sensor.EnableComponent(1, false);
      drift.DriftElectron(x0, 1.2, z0, 0.);
      double x1 = 0., y1 = gap, z1 = 0., t1 = 0.;
      int status = 0;
      drift.GetEndPoint(x1, y1, z1, t1, status);
      sensor.EnableComponent(1, true);
      sensor.EnableComponent(0, false);
      if (y1 > 0.5 * gap && y1 < 1.5 * gap) {
        // Electron drifted to a sensing wire.
        // Sample the gain from a Polya distribution.
        const int nIons = static_cast<int>(std::round(RndmPolya(theta) * gain));
        for (int j = 0; j < nIons; ++j) {
          constexpr double r = 0.01;
          const double angle = RndmGaussian(0, 1.4);
          drift.DriftIon(x1 + r * sin(angle), gap + r * cos(angle), 0, 0);
        }
      }
    }
  }

  // Plot the cell layout.
  ViewCell cellView;
  cellView.SetComponent(&cmpi);
  cellView.SetArea(xmin, 0., xmax, 5 * gap);
  cellView.Plot2d();
  if (plotDriftLines) {
    // Plot the drift lines.
    driftView.SetArea(xmin, 0., xmax, 5 * gap);
    driftView.SetCanvas(cellView.GetCanvas());
    driftView.Plot(true, false);
  }
  app.Run(true);
}
