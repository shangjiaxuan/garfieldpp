#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TH1D.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/ComponentUser.hh"


#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentGrid.hh"




#include "Garfield/ViewSignal.hh"
#include "Garfield/Plotting.hh"

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"

using namespace Garfield;


int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);
  
    double d = 199.5e-04;
    double pitch = 55e-4;
    // Define the medium.
    MediumSilicon si;
    si.SetTemperature(293.);
    
    ComponentGrid efield;
    efield.LoadElectricField("Efield.txt", "XY", true, false);
    //efield.EnablePeriodicityX();
    efield.SetMedium(&si);
   
	efield.LoadAttachment("Attachment.txt", "XY", 1, 1, 'e');
    efield.LoadAttachment("Attachment.txt", "XY", 1, 2, 'h');
    efield.ActivateTraps();


    ComponentGrid wfield;
    wfield.LoadWeightingField("Wfield.txt", "XY", true);
    wfield.EnablePeriodicityX();
    

    Sensor sensor;
    sensor.AddComponent(&efield);
    const std::string label = "pixel";
    sensor.AddElectrode(&wfield, label);

    // Set the time bins.
    const unsigned int nTimeBins = 1000;
    const double tmin = 0.;
    const double tmax = 10.;
    const double tstep = (tmax - tmin) / nTimeBins;
    sensor.SetTimeWindow(tmin, tstep, nTimeBins);

    // Set up Heed.
    TrackHeed track;
    track.SetSensor(&sensor);
    // Set the particle type and momentum [eV/c].
    track.SetParticle("pion");
    track.SetMomentum(180.e9);

    // Simulate electron/hole drift lines using MC integration.
    AvalancheMC drift;
    drift.SetSensor(&sensor);
    // Use steps of 1 micron.
    drift.SetDistanceSteps(1.e-4);
    drift.EnableSignalCalculation();
    drift.EnableGridTraps();

  
    double x0 = pitch*1.5, y0 = 5.e-5, z0 = 0., t0 = 0.;
    double dx = 0., dy = 1., dz = 0.; 
    track.NewTrack(x0, y0, z0, t0, dx, dy, dz);
    double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
    int ne = 0;
    // Retrieve the clusters along the track.
    while (track.GetCluster(xc, yc, zc, tc, ne, ec, extra)) {
      // Loop over the electrons in the cluster.
      for (int j = 0; j < ne; ++j) {
        double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
        double dxe = 0., dye = 0., dze = 0.;
        track.GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);
        // Simulate the electron and hole drift lines.
        drift.DriftElectron(xe, ye, ze, te);
        drift.DriftHole(xe, ye, ze, te);
      }
    }
    bool plotSignal = true;
    ViewSignal signalView;
    signalView.SetSensor(&sensor);
    constexpr bool plotTotalSignal = true;
    constexpr bool plotElectronSignal = false;
    constexpr bool plotHoleSignal = false;
    signalView.PlotSignal("pixel", plotTotalSignal, plotElectronSignal, plotHoleSignal);
    
    std::ofstream outfile;
    outfile.open("signal.txt", std::ios::out);
    for (unsigned int i = 0; i < nTimeBins; ++i) {
      const double t = (i + 0.5) * tstep;
      const double f = sensor.GetSignal(label, i);
      const double fe = sensor.GetElectronSignal(label, i);
      const double fh = sensor.GetIonSignal(label, i);
      outfile << t << "  " << f << "  " << fe << "  " << fh << "\n";
}
 outfile.close();
    if (plotSignal) {
    app.Run(kTRUE);
  }

}
