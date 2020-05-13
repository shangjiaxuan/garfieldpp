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

    // Define the medium.
    MediumSilicon si;
    si.SetTemperature(293.);
    
    ComponentGrid efield;
    efield.LoadElectricField("Efield.txt", "XY", false, false);
    efield.EnablePeriodicityX();
    efield.SetMedium(&si);

    ComponentGrid wfield;
    efield.LoadWeightingField("Efield.txt", "XY", false, false);
    efield.EnablePeriodicityX();
    efield.SetMedium(&si);
    
    ComponentGrid attachment;
    efield.LoadAttachmnet("Efield.txt", "XY", false, false);
    efield.EnablePeriodicityX();
    efield.SetMedium(&si);

    Sensor sensor;
    sensor.AddComponent(&efield);
    const std::string label = "pixel";
    sensor.AddElectrode(&wfield, label);
    sensor.AddComponent(&attachmentfield);

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

    // Plot the signal if requested.
    constexpr bool plotSignal = true;
    ViewSignal* signalView = nullptr;
    TCanvas* cSignal = nullptr;
    if (plotSignal) {
        cSignal = new TCanvas("cSignal", "", 600, 600);
        signalView = new ViewSignal();
        signalView->SetCanvas(cSignal);
        signalView->SetSensor(&sensor);
    }

    constexpr bool plotDrift = true;
    ViewDrift* driftView = nullptr;
    TCanvas* cDrift = nullptr;
    if (plotDrift) {
        cDrift = new TCanvas("cDrift", "", 600, 600);
        driftView = new ViewDrift();
        driftView->SetArea(-0.5 * d, 0, -0.5 * d, 0.5 * d, d, 0.5 * d);
        driftView->SetCanvas(cDrift);
        track.EnablePlotting(driftView);
    }
    // Flag to randomise the position of the track.  
    constexpr bool smearx = true;
    constexpr unsigned int nEvents = 10;
    // Flag to save the signal to a file.
    constexpr bool writeSignal = true;
    for (unsigned int i = 0; i < nEvents; ++i) {
        if (plotDrift) driftView->Clear();
        // Reset the signal.
        sensor.ClearSignal();
        if (i % 10 == 0) std::cout << i << "/" << nEvents << "\n";
        // Simulate a charged-particle track.
        double xt = 0.;
        if (smearx) xt = -0.5 * pitch + RndmUniform() * pitch;
        track.NewTrack(xt, 0, 0, 0, 0, 1, 0);
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
                if (plotDrift) {
                    drift.DisablePlotting();
                    if (RndmUniform() < 0.01) drift.EnablePlotting(driftView);
                }
                drift.DriftElectron(xe, ye, ze, te);
                drift.DriftHole(xe, ye, ze, te);
            }
        }
        if (plotSignal) {
            signalView->PlotSignal(label);
            cSignal->Update();
            gSystem->ProcessEvents();
        }
        if (plotDrift) {
            constexpr bool twod = true;
            driftView->Plot(twod);
            cDrift->Update();
            gSystem->ProcessEvents();
        }
        // Save the induced current signal to a file.
        if (writeSignal) {
            char filename[50];
            sprintf(filename, "signal_%05d.txt", i);
            std::ofstream outfile;
            outfile.open(filename, std::ios::out);
            for (unsigned int j = 0; j < nTimeBins; ++j) {
                const double t = (j + 0.5) * tstep;
                const double f = sensor.GetSignal(label, j);
                const double fe = sensor.GetElectronSignal(label, j);
                const double fh = sensor.GetIonSignal(label, j);
                outfile << t << "  " << f << "  " << fe << "  " << fh << "\n";
            }
            outfile.close();
        }
    }

    if (plotVelocity || plotSignal || plotDrift ||
        plotField || plotWeightingField) {
        app.Run(kTRUE);
    }
}
