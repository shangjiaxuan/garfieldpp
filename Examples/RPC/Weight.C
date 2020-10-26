#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <ctime>

#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TH1D.h>

#include "Garfield/ComponentParallelPlate.hh"

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"

#define LOG(x) std::cout<<x<<std::endl

using namespace Garfield;

int main(int argc, char * argv[]) {
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    const bool debug = true;
    
    double gap= 0.5; double thickness = 4; double eps = 8.; double voltage = -10000.; double sigma = 0.01;
    
    ComponentParallelPlate* RPC = new ComponentParallelPlate();
    RPC->Setup(gap,thickness,eps,voltage,sigma);
    
    const std::string label ="NaN";
    RPC->AddPlane(label);
    //RPC->AddStrip(16, 4,label);
    //RPC->AddPixel(30.,0.,10.,10.,label);
    
    
    // Setup the gas.
    MediumMagboltz gas;
    gas.SetComposition("ar", 80., "co2", 20.);
    gas.SetTemperature(293.15);
    gas.SetPressure(200.);
    gas.Initialise(true);
    
    // Set the Penning transfer efficiency.
    constexpr double rPenning = 0.56;
    constexpr double lambdaPenning = 0.;
    gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    // Load the ion mobilities.
    const std::string path = getenv("GARFIELD_HOME");
    gas.LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
    
    SolidBox box(0., 0., 0., 100., 100., gap);
    GeometrySimple geo;
    geo.AddSolid(&box, &gas);
    RPC->SetGeometry(&geo);
    
    
    RPC->SetMedium(&gas);
    
    // Create the sensor.
    Sensor sensor;
    sensor.AddComponent(RPC);
    sensor.AddElectrode(RPC, label);
    sensor.EnableInducedCharge();
  //  sensor.EnableDelayedSignal();
    
    constexpr bool plotField = true;
    if (plotField) {
        ViewField* fieldView = new ViewField();
        fieldView->SetSensor(&sensor);
        fieldView->SetPlane(0, -1, 0, 0, 0, 0);
        fieldView->SetArea(-10, 0, 10, gap);
        fieldView->PlotContour("v");
    }
    // Plot the weighting potential if requested.
    constexpr bool plotWeightingField = false;
    if (plotWeightingField) {
        ViewField* wfieldView = new ViewField();
        wfieldView->SetComponent(RPC);
        wfieldView->SetPlane(0, -1, 0, 0, 0, 0);
        wfieldView->SetArea(-10, 0, 10, gap);
        wfieldView->PlotContourWeightingField(label, "v");
    }
    
    // Set the time bins.
    const unsigned int nTimeBins = 100;
    const double tmin =  0.;
    const double tmax = 100.;
    const double tstep = (tmax - tmin) / nTimeBins;
    sensor.SetTimeWindow(tmin, tstep, nTimeBins);
    
    AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);
    aval.UseWeightingPotential();
    aval.EnableSignalCalculation();

    
    AvalancheMC drift;
    drift.SetSensor(&sensor);
    drift.SetDistanceSteps(2.e-4);
    
    constexpr bool plotSignal = true;
    ViewSignal* signalView = nullptr;
    TCanvas* cSignal = nullptr;
    if (plotSignal) {
        cSignal = new TCanvas("cSignal", "", 600, 600);
        signalView = new ViewSignal();
        signalView->SetCanvas(cSignal);
        signalView->SetSensor(&sensor);
    }
    
    ViewSignal* chargeView = nullptr;
    TCanvas* cCharge = nullptr;
    
    if (plotSignal) {
        cCharge = new TCanvas("cCharge", "", 600, 600);
        chargeView = new ViewSignal();
        chargeView->SetCanvas(cCharge);
        chargeView->SetSensor(&sensor);
    }
    
    constexpr bool plotDrift = false;
    ViewDrift* driftView = nullptr;
    TCanvas* cDrift = nullptr;
    aval.EnablePlotting(driftView);
    drift.EnablePlotting(driftView);
    if (plotDrift) {
        cDrift = new TCanvas("cDrift", "", 600, 600);
        driftView = new ViewDrift();
        driftView->SetArea(-5, -5, 0,5, 5, gap);
        driftView->SetPlane(0, -1, 0, 0, 0, 0);
        driftView->SetCanvas(cDrift);
    }
    
    const double x0 = 0.;
    const double y0 = 0.;
    const double z0 = gap*0.99;
    const double t0 = 0.;
    const double e0 = 0.1;
    const double dx = 0.;const double dy = 0.;const double dz = 0.;
    
    std::clock_t start;
       double duration;
    
    start = std::clock();
    aval.DriftElectron(x0, y0, z0, t0,e0,dx,dy,dz);
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    LOG("Script: "<<"Electrons have drifted. It has taken "<<duration<<"s to run.");
    
    if (plotSignal) {
        signalView->Plot(label);
        cSignal->Update();
        gSystem->ProcessEvents();
        
        chargeView->EnableSignalPlot();
        chargeView->Plot(label);
        cCharge->Update();
        gSystem->ProcessEvents();
    }
    
    
    if (plotDrift) {
        constexpr bool twod = true;
        driftView->Plot(twod);
        cDrift->Update();
        gSystem->ProcessEvents();
    };
    
    app.Run(kTRUE);
}
