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
#include "Garfield/Plotting.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheGrid.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"

#define LOG(x) std::cout<<x<<std::endl

using namespace Garfield;

int main(int argc, char * argv[]) {
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    const bool debug = true;
    
    // The geometry of the RPC.
    double gap= 0.2; double thickness = 0.8; double eps = 5.; double voltage = -10000; double sigma =0.5;
    
    ComponentParallelPlate* RPC = new ComponentParallelPlate();
    RPC->Setup(gap,thickness,eps,voltage,sigma);
    
    // Adding a readout structure.
    const std::string label ="NaN";
    RPC->AddPlane(label);
    
    // Setup the gas, but one can also use a gasfile.
    MediumMagboltz gas;
    gas.SetComposition("ar", 80., "co2", 20.);
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    gas.Initialise(true);
    
    // Set the Penning transfer efficiency.
    constexpr double rPenning = 0.56;
    constexpr double lambdaPenning = 0.;
    gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    // Load the ion mobilities.
    const std::string path = getenv("GARFIELD_HOME");
    gas.LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
    
    // Setting the drift medium.
    SolidBox box(0., 0., 0., 100., 100., gap);
    GeometrySimple geo;
    geo.AddSolid(&box, &gas);
    RPC->SetGeometry(&geo);
    RPC->SetMedium(&gas);
    
    // Create the sensor.
    Sensor sensor;
    sensor.AddComponent(RPC);
    sensor.AddElectrode(RPC, label);
    
    // Enable the calculation of the total induced charge on the electrode, as a function of time.
    sensor.EnableInducedCharge();
    // The resistive layer will make the weighting potential time-dependent.
    sensor.EnableDelayedSignal();
    
    // Set the time bins.
    const unsigned int nTimeBins = 100;
    const double tmin =  0.;
    const double tmax = 40;
    const double tstep = (tmax - tmin) / nTimeBins;
    sensor.SetTimeWindow(tmin, tstep, nTimeBins);
    
    // Set the times the delayed signals will be calculated.
    std::vector<double> times;
    
    for(int i=0;i<nTimeBins;i++){
        times.push_back(tmin+tstep/2+i*tstep);
    }
    sensor.SetDelayedSignalTimes(times);
    
    // Create the AvalancheMicroscopic.
    AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);
    aval.UseWeightingPotential();
    aval.EnableSignalCalculation();
    
    // Set time window where the calculations will be done microscopically.
    const double tMaxWindow =2.;
    aval.SetTimeWindow(0.,tMaxWindow);
    
    // Create the AvalancheGrid for grid-based avalanche calculations that are suited for constant drift fields. This class will take over the calculations of the microscopic class after the set time-window.
    AvalancheGrid avalgrid;
    avalgrid.SetSensor(&sensor);
    avalgrid.SetAvalancheMicroscopic(&aval);
    
    // Setting relevant paramiters of the gas. These can be obtained from the Medium class when using a gas file.
    avalgrid.SetElectronTownsend(1.3);
    avalgrid.SetElectronAttachment(.35);
    
    // Preparing the plotting of the induced charge and signal of the electrode readout.
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
    
    // Initial position of an electron in the drift gap,
    const double x0 = 0.;
    const double y0 = 0.;
    const double z0 = gap*0.999; // Starting at the top of the drift gap.
    const double t0 = 0.;
    const double e0 = 1;
    const double dx = 0.;const double dy = 0.;const double dz = 0;
    
    // Setting the timer for the running time of the algorithm.
    std::clock_t start;
       double duration;
    
    start = std::clock();
    
    // Start avalanche.
    aval.AvalancheElectron(x0, y0, z0, t0,e0,dx,dy,dz);
    // Stops calculation after tMaxWindow ns.
    
    // Start grid based avalanche calculations starting from where the microsocpic calculations stoped,
    avalgrid.StartGridAvalanche(0,gap,1000,-0.05,0.05,500);
    
    // Stop timer.
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    LOG("Script: "<<"Electrons have drifted. It took "<<duration<<"s to run.");
    
    // Plot Signals
    if (plotSignal) {
        signalView->Plot(label,true);
        cSignal->Update();
        gSystem->ProcessEvents();
        
        chargeView->Plot(label,false);
        cCharge->Update();
        gSystem->ProcessEvents();
    }
    
    // Export induced current and charge data as an csv file.
    sensor.ExportCharge(label);
    sensor.ExportSignal(label);
    
    app.Run(kTRUE);
}
