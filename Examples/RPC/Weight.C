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
#include "Garfield/Medium.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheGrid.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"

#include "Garfield/TrackHeed.hh"

#define LOG(x) std::cout<<x<<std::endl

using namespace Garfield;

int main(int argc, char * argv[]) {
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    const bool debug = true;
    
    // The geometry of the RPC.
    double gap= 0.2; double thickness = 0.8; double eps = 5.; double voltage = -10000; double sigma =0.1;
    
    ComponentParallelPlate* RPC = new ComponentParallelPlate();
    RPC->Setup(gap,thickness,eps,voltage,sigma);
    
    // Adding a readout structure.
    const std::string label ="NaN";
    RPC->AddPlane(label);
    
    // Setup the gas, but one can also use a gasfile.
    MediumMagboltz gas;
    gas.LoadGasFile("gasfiles/c2h2f4_93-5_iso_5_sf6_1-5_bis.gas");
    //gas.SetComposition("ar", 80., "co2", 20.);
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
    
    // The resistive layer will make the weighting potential time-dependent.
    sensor.EnableDelayedSignal();
    
    // Set the time bins.
    const unsigned int nTimeBins = 200;
    const double tmin =  0.;
    const double tmax = 50;
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
    const double tMaxWindow =.5;
    aval.SetTimeWindow(0.,tMaxWindow);
    
    // Create the AvalancheGrid for grid-based avalanche calculations that are suited for constant drift fields. This class will take over the calculations of the microscopic class after the set time-window.
    AvalancheGrid avalgrid;
    avalgrid.SetSensor(&sensor);
    avalgrid.SetAvalancheMicroscopic(&aval);
    
    // Setting relevant paramiters of the gas. These can be obtained from the Medium class when using a gas file.
    double alpha =0.; double eta = 0.; double evx =0.;double evy =0.;double evz =0.;
    
    double e[3],v;
    int status;
    Medium* m = nullptr;
    
    sensor.ElectricField(0.,0.,gap/2,e[0],e[1],e[2],v,m,status);
    LOG("Script: "<<"E-field= "<<e[0]<<","<<e[1]<<","<<e[2]<<".");
    m->ElectronTownsend(e[0],e[1],e[2],0.,0.,0.,alpha);
    m->ElectronAttachment(e[0],e[1],e[2],0.,0.,0.,eta);
    m->ElectronVelocity(e[0],e[1],e[2],0.,0.,0.,evx,evy,evz);
    
    LOG("Script::Parameters set as alpha = "<<alpha<<", eta = "<<eta<<",and vz = "<<evz<<".");
    
    avalgrid.SetElectronTownsend(1.33);
    avalgrid.SetElectronAttachment(.35);
    avalgrid.SetElectronVelocity(evz);
    avalgrid.SetGrid(0,gap,1000,-0.05,0.05,500);
    
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
    
    // Set up Heed.
    TrackHeed track;
    track.SetSensor(&sensor);
    // Set the particle type and momentum [eV/c].
    track.SetParticle("pion");
    track.SetMomentum(180.e9);
    
    constexpr bool plotDrift = true;
    ViewDrift* driftView = nullptr;
    TCanvas* cDrift = nullptr;
    if (plotDrift) {
        cDrift = new TCanvas("cDrift", "", 600, 600);
        driftView = new ViewDrift();
        driftView->SetArea(-0.2 , -0.2, 0 , 0.2 , 0.2, gap );
        driftView->SetPlane(0, -1, 0, 0, 0, 0);
        driftView->SetCanvas(cDrift);
        track.EnablePlotting(driftView);
    }
    
    // Setting the timer for the running time of the algorithm.
    std::clock_t start;
    double duration;
    
    start = std::clock();
    
    // Flag to randomise the position of the track.
    constexpr bool smearx = true;
    constexpr unsigned int nEvents = 10;
    // Flag to save the signal to a file.
    constexpr bool writeSignal = true;
    // Simulate a charged-particle track.
    double xt = 0.;
    track.NewTrack(-0.1, 0,1*gap, 0, 0.5, 0, -0.5);
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
            aval.EnablePlotting(driftView);
            if(ze>0.000001 && ze<gap) {aval.AvalancheElectron(xe, ye, ze, te,ee,dxe,dye,dze);
                // Stops calculation after tMaxWindow ns.
                avalgrid.ImportElectronData();
            }
        }
    }
    
    // Start grid based avalanche calculations starting from where the microsocpic calculations stoped,
    avalgrid.StartGridAvalanche();
    
    // Stop timer.
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    LOG("Script: "<<"Electrons have drifted. It took "<<duration<<"s to run.");
    
    if (plotSignal) {
        // Plot signals
        signalView->Plot(label,true);
        cSignal->Update();
        gSystem->ProcessEvents();
        // Plot induced charge
        chargeView->Plot(label,false);
        cCharge->Update();
        gSystem->ProcessEvents();
    }
    
    if (plotDrift) {
        constexpr bool twod = true;
        driftView->Plot(twod);
        cDrift->Update();
        gSystem->ProcessEvents();
    }
    
    // Export induced current and charge data as an csv file.
    sensor.ExportCharge(label);
    sensor.ExportSignal(label);
    
    app.Run(kTRUE);
}
