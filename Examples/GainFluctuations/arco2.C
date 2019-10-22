#include <iostream>
#include <cmath>

#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"

using namespace Garfield;

int main() {

  TFile outfile("arco2.root", "RECREATE");

  // Electric field [kV / cm].
  double field = 20.; 
  for (unsigned int i = 0; i < 8; ++i) {
    // Initial electron energy [eV].
    const double e0 = 1.;

    // Make a gas medium.
    MediumMagboltz gas;
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    gas.SetComposition("ar", 90., "co2", 10.);
    gas.SetMaxElectronEnergy(150.);
    constexpr double scale = 1.;
    gas.SetExcitationScaling(scale, "ar");
    gas.Initialise();
  
    // Make a drift volume.
    const double gap = 1. / (3.5 * field - 60.);
    SolidBox box(0, 0, gap, 2, 2, gap);
    GeometrySimple geo;
    geo.AddSolid(&box, &gas);
  
    // Make a component with constant drift field.
    ComponentConstant cmp;
    cmp.SetGeometry(&geo);
    cmp.SetElectricField(0, 0, field * 1.e3);

    // Make a sensor.
    Sensor sensor;
    sensor.AddComponent(&cmp);
  
    // Microscopic tracking.
    AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);
    aval.SetCollisionSteps(100000);
    
    // Histograms
    TH1::StatOverflows(true);
    TH1F hElectrons("hElectrons", "Avalanche Size", 200, 0,2000);
    TH1F hIons("hIons", "Avalanche Size", 200, 0, 2000);

    std::cout << field << " kV/cm\n";
    constexpr unsigned int nEvents = 100;
    for (unsigned int j = 0; j < nEvents; ++j) {
      aval.AvalancheElectron(0, 0, gap, 0, e0, 0, 0, 0);
      int ne = 0, ni = 0;
      aval.GetAvalancheSize(ne, ni);
      if (j % 100 == 0) {
        std::cout << j << "/" << nEvents << "\n"
                  << "    " << ne << " electrons\n";
      }
      hElectrons.Fill(ne);
      hIons.Fill(ni);
    }

    outfile.cd();
    std::string dir = std::to_string(int(field));
    outfile.mkdir(dir.c_str());
    outfile.cd(dir.c_str());
    hElectrons.Write();
    hIons.Write();

    const double mean = hElectrons.GetMean();
    const double rms = hElectrons.GetRMS();
    std::cout << " f = " << rms * rms / (mean * mean) << "\n";
    field += 5.;
  }
  outfile.Close();
}
