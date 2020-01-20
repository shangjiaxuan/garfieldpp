#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentNeBem2d.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  MediumMagboltz gas;
  gas.SetComposition("ar", 100.);
  
  ComponentNeBem2d cmp;
  const double r = 2.;
  const unsigned int n = 6;
  std::vector<double> xv(n, 0.);
  std::vector<double> yv(n, 0.);
  for (unsigned int i = 0; i < n; ++i) {
    const double phi = i * TwoPi / n; 
    xv[i] = r * cos(phi);
    yv[i] = r * sin(phi);
  }
  cmp.AddRegion(xv, yv, &gas, 1, 0.);
  cmp.AddWire(0, 0, 100.e-4, 5000.);

  cmp.SetNumberOfDivisions(100);
  cmp.Initialise();

  TCanvas canvas("c", "", 600, 600);
  ViewField fieldView;
  fieldView.SetCanvas(&canvas);
  fieldView.SetComponent(&cmp);
  fieldView.SetArea(-1.1 * r, -1.1 * r, 1.1 * r, 1.1 * r);
  fieldView.PlotContour();

  app.Run(true);
}
