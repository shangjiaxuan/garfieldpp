#include <iostream>

#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Setup the gas.
  MediumMagboltz gas;
  gas.SetComposition("ne", 85.72, "co2", 9.52, "n2", 4.76);

  // Describe the cell layout.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  cmp.SetPolarCoordinates();
  cmp.AddPlaneR(2., 0., "o");
  cmp.AddPlanePhi(0, 0., "r");
  cmp.AddPlanePhi(60., 0., "l");
  cmp.AddWire(1.5, 30., 50.e-4, 500., "w");
  cmp.PrintCell();

  // Plot the potential.
  ViewField fieldView;
  fieldView.SetComponent(&cmp);
  const double xmin = -2.1;
  const double xmax =  2.1;
  const double ymin = -2.1;
  const double ymax =  2.1;
  fieldView.SetArea(xmin, ymin, xmax, ymax);
  fieldView.PlotContour();
  ViewCell cellView;
  cellView.SetCanvas(fieldView.GetCanvas());
  cellView.SetComponent(&cmp);
  cellView.SetArea(xmin, ymin, -1., xmax, ymax, 1.);
  cellView.Plot2d();

  app.Run(true);
}
