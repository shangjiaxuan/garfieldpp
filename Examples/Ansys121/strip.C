#include <iostream>

#include <TApplication.h>
#include <TCanvas.h>

#include "Garfield/ComponentAnsys121.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Load the field map.
  ComponentAnsys121 fm;
  fm.Initialise("fieldmap/ELIST.lis", "fieldmap/NLIST.lis", 
                "fieldmap/MPLIST.lis", "fieldmap/PRNSOL.lis", "micron");
  fm.EnableMirrorPeriodicityX();
  fm.PrintRange();

  ViewField fieldView;
  fieldView.SetComponent(&fm);
  // Set the plot limits in the current viewing plane.
  fieldView.SetArea(-0.01, 0., 0.01, 0.02);
  fieldView.PlotContour();

  app.Run();
}
