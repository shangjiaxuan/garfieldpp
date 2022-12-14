#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1D.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentAnalyticField.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Define the medium.
  MediumSilicon si;
  si.SetTemperature(300.);
  si.SetLatticeMobilityModelSentaurus();
  si.SetSaturationVelocityModelCanali();
  si.SetHighFieldMobilityModelCanali();

  // Thickness of silicon [cm].
  constexpr double gap = 50.e-4;
  // Length in the x and z direction [cm].
  constexpr double width = 3 * gap;

  // Build the geometry.
  SolidBox box(0., 0.5 * gap, 0., 0.5 * width, 0.5 * gap, 0.5 * width);
  GeometrySimple geo;
  geo.AddSolid(&box, &si);

  // Make a component with analytic weighting field for a pixel.
  constexpr double pitch = 50.e-4;
  ComponentAnalyticField wField;
  wField.SetGeometry(&geo);
  wField.AddPlaneY(0., -1., "back");
  wField.AddPlaneY(gap, 0., "front");
  wField.AddPixelOnPlaneY(gap, -0.5 * pitch, 0.5 * pitch, 
                               -0.5 * pitch, 0.5 * pitch, "pixel");
  wField.AddReadout("pixel");

  const double xmin = 0.;
  const double xmax = width;
  const unsigned int nX = 150;
  const double dx = (xmax - xmin) / nX;
  const double zmin = 0.;
  const double zmax = width;
  const unsigned int nZ = 150;
  const double dz = (zmax - zmin) / nZ;
  const double ymin = 0.;
  const double ymax = gap;
  const unsigned int nY = 50;
  const double dy = (ymax - ymin) / nY;

  const std::string filename = "wfield.txt";
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out);
  for (unsigned int i = 0; i < nX; ++i) {
    const double x = xmin + (i + 0.5) * dx;
    std::cout << i << " (x = " << x << ")\n";
    for (unsigned int j = 0; j < nY; ++j) {
      const double y = ymin + (j + 0.5) * dy;
      for (unsigned int k = 0; k < nZ; ++k) {
        const double z = zmin + (k + 0.5) * dz;
        double wx = 0., wy = 0., wz = 0.;
        wField.WeightingField(x, y, z, wx, wy, wz, "pixel");
        double v = wField.WeightingPotential(x, y, z, "pixel");
        outfile << i << "  " << j << "  " << k << "  " 
                << wx << "  " << wy << "  " << wz << "  " << v << "  0\n"; 
      }
    }
  }
  outfile.close();
}
