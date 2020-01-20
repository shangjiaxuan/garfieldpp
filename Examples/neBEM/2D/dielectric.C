#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumPlastic.hh"
#include "Garfield/ComponentNeBem2d.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  MediumMagboltz gas;
  gas.SetComposition("ar", 100.);

  MediumPlastic plastic;
  plastic.SetDielectricConstant(5.);
 
  ComponentNeBem2d cmp;
  cmp.SetMedium(&gas);
  cmp.SetNumberOfDivisions(10);
  constexpr double delta = 1.;
  constexpr double v = 100.;
  // Left conducting plate.
  const double xMin = -1.5 * delta;
  cmp.AddSegment(xMin, -11., xMin, -1., v);
  cmp.AddSegment(xMin,  -1., xMin,  1., v);
  cmp.AddSegment(xMin,   1., xMin, 11., v);

  // Right conducting plate.
  const double xMax = 1.5 * delta;
  cmp.AddSegment(xMax, -11., xMax, -1., -v);
  cmp.AddSegment(xMax,  -1., xMax,  1., -v);
  cmp.AddSegment(xMax,   1., xMax, 11., -v);

  // Dielectric.
  const double xD = 0.5 * delta;
  std::vector<double> xv = {-xD, -xD, xD, xD}; 
  std::vector<double> yv = {-11., 11., 11., -11.};
  cmp.AddRegion(xv, yv, &plastic);

  // cmp.EnableDebugging();
  cmp.Initialise();

  const double eps2 = plastic.GetDielectricConstant();
  const double f1 = (2. * v / delta) / (2. + 1. / eps2);
  const double f2 = f1 / eps2;

  std::ofstream outfile;
  outfile.open("test.txt", std::ios::out);
  const int nSteps = 100;
  const double dx = 3. * delta / nSteps;
  for (int i = 1; i < nSteps; ++i) {
    const double x = xMin + i * dx;
    Medium* medium = nullptr;
    double ex = 0., ey = 0., ez = 0., p = 0.; 
    int stat = 0;
    cmp.ElectricField(x, 0., 0., ex, ey, ez, p, medium, stat);
    const double exact = x <= -xD ? f1 : x < xD ? f2 : f1;
    outfile << x << "  " << p << "  " 
            << std::setprecision(10) << ex << "  " << ey << "  " 
            << std::setprecision(10) << exact << std::endl;
  }
  outfile.close();

}
