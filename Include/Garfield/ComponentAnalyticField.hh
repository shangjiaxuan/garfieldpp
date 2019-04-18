#ifndef G_COMPONENT_ANALYTIC_FIELD_H
#define G_COMPONENT_ANALYTIC_FIELD_H

#include <cmath>
#include <complex>

#include "ComponentBase.hh"
#include "FundamentalConstants.hh"

namespace Garfield {

/// Semi-analytic calculation of two-dimensional configurations
/// consisting of wires, planes, and tubes.

class ComponentAnalyticField : public ComponentBase {
 public:
  /// Constructor
  ComponentAnalyticField();
  /// Destructor
  ~ComponentAnalyticField() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override {
    m = nullptr;
    // Calculate the field.
    double v = 0.;
    status = Field(x, y, z, ex, ey, ez, v, false);
    // If the field is ok, get the medium.
    if (status == 0) {
      m = GetMedium(x, y, z);
      if (!m) {
        status = -6;
      } else if (!m->IsDriftable()) {
        status = -5;
      }
    }
  }

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override {
    m = nullptr;
    // Calculate the field.
    status = Field(x, y, z, ex, ey, ez, v, true);
    // If the field is ok, get the medium.
    if (status == 0) {
      m = GetMedium(x, y, z);
      if (!m) {
        status = -6;
      } else if (!m->IsDriftable()) {
        status = -5;
      }
    }
  }

  bool GetVoltageRange(double& pmin, double& pmax) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override {
    wx = wy = wz = 0.;
    double volt = 0.;
    if (!m_sigset) PrepareSignals();
    Wfield(x, y, z, wx, wy, wz, volt, label, false);
  }
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override {
    double wx = 0., wy = 0., wz = 0.;
    double volt = 0.;
    if (!m_sigset) PrepareSignals();
    Wfield(x, y, z, wx, wy, wz, volt, label, true);
    return volt;
  }

  bool GetBoundingBox(double& x0, double& y0, double& z0, double& x1,
                      double& y1, double& z1) override;

  bool IsWireCrossed(const double x0, const double y0, const double z0,
                     const double x1, const double y1, const double z1,
                     double& xc, double& yc, double& zc) override;

  bool IsInTrapRadius(const double q0, const double x0, const double y0,
                      const double z0, double& xw, double& yx,
                      double& rw) override;

  /// Add a wire at (x, y) .
  void AddWire(const double x, const double y, const double diameter,
               const double voltage, const std::string& label,
               const double length = 100., const double tension = 50.,
               const double rho = 19.3, const int ntrap = 5);
  /// Add a tube.
  void AddTube(const double radius, const double voltage, const int nEdges,
               const std::string& label);
  /// Add a plane at constant x.
  void AddPlaneX(const double x, const double voltage,
                 const std::string& label);
  /// Add a plane at constant y.
  void AddPlaneY(const double y, const double voltage,
                 const std::string& label);

  void AddStripOnPlaneX(const char direction, const double x, const double smin,
                        const double smax, const std::string& label,
                        const double gap = -1.);
  void AddStripOnPlaneY(const char direction, const double y, const double smin,
                        const double smax, const std::string& label,
                        const double gap = -1.);
  void AddPixelOnPlaneX(const double x, const double ymin, const double ymax,
                        const double zmin, const double zmax,
                        const std::string& label, const double gap = -1.);
  void AddPixelOnPlaneY(const double y, const double xmin, const double xmax,
                        const double zmin, const double zmax,
                        const std::string& label, const double gap = -1.);

  /// Set the periodic length [cm] in the x-direction.
  void SetPeriodicityX(const double s);
  /// Set the periodic length [cm] in the y-direction.
  void SetPeriodicityY(const double s);
  bool GetPeriodicityX(double& s);
  bool GetPeriodicityY(double& s);

  /// Add a point charge.
  void AddCharge(const double x, const double y, const double z,
                 const double q);
  void ClearCharges();
  void PrintCharges() const;

  /** Return the cell type.
    * Cells are classified according to the number
    * and orientation of planes, the presence of
    * periodicities and the location of the wires
    * as one of the following types:
    *
    * A    non-periodic cells with at most 1 x- and 1 y-plane
    * B1X  x-periodic cells without x-planes and at most 1 y-plane
    * B1Y  y-periodic cells without y-planes and at most 1 x-plane
    * B2X  cells with 2 x-planes and at most 1 y-plane
    * B2Y  cells with 2 y-planes and at most 1 x-plane
    * C1   doubly periodic cells without planes
    * C2X  doubly periodic cells with x-planes
    * C2Y  doubly periodic cells with y-planes
    * C3   double periodic cells with x- and y-planes
    * D1   round tubes without axial periodicity
    * D2   round tubes with axial periodicity
    * D3   polygonal tubes without axial periodicity
    */
  std::string GetCellType() {
    if (!m_cellset) {
      if (CellCheck()) CellType();
    }
    return GetCellType(m_cellType);
  }

  /// Setup the weighting field for a given group of wires or planes.
  void AddReadout(const std::string& label);

  void EnableChargeCheck(const bool on = true) { m_chargeCheck = on; }
  void DisableChargeCheck() { EnableChargeCheck(false); }

  unsigned int GetNumberOfWires() const { return m_nWires; }
  bool GetWire(const unsigned int i, double& x, double& y, double& diameter,
               double& voltage, std::string& label, double& length,
               double& charge, int& ntrap) const;

  unsigned int GetNumberOfPlanesX() const;
  unsigned int GetNumberOfPlanesY() const;
  bool GetPlaneX(const unsigned int i, double& x, double& voltage,
                 std::string& label) const;
  bool GetPlaneY(const unsigned int i, double& y, double& voltage,
                 std::string& label) const;

  bool GetTube(double& r, double& voltage, int& nEdges,
               std::string& label) const;

  enum Cell {
    A00,
    B1X,
    B1Y,
    B2X,
    B2Y,
    C10,
    C2X,
    C2Y,
    C30,
    D10,
    D20,
    D30,
    D40,
    Unknown
  };

 private:
  bool m_chargeCheck = false;

  bool m_cellset = false;
  bool m_sigset = false;

  bool m_polar = false;

  // Cell type.
  Cell m_cellType;

  // Bounding box
  double m_xmin, m_xmax;
  double m_ymin, m_ymax;
  double m_zmin, m_zmax;

  // Voltage range
  double vmin, vmax;

  // Periodicities
  bool m_perx, m_pery;
  double m_sx, m_sy;

  // Signals
  int m_nFourier = 1;
  Cell m_cellTypeFourier = A00;
  bool m_fperx = false;
  bool m_fpery = false;
  int m_mxmin = 0;
  int m_mxmax = 0;
  int m_mymin = 0;
  int m_mymax = 0;
  int m_mfexp = 0;

  std::vector<std::string> m_readout;

  // Wires
  unsigned int m_nWires;
  struct Wire {
    double x, y;       //< Location.
    double d;          //< Diameter.
    double v;          //< Potential.
    double e;          //< Charge.
    std::string type;  //< Label.
    double u;          //< Length.
    int ind;           //< Readout group.
    /// Trap radius. Particle is "trapped" if within nTrap * radius of wire.
    int nTrap;
  };
  std::vector<Wire> m_w;

  // Stretching weight
  std::vector<double> weight;
  // Density
  std::vector<double> dens;
  // Mirror charges for force calculations
  std::vector<double> cnalso;

  // Option for computation of dipole terms
  bool dipole;
  // Dipole angle and amplitude
  std::vector<double> cosph2;
  std::vector<double> sinph2;
  std::vector<double> amp2;

  // Parameters for B2 type cells
  std::vector<double> m_b2sin;
  // Parameters for C type cells
  int m_mode;
  std::complex<double> m_zmult;
  double m_p1, m_p2, m_c1;
  // Parameters for D3 type cells
  // Conformal mapping in polygons
  std::vector<std::complex<double> > wmap;
  double m_kappa;

  // Reference potential
  double m_v0;
  double m_corvta, m_corvtb, m_corvtc;

  // Planes
  // Existence
  bool m_ynplan[4];
  bool m_ynplax, m_ynplay;
  // Coordinates
  double m_coplan[4];
  double m_coplax, m_coplay;
  // Voltages
  double m_vtplan[4];

  struct Strip {
    std::string type;   //< Label.
    int ind;            //< Readout group.
    double smin, smax;  //< Coordinates.
    double gap;         //< Distance to the opposite electrode.
  };

  struct Pixel {
    std::string type;   //< Label.
    int ind;            //< Readout group.
    double smin, smax;  //< Coordinates in x/y.
    double zmin, zmax;  //< Coordinates in z.
    double gap;         //< Distance to the opposite electrode.
  };

  struct Plane {
    std::string type;            //< Label.
    int ind;                     //< Readout group.
    double ewxcor, ewycor;       //< Background weighting fields
    std::vector<Strip> strips1;  //< x/y strips.
    std::vector<Strip> strips2;  //< z strips.
    std::vector<Pixel> pixels;   //< Pixels.
  };
  std::array<Plane, 5> m_planes;

  // Tube
  bool m_tube;
  int m_mtube, m_ntube;
  double m_cotube;
  double m_cotube2;
  double m_vttube;

  // Capacitance matrix
  std::vector<std::vector<double> > m_a;
  // Signal matrix
  std::vector<std::vector<std::complex<double> > > m_sigmat;
  // Induced charges on planes
  std::vector<std::vector<double> > m_qplane;

  // Point charges
  struct Charge3d {
    double x, y, z;  //< Coordinates.
    double e;        //< Charge.
  };
  std::vector<Charge3d> m_ch3d;
  unsigned int m_nTermBessel = 10;
  unsigned int m_nTermPoly = 100;

  // Gravity
  double down[3];

  void UpdatePeriodicity() override;
  void Reset() override { CellInit(); }

  void CellInit();
  bool Prepare();
  bool CellCheck();
  bool CellType();
  std::string GetCellType(const Cell) const;
  bool PrepareStrips();

  bool PrepareSignals();
  bool SetupWireSignals();
  bool SetupPlaneSignals();

  // Calculation of charges
  bool Setup();
  bool SetupA00();
  bool SetupB1X();
  bool SetupB1Y();
  bool SetupB2X();
  bool SetupB2Y();
  bool SetupC10();
  bool SetupC2X();
  bool SetupC2Y();
  bool SetupC30();
  bool SetupD10();
  bool SetupD20();
  bool SetupD30();

  bool IprA00(const int mx, const int my);
  bool IprB2X(const int my);
  bool IprB2Y(const int mx);
  bool IprC2X();
  bool IprC2Y();
  bool IprC30();
  bool IprD10();
  bool IprD30();

  bool SetupDipole() { return true; }

  // Inversion of capacitance matrix
  bool Charge();

  // Evaluation of the electric field
  int Field(const double xin, const double yin, const double zin, double& ex,
            double& ey, double& ez, double& volt, const bool opt);
  void FieldA00(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldB1X(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldB1Y(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldB2X(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldB2Y(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldC10(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldC2X(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldC2Y(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldC30(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldD10(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldD20(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;
  void FieldD30(const double xpos, const double ypos, double& ex, double& ey,
                double& volt, const bool opt) const;

  // Field due to point charges
  void Field3dA00(const double x, const double y, const double z, double& ex,
                  double& ey, double& ez, double& volt);
  void Field3dB2X(const double x, const double y, const double z, double& ex,
                  double& ey, double& ez, double& volt);
  void Field3dB2Y(const double x, const double y, const double z, double& ex,
                  double& ey, double& ez, double& volt);
  void Field3dD10(const double x, const double y, const double z, double& ex,
                  double& ey, double& ez, double& volt);
  // Evaluation of the weighting field
  bool Wfield(const double xpos, const double ypos, const double zpos,
              double& ex, double& ey, double& ez, double& volt,
              const std::string& label, const bool opt) const;
  void WfieldWireA00(const double xpos, const double ypos, double& ex,
                     double& ey, double& volt, const int mx, const int my,
                     const int sw, const bool opt) const;
  void WfieldWireB2X(const double xpos, const double ypos, double& ex,
                     double& ey, double& volt, const int my, const int sw,
                     const bool opt) const;
  void WfieldWireB2Y(const double xpos, const double ypos, double& ex,
                     double& ey, double& volt, const int mx, const int sw,
                     const bool opt) const;
  void WfieldWireC2X(const double xpos, const double ypos, double& ex,
                     double& ey, double& volt, const int sw,
                     const bool opt) const;
  void WfieldWireC2Y(const double xpos, const double ypos, double& ex,
                     double& ey, double& volt, const int sw,
                     const bool opt) const;
  void WfieldWireC30(const double xpos, const double ypos, double& ex,
                     double& ey, double& volt, const int sw,
                     const bool opt) const;
  void WfieldWireD10(const double xpos, const double ypos, double& ex,
                     double& ey, double& volt, const int sw,
                     const bool opt) const;
  void WfieldWireD30(const double xpos, const double ypos, double& ex,
                     double& ey, double& volt, const int sw,
                     const bool opt) const;
  void WfieldPlaneA00(const double xpos, const double ypos, double& ex,
                      double& ey, double& volt, const int mx, const int my,
                      const int iplane, const bool opt) const;
  void WfieldPlaneB2X(const double xpos, const double ypos, double& ex,
                      double& ey, double& volt, const int my, const int iplane,
                      const bool opt) const;
  void WfieldPlaneB2Y(const double xpos, const double ypos, double& ex,
                      double& ey, double& volt, const int mx, const int iplane,
                      const bool opt) const;
  void WfieldPlaneC2X(const double xpos, const double ypos, double& ex,
                      double& ey, double& volt, const int iplane,
                      const bool opt) const;
  void WfieldPlaneC2Y(const double xpos, const double ypos, double& ex,
                      double& ey, double& volt, const int iplane,
                      const bool opt) const;
  void WfieldPlaneC30(const double xpos, const double ypos, double& ex,
                      double& ey, double& volt, const int iplane,
                      const bool opt) const;
  void WfieldPlaneD10(const double xpos, const double ypos, double& ex,
                      double& ey, double& volt, const int iplane,
                      const bool opt) const;
  void WfieldPlaneD30(const double xpos, const double ypos, double& ex,
                      double& ey, double& volt, const int iplane,
                      const bool opt) const;
  void WfieldStripZ(const double xpos, const double ypos, double& ex,
                    double& ey, double& volt, const int ip, const Strip& strip,
                    const bool opt) const;
  void WfieldStripXy(const double xpos, const double ypos, const double zpos,
                     double& ex, double& ey, double& ez, double& volt,
                     const int ip, const Strip& strip, const bool opt) const;
  void WfieldPixel(const double xpos, const double ypos, const double zpos,
                   double& ex, double& ey, double& ez, double& volt,
                   const int ip, const Pixel& pixel, const bool opt) const;

  // Auxiliary functions for C type cells
  double Ph2(const double xpos, const double ypos) const;
  double Ph2Lim(const double radius) const {
    return -log(abs(m_zmult) * radius * (1. - 3. * m_p1 + 5. * m_p2));
  }
  void E2Sum(const double xpos, const double ypos, double& ex,
             double& ey) const;

  // Mapping function for D30 type cells
  void ConformalMap(const std::complex<double>& z, std::complex<double>& ww,
                    std::complex<double>& wd) const;

  bool InTube(const double x0, const double y0, const double a,
              const int n) const;

  // Transformation between cartesian and polar coordinates
  void Cartesian2Polar(const double x0, const double y0, double& r,
                       double& theta) {
    if (x0 == 0. && y0 == 0.) {
      r = theta = 0.;
      return;
    }
    r = sqrt(x0 * x0 + y0 * y0);
    theta = atan2(y0, x0) * RadToDegree;
  }

  void Polar2Cartesian(const double r, const double theta, double& x0,
                       double& y0) const {
    const double thetap = theta * DegreeToRad;
    x0 = r * cos(thetap);
    y0 = r * sin(thetap);
  }

  // Transformation (r, theta) to (rho, phi) via the map
  // (r, theta) = (exp(rho), 180 * phi / Pi).
  void RTheta2RhoPhi(const double rho, const double phi, double& r,
                     double& theta) const {
    r = exp(rho);
    theta = RadToDegree * phi;
  }
};
}

#endif
