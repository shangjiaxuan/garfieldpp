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

  Medium* GetMedium(const double x, const double y, const double z) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override {
    m = nullptr;
    // Calculate the field.
    double v = 0.;
    status = Field(x, y, z, ex, ey, ez, v, false);
    // If the field is ok, get the medium.
    if (status == 0) {
      m = m_geometry ? m_geometry->GetMedium(x, y, z) : m_medium;
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
      m = m_geometry ? m_geometry->GetMedium(x, y, z) : m_medium;
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
                     double& xc, double& yc, double& zc, const bool centre,
                     double& rc) override;

  bool IsInTrapRadius(const double q0, const double x0, const double y0,
                      const double z0, double& xw, double& yx,
                      double& rw) override;


  /// Set the medium inside the cell.
  void SetMedium(Medium* medium) { m_medium = medium; }
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
  /// Add a plane at constant radius.
  void AddPlaneR(const double r, const double voltage,
                 const std::string& label);
  /// Add a plane at constant phi.
  void AddPlanePhi(const double phi, const double voltage,
                   const std::string& label);
  /// Add a strip in the y or z direction on an existing plane at constant x.
  void AddStripOnPlaneX(const char direction, const double x, const double smin,
                        const double smax, const std::string& label,
                        const double gap = -1.);
  /// Add a strip in the x or z direction on an existing plane at constant y.
  void AddStripOnPlaneY(const char direction, const double y, const double smin,
                        const double smax, const std::string& label,
                        const double gap = -1.);
  /// Add a strip in the phi or z direction on an existing plane at constant radius.
  void AddStripOnPlaneR(const char direction, const double r, const double smin,
                        const double smax, const std::string& label,
                        const double gap = -1.);
  /// Add a strip in the r or z direction on an existing plane at constant phi.
  void AddStripOnPlanePhi(const char direction, const double phi, const double smin,
                          const double smax, const std::string& label,
                          const double gap = -1.);
  /// Add a pixel on an existing plane at constant x.
  void AddPixelOnPlaneX(const double x, const double ymin, const double ymax,
                        const double zmin, const double zmax,
                        const std::string& label, const double gap = -1.);
  /// Add a pixel on an existing plane at constant y.
  void AddPixelOnPlaneY(const double y, const double xmin, const double xmax,
                        const double zmin, const double zmax,
                        const std::string& label, const double gap = -1.);
  /// Add a pixel on an existing plane at constant radius.
  void AddPixelOnPlaneR(const double r, 
                        const double phimin, const double phimax,
                        const double zmin, const double zmax,
                        const std::string& label, const double gap = -1.);
  /// Add a pixel on an existing plane at constant phi.
  void AddPixelOnPlanePhi(const double phi, 
                          const double rmin, const double rmax,
                          const double zmin, const double zmax,
                          const std::string& label, const double gap = -1.);

  /// Set the periodic length [cm] in the x-direction.
  void SetPeriodicityX(const double s);
  /// Set the periodic length [cm] in the y-direction.
  void SetPeriodicityY(const double s);
  /// Set the periodicity [degree] in phi.
  void SetPeriodicityPhi(const double phi);
  /// Get the periodic length in the x-direction.
  bool GetPeriodicityX(double& s);
  /// Get the periodic length in the y-direction.
  bool GetPeriodicityY(double& s);
  /// Get the periodicity [degree] in phi.
  bool GetPeriodicityPhi(double& s);

  /// Use Cartesian coordinates (default).
  void SetCartesianCoordinates();
  /// Use polar coordinates.
  void SetPolarCoordinates();
  /// Are polar coordinates being used?  
  bool IsPolar() const { return m_polar; }

  /// Print all available information on the cell.
  void PrintCell();

  /// Add a point charge.
  void AddCharge(const double x, const double y, const double z,
                 const double q);
  /// Remove all point charges.
  void ClearCharges();
  /// Print a list of the point charges.
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

  /** Calculate multipole moments for a given wire.
    * \param iw Index of the wire.
    * \param order Order of the highest multipole moment.
    * \param print Print information about the fitting process.
    * \param plot Plot the potential and multipole fit around the wire. 
    * \param rmult Distance in multiples of the wire radius
    *              at which the decomposition is to be carried out.
    * \param eps Used in the fit for calculating the covariance matrix.
    * \param nMaxIter Maximum number of iterations in the fit.
    **/
  bool MultipoleMoments(const unsigned int iw, const unsigned int order = 4,
                        const bool print = false, const bool plot = false,
                        const double rmult = 1., const double eps = 1.e-4,
                        const unsigned int nMaxIter = 20); 
  /// Request dipole terms be included for each of the wires (default: off).
  void EnableDipoleTerms(const bool on = true);
  /// Check the quality of the capacitance matrix inversion (default: off).
  void EnableChargeCheck(const bool on = true) { m_chargeCheck = on; }

  /// Get the number of wires.
  unsigned int GetNumberOfWires() const { return m_nWires; }
  /// Retrieve the parameters of a wire.
  bool GetWire(const unsigned int i, double& x, double& y, double& diameter,
               double& voltage, std::string& label, double& length,
               double& charge, int& ntrap) const;

  /// Get the number of equipotential planes at constant x.
  unsigned int GetNumberOfPlanesX() const;
  /// Get the number of equipotential planes at constant y.
  unsigned int GetNumberOfPlanesY() const;
  /// Get the number of equipotential planes at constant radius.
  unsigned int GetNumberOfPlanesR() const;
  /// Get the number of equipotential planes at constant phi.
  unsigned int GetNumberOfPlanesPhi() const;
  /// Retrieve the parameters of a plane at constant x.
  bool GetPlaneX(const unsigned int i, double& x, double& voltage,
                 std::string& label) const;
  /// Retrieve the parameters of a plane at constant y.
  bool GetPlaneY(const unsigned int i, double& y, double& voltage,
                 std::string& label) const;
  /// Retrieve the parameters of a plane at constant radius.
  bool GetPlaneR(const unsigned int i, double& r, double& voltage,
                 std::string& label) const;
  /// Retrieve the parameters of a plane at constant phi.
  bool GetPlanePhi(const unsigned int i, double& phi, double& voltage,
                   std::string& label) const;
  /// Retrieve the tube parameters.
  bool GetTube(double& r, double& voltage, int& nEdges,
               std::string& label) const;

  /// Calculate the electric field at a given wire position, as if the wire
  /// itself were not there, but with the presence of its mirror images.
  bool ElectricFieldAtWire(const unsigned int iw, double& ex, double& ey);

  /// Set the number of grid lines at which the electrostatic force
  /// as function of the wire displacement is computed.
  void SetScanningGrid(const unsigned int nX, const unsigned int nY);
  /// Specify explicitly the boundaries of the the scanning area (i. e. the
  /// area in which the electrostatic force acting on a wire is computed).
  void SetScanningArea(const double xmin, const double xmax, const double ymin,
                       const double ymax);
  /// Set the scanning area to the largest area around each wire
  /// which is free from other cell elements.
  void SetScanningAreaLargest() { m_scanRange = ScanningRange::Largest; }
  /// Set the scanning area based on the zeroth-order estimates of the
  /// wire shift, enlarged by a scaling factor. This is the default behaviour.
  void SetScanningAreaFirstOrder(const double scale = 2.);
  /// Switch on/off extrapolation of electrostatic forces beyond the 
  /// scanning area (default: off).
  void EnableExtrapolation(const bool on = true) { m_extrapolateForces = on; }

  /// Set the gravity orientation.
  void SetGravity(const double dx, const double dy, const double dz);
  /// Get the gravity orientation.
  void GetGravity(double& dx, double& dy, double& dz) const;

  /** Calculate a table of the forces acting on a wire.
    * \param iw index of the wire
    * \param xMap coordinates of the grid lines in x
    * \param yMap coordinates of the grid lines in y
    * \param fxMap x-components of the force at the grid points
    * \param fyMap y-components of the force at the grid points 
    **/ 
  bool ForcesOnWire(const unsigned int iw, std::vector<double>& xMap,
                    std::vector<double>& yMap,
                    std::vector<std::vector<double> >& fxMap,
                    std::vector<std::vector<double> >& fyMap);
  /** Compute the sag profile of a wire.
    * \param iw index of the wire
    * \param detailed flag to request a detailed calculation of the profile
                      or only a fast one.
    * \param csag coordinate along the wire.
    * \param xsag x components of the sag profile.
    * \param ysag y components of the sag profile.
    * \param stretch relative elongation.
    * \param print flag to print the calculation results or not.
    **/ 
  bool WireDisplacement(const unsigned int iw, const bool detailed,
                        std::vector<double>& csag, std::vector<double>& xsag,
                        std::vector<double>& ysag, double& stretch,
                        const bool print = true);
  /// Set the number of shots used for numerically solving the wire sag
  /// differential equation.
  void SetNumberOfShots(const unsigned int n) { m_nShots = n; }
  /// Set the number of integration steps within each shot (must be >= 1).
  void SetNumberOfSteps(const unsigned int n);

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
  Medium* m_medium = nullptr;

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
  double m_vmin, m_vmax;

  // Periodicities
  bool m_perx = false;
  bool m_pery = false;
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
    double x, y;       ///< Location.
    double r;          ///< Radius.
    double v;          ///< Potential.
    double e;          ///< Charge.
    std::string type;  ///< Label.
    double u;          ///< Length.
    int ind;           ///< Readout group.
    /// Trap radius. Particle is "trapped" if within nTrap * radius of wire.
    int nTrap;
    double tension;    ///< Stretching weight.
    double density;    ///< Density.
  };
  std::vector<Wire> m_w;

  // Option for computation of dipole terms
  bool m_dipole;
  // Dipole angle and amplitude
  std::vector<double> m_cosph2;
  std::vector<double> m_sinph2;
  std::vector<double> m_amp2;

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
    std::string type;   ///< Label.
    int ind;            ///< Readout group.
    double smin, smax;  ///< Coordinates.
    double gap;         ///< Distance to the opposite electrode.
  };

  struct Pixel {
    std::string type;   ///< Label.
    int ind;            ///< Readout group.
    double smin, smax;  ///< Coordinates in x/y.
    double zmin, zmax;  ///< Coordinates in z.
    double gap;         ///< Distance to the opposite electrode.
  };

  struct Plane {
    std::string type;            ///< Label.
    int ind;                     ///< Readout group.
    double ewxcor, ewycor;       ///< Background weighting fields
    std::vector<Strip> strips1;  ///< x/y strips.
    std::vector<Strip> strips2;  ///< z strips.
    std::vector<Pixel> pixels;   ///< Pixels.
  };
  std::array<Plane, 5> m_planes;

  // Tube
  bool m_tube = false;
  int m_mtube = 0;
  int m_ntube = 1;
  double m_cotube = 1.;
  double m_cotube2 = 1.;
  double m_vttube = 0.;

  // Capacitance matrix
  std::vector<std::vector<double> > m_a;
  // Signal matrix
  std::vector<std::vector<std::complex<double> > > m_sigmat;
  // Induced charges on planes
  std::vector<std::vector<double> > m_qplane;

  // Point charges
  struct Charge3d {
    double x, y, z;  ///< Coordinates.
    double e;        ///< Charge.
  };
  std::vector<Charge3d> m_ch3d;
  unsigned int m_nTermBessel = 10;
  unsigned int m_nTermPoly = 100;

  bool m_useElectrostaticForce = true;
  bool m_useGravitationalForce = true;
  // Gravity
  std::array<double, 3> m_down{{0, 0, 1}};
  // Number of shots used for solving the wire sag differential equations
  unsigned int m_nShots = 2;
  // Number of integration steps within each shot.
  unsigned int m_nSteps = 20;
  // Options for setting the range of wire shifts
  // for which the forces are computed.
  enum class ScanningRange { Largest = 0, FirstOrder, User };
  ScanningRange m_scanRange = ScanningRange::FirstOrder;
  // Limits of the user-specified scanning range.
  double m_xScanMin = 0.;
  double m_xScanMax = 0.;
  double m_yScanMin = 0.;
  double m_yScanMax = 0.;
  // Scaling factor for first-order estimate of the scanning range.
  double m_scaleRange = 2.;
  // Number of grid lines at which the forces are stored.
  unsigned int m_nScanX = 11;
  unsigned int m_nScanY = 11;
  // Extrapolate beyond the scanning range or not.
  bool m_extrapolateForces = false;

  void UpdatePeriodicity() override;
  void Reset() override { 
    CellInit(); 
    m_medium = nullptr;
  }

  void CellInit();
  bool Prepare();
  bool CellCheck();
  bool WireCheck() const;
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

  bool SetupDipoleTerms();

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

  // Functions for calculating the electric field at a given wire position,
  // as if the wire itself were not there but with the presence
  // of its mirror images.
  void FieldAtWireA00(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireB1X(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireB1Y(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireB2X(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireB2Y(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireC10(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireC2X(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireC2Y(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireC30(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireD10(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireD20(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
  void FieldAtWireD30(const double xpos, const double ypos, double& ex,
                      double& ey, const std::vector<bool>& cnalso) const;
 
  void DipoleFieldA00(const double xpos, const double ypos, double& ex, 
                      double& ey, double& volt, const bool opt) const; 
  void DipoleFieldB1X(const double xpos, const double ypos, double& ex, 
                      double& ey, double& volt, const bool opt) const; 
  void DipoleFieldB1Y(const double xpos, const double ypos, double& ex, 
                      double& ey, double& volt, const bool opt) const; 
  void DipoleFieldB2X(const double xpos, const double ypos, double& ex, 
                      double& ey, double& volt, const bool opt) const; 
  void DipoleFieldB2Y(const double xpos, const double ypos, double& ex, 
                      double& ey, double& volt, const bool opt) const; 

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

  bool SagDetailed(const Wire& wire, const std::vector<double>& xMap,
                   const std::vector<double>& yMap,
                   const std::vector<std::vector<double> >& fxMap,
                   const std::vector<std::vector<double> >& fyMap,
                   std::vector<double>& csag, std::vector<double>& xsag,
                   std::vector<double>& ysag) const;
  bool GetForceRatio(const Wire& wire, const double coor,
                     const std::array<double, 2>& bend,
                     const std::array<double, 2>& dbend,
                     std::array<double, 2>& f, const std::vector<double>& xMap,
                     const std::vector<double>& yMap,
                     const std::vector<std::vector<double> >& fxMap,
                     const std::vector<std::vector<double> >& fyMap) const;
  bool FindZeroes(const Wire& wire, const double h, std::vector<double>& x,
                  const std::vector<double>& xMap,
                  const std::vector<double>& yMap,
                  const std::vector<std::vector<double> >& fxMap,
                  const std::vector<std::vector<double> >& fyMap) const;
  bool StepRKN(const Wire& wire, const double h, double& x,
               std::array<double, 2>& y, std::array<double, 2>& yp,
               const std::vector<double>& xMap, const std::vector<double>& yMap,
               const std::vector<std::vector<double> >& fxMap,
               const std::vector<std::vector<double> >& fyMap) const;
  bool Trace(const Wire& wire, const double h, const std::vector<double>& xx,
             std::vector<double>& f, const std::vector<double>& xMap,
             const std::vector<double>& yMap,
             const std::vector<std::vector<double> >& fxMap,
             const std::vector<std::vector<double> >& fyMap) const;
};
}  // namespace Garfield

#endif
