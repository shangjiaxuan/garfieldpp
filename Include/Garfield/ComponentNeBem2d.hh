#ifndef G_COMPONENT_NEBEM_2D_H
#define G_COMPONENT_NEBEM_2D_H

#include "ComponentBase.hh"

namespace Garfield {

/// Two-dimensional implementation of the nearly exact Boundary Element Method.

class ComponentNeBem2d : public ComponentBase {
 public:
  /// Constructor
  ComponentNeBem2d();
  /// Destructor
  ~ComponentNeBem2d() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  bool GetVoltageRange(double& vmin, double& vmax) override;

  Medium* GetMedium(const double x, const double y, const double z) override;

  /** Add a conducting straight-line segment.
    * \param x0,y0,x1,y1 coordinates of start and end point.
    * \param v applied potential.
    */
  bool AddSegment(const double x0, const double y0, const double x1,
                  const double y1, const double v, const int ndiv = -1);
  /** Add a wire.
    * \param x,y centre of the wire.
    * \param d wire diameter.
    * \param v applied potential.
    */
  bool AddWire(const double x, const double y, const double d, const double v);
  /** Add a region bounded by a polygon.
    * \param xp,yp x/y-coordinates of the vertices of the polygon.
    * \param medium pointer to the medium associated to the region. 
    * \param bctype 1: fixed voltage, 4: dielectric-dielectric interface.
    * \param v applied potential.
    */ 
  bool AddRegion(const std::vector<double>& xp,
                 const std::vector<double>& yp, Medium* medium,
                 const unsigned int bctype = 4, const double v = 0.,
                 const int ndiv = -1);
 
  /// Discretise the geometry and compute the solution.
  bool Initialise();

  /// Set the default number of elements per segment.
  void SetNumberOfDivisions(const unsigned int ndiv);

  void SetNumberOfCollocationPoints(const unsigned int ncoll);
  void EnableAutoResizing(const bool on = true) { m_autoSize = on; }
  void EnableRandomCollocation(const bool on = true) {
    m_randomCollocation = on;
  }
  void SetMaxNumberOfIterations(const unsigned int niter);

  /// Return the number of conducting straight-line segments.
  unsigned int GetNumberOfSegments() const { return m_segments.size(); }
  /// Return the coordinates and voltage of a given straight-line segment.
  bool GetSegment(const unsigned int i, double& x0, double& y0, 
                  double& x1, double& x2, double& v) const;
  /// Return the number of wires.
  unsigned int GetNumberOfWires() const { return m_wires.size(); }
  /// Return the coordinates, diameter, potential and charge of a given wire.
  bool GetWire(const unsigned int i, double& x, double& y, double& d,
               double& v, double& q) const; 
  /// Return the number of boundary elements.
  unsigned int GetNumberOfElements() const { return m_elements.size(); }
  /// Return the coordinates and charge of a given boundary element.
  bool GetElement(const unsigned int i, double& x0, double& y0,
                  double& x1, double& y1, double& q) const;
 private:
  static const double InvEpsilon0;
  static const double InvTwoPiEpsilon0;

  unsigned int m_nDivisions = 5;
  unsigned int m_nCollocationPoints = 1;
  bool m_autoSize = false;
  bool m_randomCollocation = false;
  unsigned int m_nMaxIterations = 3;

  // Boundary condition types
  enum BC {
    Voltage = 1,
    Charge,
    Floating,
    Dielectric
  };
  struct Region {
    std::vector<double> xv;    //< x-coordinates of the vertices.
    std::vector<double> yv;    //< y-coordinates of the vertices.
    Medium* medium;            //< Medium associated to the region.
    std::pair<BC, double> bc;  //< Applied boundary condition.
    unsigned int depth;        //< Level in the hierarchy.
    int ndiv;                  //< Number of elements per segment.
  };
  std::vector<Region> m_regions;

  struct Segment {
    std::array<double, 2> x0;  //< Coordinates of the start point.
    std::array<double, 2> x1;  //< Coordinates of the end point.
    int region1;               //< Inner region. 
    int region2;               //< Outer region.
    std::pair<BC, double> bc;  //< Applied boundary condition.
    int ndiv;                  //< Number of elements.
  };
  std::vector<Segment> m_segments;

  struct Wire {
    double x, y; //< Coordinates of the centre.
    double r;    //< Radius.
    double v;    //< Potential.
    double q;    //< Charge.
  };
  std::vector<Wire> m_wires;

  struct Element {
    double x, y; //< Coordinates of the element centre (collocation point).
    double a;    //< Half-length.
    double cphi; //< Rotation.
    double sphi; //< Rotation.
    double q;    //< Charge density (solution).
    std::pair<BC, double> bc; //< Boundary condition.
    double lambda;            //< Ratio of dielectric permittivities.
  };
  std::vector<Element> m_elements;

  /// Split/merge overlapping segments.
  void EliminateOverlaps(std::vector<Segment>& segments);
  /// Create elements from a straight-line segment.
  bool Discretise(const Segment& segment, std::vector<Element>& elements,
                  const double lambda, const unsigned int ndiv);

  bool ComputeInfluenceMatrix(std::vector<std::vector<double> >& infmat) const;
  void SplitElement(Element& oldElement, std::vector<Element>& elements);
  bool InvertMatrix(std::vector<std::vector<double> >& influenceMatrix,
                    std::vector<std::vector<double> >& inverseMatrix) const;
  bool LUDecomposition(std::vector<std::vector<double> >& mat,
                       std::vector<int>& index) const;
  void LUSubstitution(const std::vector<std::vector<double> >& mat,
                      const std::vector<int>& index,
                      std::vector<double>& col) const;

  bool Solve(const std::vector<std::vector<double> >& inverseMatrix,
             const std::vector<double>& bc);
  bool CheckConvergence(const double tol, std::vector<bool>& ok);

  double LinePotential(const double a, const double x, const double y) const;
  double WirePotential(const double r0, const double x, const double y) const;
  void LineFlux(const double a, const double x, const double y, double& ex,
                double& ey) const;
  void WireFlux(const double r0, const double x, const double y, double& ex,
                double& ey) const;

  void Reset() override;
  void UpdatePeriodicity() override;
  void ToLocal(const double xIn, const double yIn,
               const double cphi, const double sphi,
               double& xOut, double& yOut) const;
  void ToGlobal(const double xIn, const double yIn, 
                const double cphi, const double sphi,
                double& xOut, double& yOut) const;

};
}

#endif
