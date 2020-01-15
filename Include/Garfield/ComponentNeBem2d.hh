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

  /** Add a conducting line segment.
    * \param x0,y0,x1,y1 coordinates of start and end point.
    * \param v applied potential.
    */
  bool AddSegment(const double x0, const double y0, const double x1,
                  const double y1, const double v);
  /** Add a wire.
    * \param x,y centre of the wire.
    * \param d wire diameter.
    * \param v applied potential.
    */
  bool AddWire(const double x, const double y, const double d, const double v);
  /** Add an area bounded by a polygon.
    * \param xp,yp x/y-coordinates of the vertices.
    * \param medium pointer to the medium associated to the region. 
    * \param bctype 1: fixed voltage, 2: fixed charge density, 
                    3: floating conductor, 4: dielectric-dielectric interface.
    * \param bcval applied potential or charge density.
    */ 
  bool AddPolygon(const std::vector<double>& xp,
                  const std::vector<double>& yp, Medium* medium,
                  const unsigned int bctype = 4, const double bcval = 0.);

  bool Initialise();

  void SetNumberOfDivisions(const unsigned int ndiv);
  void SetNumberOfCollocationPoints(const unsigned int ncoll);
  void SetMinimumElementSize(const double min);
  void EnableAutoResizing(const bool on = true) { m_autoSize = on; }
  void EnableRandomCollocation(const bool on = true) {
    m_randomCollocation = on;
  }
  void SetMaxNumberOfIterations(const unsigned int niter);

  unsigned int GetNumberOfRegions() const { return m_regions.size(); }
  unsigned int GetNumberOfSegments() const { return m_segments.size(); }
  unsigned int GetNumberOfWires() const { return m_wires.size(); }
  unsigned int GetNumberOfElements() const { return m_elements.size(); }

 private:
  static const double InvEpsilon0;
  static const double InvTwoPiEpsilon0;

  unsigned int m_nDivisions = 5;
  unsigned int m_nCollocationPoints = 3;
  double m_minSize = 1.e-3;
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
  };
  std::vector<Region> m_regions;

  struct Segment {
    std::array<double, 2> x0;  //< Coordinates of the start point.
    std::array<double, 2> x1;  //< Coordinates of the end point.
    int region1;               //< Inner region. 
    int region2;               //< Outer region.
    std::pair<BC, double> bc;  //< Applied boundary condition.
  };
  std::vector<Segment> m_segments;

  struct Wire {
    double x, y; //< Coordinates of the centre.
    double r;    //< Radius.
    double v;    //< Potential.
  };
  std::vector<Wire> m_wires;

  struct Element {
    bool wire;   //< Flag whether this element is a wire or a straight line.
    double x, y; //< Coordinates of the element centre (collocation point).
    double a;    //< Half-length or radius.
    double cphi; //< Rotation.
    double sphi; //< Rotation.
    double q;    //< Charge density (solution).
    std::pair<BC, double> bc; //< Boundary condition.
    double lambda;            //< Ratio of dielectric permittivities.
  };
  std::vector<Element> m_elements;

  bool m_matrixInversionFlag = false;

  /// Split/merge overlapping segments.
  void EliminateOverlaps(std::vector<Segment>& segments);
  /// Create elements from a straight-line segment.
  bool Discretise(const Segment& segment, std::vector<Element>& elements,
                  const double lambda);

  bool ComputeInfluenceMatrix(std::vector<std::vector<double> >& infmat) const;
  void SplitElement(Element& element);
  bool InvertMatrix(std::vector<std::vector<double> >& influenceMatrix,
                    std::vector<std::vector<double> >& inverseMatrix) const;
  bool LUDecomposition(std::vector<std::vector<double> >& mat,
                       std::vector<int>& index) const;
  void LUSubstitution(const std::vector<std::vector<double> >& mat,
                      const std::vector<int>& index,
                      std::vector<double>& col) const;

  bool Solve(const std::vector<std::vector<double> >& inverseMatrix,
             const std::vector<double>& bc);
  bool CheckConvergence() const;

  // Compute the field
  bool Flux(const int gt, const double len, const double cphi,
            const double sphi, const double x, const double y, 
            double& ex, double& ey) const {
    double fx = 0., fy = 0.;
    switch (gt) {
      case 0:
        LineFlux(len, x, y, fx, fy);
        break;
      case 1:
        WireFlux(len, x, y, fx, fy);
        break;
      default:
        return false;
    }
    // Transformation to global coordinate system
    ToGlobal(fx, fy, cphi, sphi, ex, ey);
    return true;
  }

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
