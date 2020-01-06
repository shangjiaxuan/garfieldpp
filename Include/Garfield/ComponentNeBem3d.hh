#ifndef G_COMPONENT_NEBEM_3D_H
#define G_COMPONENT_NEBEM_3D_H
#include "ComponentBase.hh"

namespace Garfield {

class ComponentNeBem3d : public ComponentBase {
 public:
  /// Constructor
  ComponentNeBem3d();
  /// Destructor
  ~ComponentNeBem3d() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  bool GetVoltageRange(double& vmin, double& vmax) override;

  unsigned int GetNumberOfPrimitives() const { return m_primitives.size(); }
  bool GetPrimitive(const unsigned int i, double& a, double& b, double& c,
                    std::vector<double>& xv, std::vector<double>& yv,
                    std::vector<double>& zv, int& interface, double& v,
                    double& q, double& lambda) const;
  unsigned int GetNumberOfElements() const { return m_elements.size(); }
  bool GetElement(const unsigned int i, 
                  std::vector<double>& xv, std::vector<double>& yv,
                  std::vector<double>& zv, int& interface, double& bc,
                  double& lambda) const;

  void Reset() override;
  void UpdatePeriodicity() override;

  /// Retrieve surface panels, remove contacts and cut polygons to rectangles
  /// and right-angle triangles.
  bool Initialise();

  /// Set the default value of the target linear size of the elements
  /// produced by neBEM's discretisation process.
  void SetTargetElementSize(const double length);

 private:
  struct Primitive {
    /// Perpendicular vector
    double a, b, c;
    /// X-coordinates of vertices
    std::vector<double> xv;
    /// Y-coordinates of vertices
    std::vector<double> yv;
    /// Z-coordinates of vertices
    std::vector<double> zv;
    /// Interface type.
    int interface;
    /// Potential
    double v;
    /// Charge
    double q;
    /// Ratio of dielectric constants
    double lambda;
    /// Target element size.
    double elementSize;
  };
  /// List of primitives.
  std::vector<Primitive> m_primitives;

  struct Element {
    /// Local origin.
    std::array<double, 3> origin;
    double lx;
    double lz;
    /// Area.
    double dA;
    /// Direction cosines.
    std::array<std::array<double, 3>, 3> dcos;
    /// X-coordinates of vertices
    std::vector<double> xv;
    /// Y-coordinates of vertices
    std::vector<double> yv;
    /// Z-coordinates of vertices
    std::vector<double> zv;
    /// Interface type.
    int interface;
    /// Ratio of dielectric permittivities.
    double lambda;
    /// Collocation point.
    std::array<double, 3> collocationPoint;
    /// Boundary condition.
    double bc;
    /// Fixed charge density.
    double assigned;
    /// Solution (accumulated charge).
    double solution;
  };
  /// List of elements.
  std::vector<Element> m_elements;

  static constexpr double MinDist = 1.e-6;
  /// Target size of elements [cm].
  double m_targetElementSize = 50.0e-4;
  /// Smallest number of elements produced along the axis of a primitive. 
  unsigned int m_minNbElementsOnLength = 1;
  /// Largest number of elements produced along the axis of a primitive. 
  unsigned int m_maxNbElementsOnLength = 100; 

  /// Isolate the parts of polygon 1 that are not hidden by 2 and vice versa.
  bool EliminateOverlaps(const Panel& panel1, const Panel& panel2,
                         std::vector<Panel>& panelsOut,
                         std::vector<int>& itypo);

  bool TraceEnclosed(const std::vector<double>& xl1,
                     const std::vector<double>& yl1,
                     const std::vector<double>& xl2,
                     const std::vector<double>& yl2, const Panel& originalPanel,
                     std::vector<Panel>& newPanels) const;

  void TraceNonOverlap(
      const std::vector<double>& xp1, const std::vector<double>& yp1,
      const std::vector<double>& xl1, const std::vector<double>& yl1,
      const std::vector<double>& xl2, const std::vector<double>& yl2,
      const std::vector<int>& flags1, const std::vector<int>& flags2,
      const std::vector<int>& links1, const std::vector<int>& links2,
      std::vector<bool>& mark1, int ip1, const Panel& originalPanel,
      std::vector<Panel>& newPanels) const;

  void TraceOverlap(
      const std::vector<double>& xp1, const std::vector<double>& yp1,
      const std::vector<double>& xp2, const std::vector<double>& yp2,
      const std::vector<double>& xl1, const std::vector<double>& yl1,
      const std::vector<double>& xl2, const std::vector<double>& yl2,
      const std::vector<int>& flags1, const std::vector<int>& links1,
      const std::vector<int>& links2, std::vector<bool>& mark1, int ip1,
      int ip2, const Panel& originalPanel, std::vector<Panel>& newPanels) const;

  /// Split a polygon into rectangles and right-angled triangles.
  bool MakePrimitives(const Panel& panelIn,
                      std::vector<Panel>& panelsOut) const;

  /// Check whether a polygon contains parallel lines.
  /// If it does, split it in rectangular and non-rectangular parts.
  bool SplitTrapezium(const Panel panelIn, std::vector<Panel>& stack,
                      std::vector<Panel>& panelsOut, const double epsang) const;

  unsigned int NbOfSegments(const double length, const double target) const;
  bool DiscretizeWire(const Primitive& primitive, const double targetSize,
                      std::vector<Element>& elements) const;
  bool DiscretizeTriangle(const Primitive& primitive, const double targetSize,
                          std::vector<Element>& elements) const;
  bool DiscretizeRectangle(const Primitive& prim, const double targetSize,
                           std::vector<Element>& elements) const;
  int InterfaceType(const Solid::BoundaryCondition bc) const;
};
}

#endif
