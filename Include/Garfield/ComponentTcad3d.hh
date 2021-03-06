#ifndef G_COMPONENT_TCAD_3D_H
#define G_COMPONENT_TCAD_3D_H

#include <array>
#include <memory>

#include "Component.hh"
#include "TetrahedralTree.hh"

namespace Garfield {

/// Interpolation in a three-dimensional field map created by Sentaurus Device.

class ComponentTcad3d : public Component {
 public:
  /// Constructor
  ComponentTcad3d();
  /// Destructor
  ~ComponentTcad3d() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  Medium* GetMedium(const double x, const double y, const double z) override;

  bool GetVoltageRange(double& vmin, double& vmax) override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override;

  /** Import mesh and field map from files.
    * \param gridfilename name of the .grd file containing the mesh 
    * \param datafilename name of the .dat file containing the nodal solution
    */
  bool Initialise(const std::string& gridfilename,
                  const std::string& datafilename);

  /** Import field maps defining the weighting field and potential.
    * \param datfile1 .dat file containing the field map at nominal bias.
    * \param datfile2 .dat file containing the field map for a configuration 
                      with the potential at the electrode to be read out
                      increased by a small voltage dv.
    * \param dv increase in electrode potential between the two field maps. 
    *
    * The field maps must use the same mesh as the drift field. 
    */ 
  bool SetWeightingField(const std::string& datfile1,
                         const std::string& datfile2, const double dv);

  /// List all currently defined regions.
  void PrintRegions();
  /// Get the number of regions in the device.
  size_t GetNumberOfRegions() const { return m_regions.size(); }
  void GetRegion(const size_t ireg, std::string& name,
                 bool& active) const;
  void SetDriftRegion(const size_t ireg);
  void UnsetDriftRegion(const size_t ireg);
  /// Set the medium for a given region
  void SetMedium(const size_t ireg, Medium* m);
  /// Get the medium for a given region
  bool GetMedium(const size_t ireg, Medium*& m) const;

  size_t GetNumberOfElements() const { return m_elements.size(); }
  bool GetElement(const size_t i, double& vol, double& dmin, double& dmax,
                  int& type) const;
  bool GetElement(const size_t i, double& vol, double& dmin, double& dmax,
                  int& type, std::vector<size_t>& nodes, int& reg) const;
  size_t GetNumberOfNodes() const { return m_vertices.size(); }
  bool GetNode(const size_t i, double& x, double& y, double& z, double& v,
               double& ex, double& ey, double& ez) const;

 private:
  // Max. number of vertices per element
  static constexpr size_t nMaxVertices = 7;

  // Regions
  struct Region {
    // Name of region (from Tcad)
    std::string name;
    // Flag indicating if the region is active (i. e. a drift medium)
    bool drift;
    Medium* medium;
  };
  std::vector<Region> m_regions;

  // Vertex coordinates [cm].
  std::vector<std::array<double, 3> > m_vertices;
  // Potential [V] at each vertex.
  std::vector<double> m_potential;
  // Electric field [V / cm].
  std::vector<std::array<double, 3> > m_efield;
  // Weighting field and potential at each vertex.
  std::vector<std::array<double, 3> > m_wf;
  std::vector<double> m_wp;

  // Elements
  struct Element {
    // Indices of vertices
    size_t vertex[nMaxVertices];
    // Type of element
    // 1: Segment (line)
    // 2: Triangle
    // 3: Rectangle
    // 4: Polygon
    // 5: Tetrahedron
    // 6: Pyramid
    // 7: Prism
    // 8: Brick
    // 9: Tetrabrick
    // 10: Polyhedron
    // Only types 2 and 5 are supported by this class.
    int type;
    // Associated region
    unsigned int region;
    // Bounding box
    float xmin, xmax;
    float ymin, ymax;
    float zmin, zmax;
  };
  std::vector<Element> m_elements;

  // Face
  struct Face {
    // Indices of edges
    int edge[4];
    int type;
  };

  // Voltage range
  double m_pMin = 0.;
  double m_pMax = 0.;

  // Bounding box
  std::array<double, 3> m_bbMin = {{0., 0., 0.}};
  std::array<double, 3> m_bbMax = {{0., 0., 0.}};

  // Tetrahedral tree.
  std::unique_ptr<TetrahedralTree> m_tree;

  void Reset() override;
  void UpdatePeriodicity() override;

  size_t FindElement(const double x, const double y, const double z,
                     std::array<double, nMaxVertices>& w) const;
  bool InElement(const double x, const double y, const double z,
                 const Element& element, 
                 std::array<double, nMaxVertices>& w) const {
    if (x < element.xmin || x > element.xmax || 
        y < element.ymin || y > element.ymax || 
        z < element.zmin || z > element.zmax) {
      return false;
    }
    bool inside = false;
    switch (element.type) {
      case 2:
        if (InTriangle(x, y, z, element, w)) inside = true;
        break;
      case 5:
        if (InTetrahedron(x, y, z, element, w)) inside = true;
        break;
      default:
        std::cerr << m_className << "::InElement:\n"
                  << "    Invalid element type (" << element.type << ").\n";
        break;
    }
    return inside;
  }
  bool InTetrahedron(const double x, const double y, const double z,
                     const Element& element, 
                     std::array<double, nMaxVertices>& w) const;
  bool InTriangle(const double x, const double y, const double z,
                  const Element& element, 
                  std::array<double, nMaxVertices>& w) const;
  
  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<double>& field, double& f);
  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<std::array<double, 3> >& field, 
                   double& fx, double& fy, double& fz);

  bool LoadGrid(const std::string& gridfilename);
  bool LoadData(const std::string& datafilename);
  bool ReadDataset(std::ifstream& datafile, const std::string& dataset);
  bool LoadWeightingField(const std::string& datafilename,
                          std::vector<std::array<double, 3> >& wf, 
                          std::vector<double>& wp);
  void Cleanup();

  void MapCoordinates(double& x, double& y, double& z, bool& xmirr, bool& ymirr,
                      bool& zmirr) const;
  bool InBoundingBox(const double x, const double y, const double z) const {
    bool inside = true;
    if (x < m_bbMin[0] || x > m_bbMax[0] || 
        y < m_bbMin[1] || y > m_bbMax[1] ||
        z < m_bbMin[2] || z > m_bbMax[2]) {
      inside = false;
    }
    return inside;
  }
  size_t FindRegion(const std::string& name) const;
};
}
#endif
