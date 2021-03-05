#ifndef G_COMPONENT_TCAD_2D_H
#define G_COMPONENT_TCAD_2D_H

#include <array>
#include <memory>

#include "Component.hh"
#include "QuadTree.hh"

namespace Garfield {

/// Interpolation in a two-dimensional field map created by Sentaurus Device.

class ComponentTcad2d : public Component {
 public:
  /// Constructor
  ComponentTcad2d();
  /// Destructor
  ~ComponentTcad2d() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override {
    double v = 0.;
    ElectricField(x, y, z, ex, ey, ez, v, m, status);
  }

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  Medium* GetMedium(const double x, const double y, const double z) override;

  bool GetVoltageRange(double& vmin, double& vmax) override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override;
  void SetRangeZ(const double zmin, const double zmax);

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
  void PrintRegions() const;
  /// Get the number of regions in the device.
  size_t GetNumberOfRegions() const { return m_regions.size(); }
  void GetRegion(const size_t i, std::string& name, bool& active) const;
  void SetDriftRegion(const size_t ireg);
  void UnsetDriftRegion(const size_t ireg);
  /// Set the medium for a given region.
  void SetMedium(const size_t ireg, Medium* m);
  /// Get the medium for a given region.
  Medium* GetMedium(const size_t ireg) const;

  // Retrieve information about the mesh.
  size_t GetNumberOfElements() const { return m_elements.size(); }
  bool GetElement(const size_t i, double& vol, double& dmin, double& dmax,
                  int& type) const;
  bool GetElement(const size_t i, double& vol, double& dmin, double& dmax,
                  int& type, int& node1, int& node2, int& node3, int& node4,
                  int& reg) const;
  size_t GetNumberOfNodes() const { return m_vertices.size(); }
  bool GetNode(const size_t i, double& x, double& y, double& v,
               double& ex, double& ey) const;

  // Mobilities
  bool GetElectronMobility(const double x, const double y, const double z, 
                           double& mob);
  bool GetHoleMobility(const double x, const double y, const double z, 
                       double& mob);

  /// Switch use of the imported velocity map on/off.
  void EnableVelocityMap(const bool on) { m_useVelocityMap = on; }
  bool HasVelocityMap() const override { return m_useVelocityMap; }
  bool ElectronVelocity(const double x, const double y, const double z,
                        double& vx, double& vy, double& vz) override;
  bool HoleVelocity(const double x, const double y, const double z, 
                    double& vx, double& vy, double& vz) override;

  // Lifetime field maps
  bool GetElectronLifetime(const double x, const double y, const double z,
                           double& etau) override;
  bool GetHoleLifetime(const double x, const double y, const double z,
                       double& htau) override;

  // Trapping
  int GetNumberOfDonors() { return m_donors.size(); }
  int GetNumberOfAcceptors() { return m_acceptors.size(); }

  bool SetDonor(const size_t donorNumber, const double eXsec,
                const double hxSec, const double concentration);
  bool SetAcceptor(const size_t acceptorNumber, const double eXsec,
                   const double hxSec, const double concentration);

  /// Switch use of the imported trapping map on/off.
  void EnableAttachmentMap(const bool on) { m_useAttachmentMap = on; }
  bool HasAttachmentMap() const override;
  bool ElectronAttachment(const double x, const double y, const double z,
                          double& eta) override;
  bool HoleAttachment(const double x, const double y, const double z,
                      double& eta) override;

 private:
  // Max. number of vertices per element
  static constexpr size_t nMaxVertices = 4;

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
  std::vector<std::array<double, 2> > m_vertices;

  // Potential [V] at each vertex.
  std::vector<double> m_potential;
  // Electric field [V / cm].
  std::vector<std::array<double, 2> > m_efield;
  // Weighting field and potential at each vertex.
  std::vector<std::array<double, 2> > m_wf;
  std::vector<double> m_wp;
  // Velocities [cm / ns]
  std::vector<std::array<double, 2> > m_eVelocity; 
  std::vector<std::array<double, 2> > m_hVelocity;
  // Mobilities [cm2 / (V ns)]
  std::vector<double> m_eMobility;
  std::vector<double> m_hMobility; 
  // Trap occupations [dimensionless]
  std::vector<std::vector<float> > m_donorOcc;
  std::vector<std::vector<float> > m_acceptorOcc;
  // Lifetimes [ns]
  std::vector<double> m_eLifetime;
  std::vector<double> m_hLifetime;
  // Attachment coefficients [1 / cm]
  std::vector<double> m_eAttachment;
  std::vector<double> m_hAttachment;

  // Elements
  struct Element {
    // Indices of vertices
    int vertex[nMaxVertices];
    // Type of element
    // 0: Point
    // 1: Segment (line)
    // 2: Triangle
    // 3: Rectangle
    // 4: Polygon
    // Types 1 - 3 are supported by this class.
    int type;
    // Associated region
    unsigned int region;
    // Bounding box
    float xmin, xmax;
    float ymin, ymax;
  };
  std::vector<Element> m_elements;

  struct Defect {
    // Electron cross-section
    double xsece;
    // Hole cross-section
    double xsech;
    // Concentration
    double conc;
  };
  std::vector<Defect> m_donors;
  std::vector<Defect> m_acceptors;

  // Use velocity and trapping maps or not.
  bool m_useVelocityMap = false;
  bool m_useAttachmentMap = false;

  // Voltage range
  double m_pMin = 0.;
  double m_pMax = 0.;

  // Bounding box
  bool m_hasRangeZ = false;
  std::array<double, 3> m_bbMin = {{0., 0., 0.}};
  std::array<double, 3> m_bbMax = {{0., 0., 0.}};

  // Tetrahedral tree.
  std::unique_ptr<QuadTree> m_tree;

  // Element from the previous call
  int m_lastElement = 0;

  void Reset() override;
  void UpdatePeriodicity() override;
  size_t FindElement(const double x, const double y,
                     std::array<double, nMaxVertices>& w) const;
  // Check whether a point is inside a given element and calculate the
  // shape functions if it is.
  bool InElement(const double x, const double y, const Element& element,
                 std::array<double, nMaxVertices>& w) const;
  bool InRectangle(const double x, const double y, const Element& element,
                   std::array<double, nMaxVertices>& w) const;
  bool InTriangle(const double x, const double y, const Element& element,
                  std::array<double, nMaxVertices>& w) const;
  bool OnLine(const double x, const double y, const Element& element,
              std::array<double, nMaxVertices>& w) const;
  bool AtPoint(const double x, const double y, const Element& element,
               std::array<double, nMaxVertices>& w) const;

  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<double>& field, double& f);
  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<std::array<double, 2> >& field, 
                   double& fx, double& fy);

  bool LoadGrid(const std::string& gridfilename);
  bool LoadData(const std::string& datafilename);
  bool ReadDataset(std::ifstream& datafile, const std::string& dataset);
  bool LoadWeightingField(const std::string& datafilename,
                          std::vector<std::array<double, 2> >& wf,
                          std::vector<double>& wp);
  void Cleanup();

  size_t FindRegion(const std::string& name) const;

  void MapCoordinates(double& x, double& y, bool& xmirr, bool& ymirr) const;
  bool InsideBoundingBox(const double x, const double y) const { 
    bool inside = true;
    if (x < m_bbMin[0] || x > m_bbMax[0] || 
        y < m_bbMin[1] || y > m_bbMax[1]) {
      inside = false;
    }
    return inside;
  }
  void UpdateAttachment();
};
}
#endif
