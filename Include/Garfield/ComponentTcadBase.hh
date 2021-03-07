#ifndef G_COMPONENT_TCAD_BASE_H
#define G_COMPONENT_TCAD_BASE_H

#include <array>

#include "Component.hh"

namespace Garfield {

/// Interpolation in a field map created by Sentaurus Device.

template<size_t N> 
class ComponentTcadBase : public Component {
 public:
  /// Default constructor
  ComponentTcadBase() = delete;
  /// Constructor
  ComponentTcadBase(const std::string& name) : Component(name) {
    m_regions.reserve(10);
    m_vertices.reserve(10000);
    m_elements.reserve(10000);
  }
  /// Destructor
  virtual ~ComponentTcadBase() {}

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  bool GetVoltageRange(double& vmin, double& vmax) override;

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
  /// Get the name and "active volume" flag of a region.
  void GetRegion(const size_t ireg, std::string& name, bool& active) const;
  /// Make a region active ("driftable").
  void SetDriftRegion(const size_t ireg);
  /// Make a region inactive.
  void UnsetDriftRegion(const size_t ireg);
  /// Set the medium to be associated to a given region.
  void SetMedium(const size_t ireg, Medium* m);

  /// Get the number of elements in the mesh.
  size_t GetNumberOfElements() const { return m_elements.size(); }
  /// Get the number of vertices in the mesh.
  size_t GetNumberOfNodes() const { return m_vertices.size(); }
  
  /// Switch use of the imported velocity map on/off.
  void EnableVelocityMap(const bool on) { m_useVelocityMap = on; }
  bool HasVelocityMap() const override { return m_useVelocityMap; }
  bool ElectronVelocity(const double x, const double y, const double z,
                        double& vx, double& vy, double& vz) override;
  bool HoleVelocity(const double x, const double y, const double z, 
                    double& vx, double& vy, double& vz) override;

  /// Get the number of donor states found in the map.
  size_t GetNumberOfDonors() { return m_donors.size(); }
  /// Get the number of acceptor states found in the map.
  size_t GetNumberOfAcceptors() { return m_acceptors.size(); }

  bool SetDonor(const size_t donorNumber, const double eXsec,
                const double hxSec, const double concentration);
  bool SetAcceptor(const size_t acceptorNumber, const double eXsec,
                   const double hxSec, const double concentration);

  /// Switch use of the imported trapping map on/off.
  void EnableAttachmentMap(const bool on) { m_useAttachmentMap = on; }
  bool HasAttachmentMap() const override {
    return (m_useAttachmentMap && !(m_acceptors.empty() && m_donors.empty()));
  }
  bool ElectronAttachment(const double x, const double y, const double z,
                          double& eta) override;
  bool HoleAttachment(const double x, const double y, const double z,
                      double& eta) override;

  bool GetElectronLifetime(const double x, const double y, const double z,
                           double& etau) override;
  bool GetHoleLifetime(const double x, const double y, const double z,
                       double& htau) override;
  
  // Mobilities
  bool GetElectronMobility(const double x, const double y, const double z, 
                           double& mob);
  bool GetHoleMobility(const double x, const double y, const double z, 
                       double& mob);
 protected:
  // Max. number of vertices per element
  static constexpr size_t nMaxVertices = N == 2 ? 4 : 7;

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
  std::vector<std::array<double, N> > m_vertices;

  // Elements
  struct Element {
    // Indices of vertices
    unsigned int vertex[nMaxVertices];
    // Type of element
    // 0: Point
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
    // In 2D, types 1 - 3 are supported.
    // In 3D, only types 2 and 5 are supported.
    int type;
    // Associated region
    unsigned int region;
    // Bounding box
    std::array<float, N> bbMin;
    std::array<float, N> bbMax;
  };
  std::vector<Element> m_elements;

  // Potential [V] at each vertex.
  std::vector<double> m_potential;
  // Electric field [V / cm].
  std::vector<std::array<double, N> > m_efield;

  // Weighting field and potential at each vertex.
  std::vector<std::array<double, N> > m_wf;
  std::vector<double> m_wp;
  // Velocities [cm / ns]
  std::vector<std::array<double, N> > m_eVelocity; 
  std::vector<std::array<double, N> > m_hVelocity;
  // Mobilities [cm2 / (V ns)]
  std::vector<double> m_eMobility;
  std::vector<double> m_hMobility; 
  // Lifetimes [ns]
  std::vector<double> m_eLifetime;
  std::vector<double> m_hLifetime;
  // Trap occupations [dimensionless]
  std::vector<std::vector<float> > m_donorOcc;
  std::vector<std::vector<float> > m_acceptorOcc;
  // Attachment coefficients [1 / cm]
  std::vector<double> m_eAttachment;
  std::vector<double> m_hAttachment;
  
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
 
  // Use velocity map or not.
  bool m_useVelocityMap = false;
  // Use trapping map or not.
  bool m_useAttachmentMap = false;

  // Bounding box.
  std::array<double, 3> m_bbMin = {{0., 0., 0.}};
  std::array<double, 3> m_bbMax = {{0., 0., 0.}};
  
  // Voltage range
  double m_pMin = 0.;
  double m_pMax = 0.;

  void UpdatePeriodicity() override;

  void Cleanup();

  virtual bool Interpolate(const double x, const double y, const double z,
                           const std::vector<double>& field, double& f) = 0;
  virtual bool Interpolate(const double x, const double y, const double z,
                           const std::vector<std::array<double, N> >& field,
                           double& fx, double& fy, double& fz) = 0;

  size_t FindRegion(const std::string& name) const;
  void MapCoordinates(std::array<double, N>& x, 
                      std::array<bool, N>& mirr) const;
  bool InBoundingBox(const std::array<double, N>& x) const {
    for (size_t i = 0; i < N; ++i) {
      if (x[i] < m_bbMin[i] || x[i] > m_bbMax[i]) return false;
    }
    return true;
  }
  void UpdateAttachment();
  
  bool LoadWeightingField(const std::string& datafilename,
                          std::vector<std::array<double, N> >& wf,
                          std::vector<double>& wp);
};
}
#endif
