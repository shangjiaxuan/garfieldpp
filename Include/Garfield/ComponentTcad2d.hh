#ifndef G_COMPONENT_TCAD_2D_H
#define G_COMPONENT_TCAD_2D_H

#include <array>

#include "Component.hh"

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
  unsigned int GetNumberOfRegions() const { return m_regions.size(); }
  void GetRegion(const unsigned int i, std::string& name, bool& active) const;
  void SetDriftRegion(const unsigned int ireg);
  void UnsetDriftRegion(const unsigned int ireg);
  /// Set the medium for a given region.
  void SetMedium(const unsigned int ireg, Medium* m);
  /// Get the medium for a given region.
  Medium* GetMedium(const unsigned int ireg) const;

  // Retrieve information about the mesh.
  unsigned int GetNumberOfElements() const { return m_elements.size(); }
  bool GetElement(const unsigned int i, double& vol, double& dmin, double& dmax,
                  int& type) const;
  bool GetElement(const unsigned int i, double& vol, double& dmin, double& dmax,
                  int& type, int& node1, int& node2, int& node3, int& node4,
                  int& reg) const;
  unsigned int GetNumberOfNodes() const { return m_vertices.size(); }
  bool GetNode(const unsigned int i, double& x, double& y, double& v,
               double& ex, double& ey) const;

  // Mobilities
  bool GetMobility(const double x, const double y, const double z, double& emob,
                   double& hmob);

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

  bool GetDonorOccupation(const double x, const double y, const double z,
                          const unsigned int donorNumber,
                          double& occupationFraction);
  bool GetAcceptorOccupation(const double x, const double y, const double z,
                             const unsigned int acceptorNumber,
                             double& occupationFraction);
  bool SetDonor(const unsigned int donorNumber, const double eXsec,
                const double hxSec, const double concentration);
  bool SetAcceptor(const unsigned int acceptorNumber, const double eXsec,
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
  static constexpr unsigned int nMaxVertices = 4;

  // Regions
  struct Region {
    // Name of region (from Tcad)
    std::string name;
    // Flag indicating if the region is active (i. e. a drift medium)
    bool drift;
    Medium* medium;
  };
  std::vector<Region> m_regions;

  // Vertices
  struct Vertex {
    // Coordinates [cm]
    double x, y;
    // Potential [V] and electric field [V / cm]
    double p, ex, ey;
    // Mobilities [cm2 / (V ns)]
    double emob, hmob;
    // Velocities [cm/ns]
    double eVx, eVy;
    double hVx, hVy;
    // Lifetimes [1/ns]
    double eTau, hTau;
    // Trap occupations [dimensionless]
    std::vector<float> donorOcc;
    std::vector<float> acceptorOcc;
  };
  std::vector<Vertex> m_vertices;

  // Weighting field and potential at each vertex.
  std::vector<std::array<double, 2> > m_wf;
  std::vector<double> m_wp;

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
    int region;
    std::vector<int> neighbours;
    // Bounding box
    double xmin, xmax;
    double ymin, ymax;
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

  // Available data.
  bool m_hasPotential = false;
  bool m_hasField = false;
  bool m_hasElectronMobility = false;
  bool m_hasHoleMobility = false;
  bool m_hasElectronVelocity = false;
  bool m_hasHoleVelocity = false;
  bool m_hasElectronLifetime = false;
  bool m_hasHoleLifetime = false;

  // Use velocity and trapping maps or not.
  bool m_useVelocityMap = false;
  bool m_useAttachmentMap = false;

  // Are all the cross-sections and concentrations valid and set.
  bool m_validTraps = false;

  // Voltage range
  double m_pMin = 0.;
  double m_pMax = 0.;

  // Bounding box
  bool m_hasRangeZ = false;
  double m_xMinBB, m_yMinBB, m_zMinBB;
  double m_xMaxBB, m_yMaxBB, m_zMaxBB;

  // Element from the previous call
  int m_lastElement = 0;

  void Reset() override;
  void UpdatePeriodicity() override;
  unsigned int FindElement(const double x, const double y,
                           std::array<double, nMaxVertices>& w) const;
  // Check whether a point is inside a given element and calculate the
  // shape functions if it is.
  bool CheckElement(const double x, const double y, const Element& element,
                    std::array<double, nMaxVertices>& w) const;
  bool CheckRectangle(const double x, const double y, const Element& element,
                      std::array<double, nMaxVertices>& w) const;
  bool CheckTriangle(const double x, const double y, const Element& element,
                     std::array<double, nMaxVertices>& w) const;
  bool CheckLine(const double x, const double y, const Element& element,
                 std::array<double, nMaxVertices>& w) const;
  bool CheckPoint(const double x, const double y, const Element& element,
                  std::array<double, nMaxVertices>& w) const;

  bool LoadGrid(const std::string& gridfilename);
  bool LoadData(const std::string& datafilename);
  bool ReadDataset(std::ifstream& datafile, const std::string& dataset);
  bool LoadWeightingField(const std::string& datafilename,
                          std::vector<std::array<double, 2> >& wf,
                          std::vector<double>& wp);
  void FindNeighbours();
  void Cleanup();

  int FindRegion(const std::string& name) const;

  void MapCoordinates(double& x, double& y, bool& xmirr, bool& ymirr) const;
  bool InsideBoundingBox(const double x, const double y, 
                         const double z) const {
    bool inside = true;
    if (x < m_xMinBB || x > m_xMaxBB || y < m_yMinBB || y > m_yMaxBB ||
        (m_hasRangeZ && (z < m_zMinBB || z > m_zMaxBB))) {
      inside = false;
    }
    return inside;
  }
  bool CheckTraps() const;
};
}
#endif
