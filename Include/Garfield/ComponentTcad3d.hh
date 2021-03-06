#ifndef G_COMPONENT_TCAD_3D_H
#define G_COMPONENT_TCAD_3D_H

#include <memory>

#include "ComponentTcadBase.hh"
#include "TetrahedralTree.hh"

namespace Garfield {

/// Interpolation in a three-dimensional field map created by Sentaurus Device.

class ComponentTcad3d : public ComponentTcadBase<3> {
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

  /** Retrieve the properties of an element.
    * \param i index of the element
    * \param vol volume
    * \param dmin smallest length in the element
    * \param dmax largest length in the element
    * \param type element type
    * \param nodes indices of the constituent vertices
    * \param reg region
    */
  bool GetElement(const size_t i, double& vol, double& dmin, double& dmax,
                  int& type, std::vector<size_t>& nodes, int& reg) const;
  /// Get the coordinates of a mesh node and the potential and 
  /// electric field at this node.
  bool GetNode(const size_t i, double& x, double& y, double& z, double& v,
               double& ex, double& ey, double& ez) const;

 private:
  // Face
  struct Face {
    // Indices of edges
    int edge[4];
    int type;
  };

  // Tetrahedral tree.
  std::unique_ptr<TetrahedralTree> m_tree;

  void Reset() override;

  size_t FindElement(const double x, const double y, const double z,
                     std::array<double, nMaxVertices>& w) const;
  bool InElement(const double x, const double y, const double z,
                 const Element& element, 
                 std::array<double, nMaxVertices>& w) const {
    if (x < element.bbMin[0] || x > element.bbMax[0] || 
        y < element.bbMin[1] || y > element.bbMax[1] || 
        z < element.bbMin[2] || z > element.bbMax[2]) {
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
};
}
#endif
