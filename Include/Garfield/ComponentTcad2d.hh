#ifndef G_COMPONENT_TCAD_2D_H
#define G_COMPONENT_TCAD_2D_H

#include <memory>

#include "ComponentTcadBase.hh"
#include "QuadTree.hh"

namespace Garfield {

/// Interpolation in a two-dimensional field map created by Sentaurus Device.

class ComponentTcad2d : public ComponentTcadBase<2> {
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

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override;
  /// Set the z-extent of the bounding box (default: unlimited).
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
  /// Get the coordinates of a mesh node and the potential 
  /// and electric field at this node.
  bool GetNode(const size_t i, double& x, double& y, double& v,
               double& ex, double& ey) const;

 private:
  // Bounding box
  bool m_hasRangeZ = false;

  // Tetrahedral tree.
  std::unique_ptr<QuadTree> m_tree;

  // Element from the previous call
  size_t m_lastElement = 0;

  void Reset() override;

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
                   const std::vector<double>& field, double& f) override;
  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<std::array<double, 2> >& field, 
                   double& fx, double& fy, double& fz) override;

  bool LoadGrid(const std::string& gridfilename);
  bool LoadData(const std::string& datafilename);
  bool ReadDataset(std::ifstream& datafile, const std::string& dataset);
  bool LoadWeightingField(const std::string& datafilename,
                          std::vector<std::array<double, 2> >& wf,
                          std::vector<double>& wp);
};
}
#endif
