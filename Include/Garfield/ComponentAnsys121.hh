#ifndef G_COMPONENT_ANSYS121_H
#define G_COMPONENT_ANSYS121_H

#include "ComponentFieldMap.hh"

namespace Garfield {

/// Component for importing and interpolating two-dimensional ANSYS field maps.

class ComponentAnsys121 : public ComponentFieldMap {
 public:
  /// Constructor
  ComponentAnsys121();
  /// Destructor
  ~ComponentAnsys121() {}

  Medium* GetMedium(const double x, const double y, const double z) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;

  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  /** Import a field map.
    * \param elist name of the file containing the list of elements
    * \param nlist name of the file containing the list of nodes
    * \param mplist name of the file containing the list of materials
    * \param prnsol name of the file containing the nodal solutions
    * \param unit length unit
    */ 
  bool Initialise(std::string elist = "ELIST.lis",
                  std::string nlist = "NLIST.lis",
                  std::string mplist = "MPLIST.lis",
                  std::string prnsol = "PRNSOL.lis", std::string unit = "cm");
  /// Import a weighting field map.
  bool SetWeightingField(std::string prnsol, std::string label);
  /// Set the limits of the active region along z.
  void SetRangeZ(const double zmin, const double zmax);

};
}

#endif
