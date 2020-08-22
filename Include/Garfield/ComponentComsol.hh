#pragma once

#include "ComponentFieldMap.hh"

namespace Garfield {

/// Component for importing and interpolating Comsol field maps.

class ComponentComsol : public ComponentFieldMap {
 public:
  /// Default constructor.
  ComponentComsol();
  /// Constructor from file names.
  ComponentComsol(const std::string& mesh, 
                  const std::string& mplist, const std::string& field,
                  const std::string& unit = "m");
  /// Destructor.
  ~ComponentComsol() {}

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

  Medium* GetMedium(const double x, const double y, const double z) override;

  /** Import a field map.
    * \param header name of the file containing the list of nodes
    * \param mplist name of the file containing the material properties
    * \param field name of the file containing the potentials at the nodes
    * \param unit length unit
    */ 
  bool Initialise(const std::string& header = "mesh.mphtxt",
                  const std::string& mplist = "dielectrics.dat",
                  const std::string& field = "field.txt",
                  const std::string& unit = "m");
  /// Import weighting field maps.
  bool SetWeightingField(const std::string& file, const std::string& label);

 protected:
  void UpdatePeriodicity() override { UpdatePeriodicityCommon(); }

  double GetElementVolume(const unsigned int i) override;
  void GetAspectRatio(const unsigned int i, double& dmin,
                      double& dmax) override;
 private:
  double m_unit = 100.;

};
}
