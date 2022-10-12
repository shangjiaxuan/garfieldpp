#ifndef G_COMPONENT_ANSYS123_H
#define G_COMPONENT_ANSYS123_H

#include "ComponentFieldMap.hh"

namespace Garfield {

/// Component for importing and interpolating three-dimensional ANSYS field maps.

class ComponentAnsys123 : public ComponentFieldMap {
 public:
  /// Constructor
  ComponentAnsys123();
  /// Destructor
  ~ComponentAnsys123() {}

  bool Initialise(const std::string& elist = "ELIST.lis",
                  const std::string& nlist = "NLIST.lis",
                  const std::string& mplist = "MPLIST.lis",
                  const std::string& prnsol = "PRNSOL.lis", 
                  const std::string& unit = "cm");

  bool SetWeightingField(const std::string& prnsol, const std::string& label);

};
}
#endif
