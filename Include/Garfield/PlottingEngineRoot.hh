#ifndef G_PLOTTING_ENGINE_ROOT_H
#define G_PLOTTING_ENGINE_ROOT_H

#include <TStyle.h>

#include "PlottingEngine.hh"

namespace Garfield {

/// Definition of styles.

class PlottingEngineRoot : public PlottingEngine {
 public:
  /// Constructor
  PlottingEngineRoot();
  /// Destructor
  virtual ~PlottingEngineRoot();

  /// Apply the default Garfield ROOT style.
  void SetDefaultStyle();

 private:
  TStyle m_garfieldStyle;
};
}

#endif
