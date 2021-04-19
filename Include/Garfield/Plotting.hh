#ifndef G_PLOTTING_H
#define G_PLOTTING_H

#include "PlottingEngineRoot.hh"

namespace Garfield {

extern PlottingEngineRoot plottingEngine;

inline void SetDefaultStyle() { 
  plottingEngine.SetDefaultStyle();
}

}

#endif
