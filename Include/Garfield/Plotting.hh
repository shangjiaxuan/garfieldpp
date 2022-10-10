#ifndef G_PLOTTING_H
#define G_PLOTTING_H

#include "PlottingEngine.hh"
#include "exports.h"

namespace Garfield {

GARFIELD_EXTERNAL_SYMBOL PlottingEngine plottingEngine;

inline void SetDefaultStyle() { 
  plottingEngine.SetDefaultStyle();
}

inline void SetSerif() {
  plottingEngine.SetSerif();
}

inline void SetSansSerif() {
  plottingEngine.SetSansSerif();
}

}

#endif
