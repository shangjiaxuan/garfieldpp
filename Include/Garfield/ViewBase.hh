#ifndef G_VIEW_BASE
#define G_VIEW_BASE

#include <TCanvas.h>

namespace Garfield {

/// Base class for visualization classes.

class ViewBase {
 public:
  /// Default constructor.
  ViewBase() = delete;
  /// Constructor.
  ViewBase(const std::string& name);
  /// Destructor.
  virtual ~ViewBase();

  /// Set the canvas to be painted on.
  void SetCanvas(TCanvas* c);
  /// Retrieve the canvas.
  TCanvas* GetCanvas();

  /// Switch on/off debugging output.
  void EnableDebugging(const bool on = true) { m_debug = on; }

 protected:
  std::string m_className = "ViewBase";

  // Options
  bool m_debug = false;

  // Viewing plane
  double m_proj[3][3];

  // Find an unused function name.
  std::string FindUnusedFunctionName(const std::string& s) const;
  // Find an unused histogram name.
  std::string FindUnusedHistogramName(const std::string& s) const;
 private:
  // Canvas
  TCanvas* m_canvas = nullptr;
  bool m_hasExternalCanvas = false;

};
}
#endif
