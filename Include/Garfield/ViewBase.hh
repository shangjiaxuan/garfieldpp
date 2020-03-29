#ifndef G_VIEW_BASE
#define G_VIEW_BASE

#include <array>
#include <string>

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

  /** Set the projection (viewing plane), if applicable.
    * \param fx,fy,fz normal vector
    * \param x0,y0,z0 in-plane point
    */
  virtual void SetPlane(const double fx, const double fy, const double fz,
                        const double x0, const double y0, const double z0);
  /// Rotate the viewing plane (angle in radian).
  void Rotate(const double angle);

  /// Switch on/off debugging output.
  void EnableDebugging(const bool on = true) { m_debug = on; }

 protected:
  std::string m_className = "ViewBase";

  // Options
  bool m_debug = false;

  // Viewing plane for field plots (FPROJ).
  // Default projection: x-y at z = 0.
  std::array<std::array<double, 3>, 3> m_proj = {{
    {1, 0, 0}, {0, 1, 0}, {0, 0, 0}
  }};
  std::array<double, 4> m_plane = {0, 0, 1, 0};

  // Find an unused function name.
  std::string FindUnusedFunctionName(const std::string& s) const;
  // Find an unused histogram name.
  std::string FindUnusedHistogramName(const std::string& s) const;

  // X-axis label for the current viewing plane.
  std::string LabelX();
  // Y-axis label for the current viewing plane.
  std::string LabelY();
  // Description of the current viewing plane.
  std::string PlaneDescription();

 private:
  // Canvas
  TCanvas* m_canvas = nullptr;
  bool m_hasExternalCanvas = false;

};
}
#endif
