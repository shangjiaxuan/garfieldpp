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

  /// Set the x- and y-axis limits
  /// (in local coordinates of the current viewing plane, if applicable).
  void SetArea(const double xmin, const double ymin, const double xmax,
               const double ymax);
  /// Set a bounding box (if applicable).
  virtual void SetArea(const double xmin, const double ymin, const double zmin,
                       const double xmax, const double ymax, const double zmax);
  /// Use default x- and y-axis limits
  /// (based on the bounding box of the sensor/component, if applicable).
  void SetArea() {
    m_userBox = false; 
    m_userPlotLimits = false; 
  }

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

  // Plot axis limits. 
  bool m_userPlotLimits = false;
  double m_xMinPlot = -1., m_xMaxPlot = 1.;
  double m_yMinPlot = -1., m_yMaxPlot = 1.;

  // Bounding box.
  bool m_userBox = false;
  double m_xMinBox = -1., m_xMaxBox = 1.;
  double m_yMinBox = -1., m_yMaxBox = 1.; 
  double m_zMinBox = -1., m_zMaxBox = 1.;

  // Viewing plane (FPROJ).
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
