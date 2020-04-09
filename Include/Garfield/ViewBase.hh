#ifndef G_VIEW_BASE
#define G_VIEW_BASE

#include <array>
#include <string>

#include <TCanvas.h>

namespace Garfield {

class Sensor;
class ComponentBase;

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
  /// Set the viewing plane to x-y. 
  void SetPlaneXY();
  /// Set the viewing plane to x-z. 
  void SetPlaneXZ();
  /// Set the viewing plane to y-z. 
  void SetPlaneYZ();

  /// Switch on/off debugging output.
  void EnableDebugging(const bool on = true) { m_debug = on; }

  /// Find an unused function name.
  static std::string FindUnusedFunctionName(const std::string& s);
  /// Find an unused histogram name.
  static std::string FindUnusedHistogramName(const std::string& s);
  /// Find an unused canvas name.
  static std::string FindUnusedCanvasName(const std::string& s);

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
  std::array<std::array<double, 3>, 3> m_proj{{
    {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 0}}
  }};
  std::array<double, 4> m_plane{{0, 0, 1, 0}};
  // Matrix used for projections (FPRMAT).
  std::array<std::array<double, 3>, 3> m_prmat{{
    {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}
  }};

  // Update and invert the projection matrix.
  void UpdateProjectionMatrix();
  // Determine plane coordinates.
  void ToPlane(const double x, const double y, const double z,
               double& xp, double& yp) const;
  // X-axis label for the current viewing plane.
  std::string LabelX();
  // Y-axis label for the current viewing plane.
  std::string LabelY();
  // Description of the current viewing plane.
  std::string PlaneDescription();

  bool PlotLimits(Sensor* sensor, double& xmin, double& ymin,
                  double& xmax, double& ymax) const;
  bool PlotLimits(ComponentBase* cmp, double& xmin, double& ymin,
                  double& xmax, double& ymax) const;
  bool PlotLimitsFromUserBox(double& xmin, double& ymin, 
                             double& xmax, double& ymax) const;
  bool PlotLimits(std::array<double, 3>& bbmin,
                  std::array<double, 3>& bbmax,
                  double& xmin, double& ymin,
                  double& xmax, double& ymax) const;
 private:
  // Canvas
  TCanvas* m_canvas = nullptr;
  bool m_hasExternalCanvas = false;

};
}
#endif
