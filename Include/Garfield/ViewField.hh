#ifndef G_VIEW_FIELD
#define G_VIEW_FIELD

#include "ViewBase.hh"

namespace Garfield {

class Sensor;
class ComponentBase;

/// Visualize the potential or electric field of a component or sensor.

class ViewField : public ViewBase {
 public:
  /// Constructor.
  ViewField();
  /// Destructor.
  ~ViewField() = default;

  /// Set the sensor from which to retrieve the field.
  void SetSensor(Sensor* s);
  /// Set the component from which to retrieve the field.
  void SetComponent(ComponentBase* c);

  /// Set the plot limits for the potential.
  void SetVoltageRange(const double vmin, const double vmax);
  /// Set the plot limits for the electric field.
  void SetElectricFieldRange(const double emin, const double emax);
  /// Set the plot limits for the weighting field.
  void SetWeightingFieldRange(const double wmin, const double wmax);

  /// Set the viewing area (in local coordinates of the current viewing plane).
  void SetArea(const double xmin, const double ymin, const double xmax,
               const double ymax);
  /// Set the viewing area based on the bounding box of the sensor/component.
  void SetArea() { m_hasUserArea = false; }

  /** Set the projection (viewing plane).
    * \param fx,fy,fz normal vector
    * \param x0,y0,z0 in-plane point
    */
  void SetPlane(const double fx, const double fy, const double fz,
                const double x0, const double y0, const double z0);
  /// Set the default viewing plane (\f$x\f$-\f$y\f$ at \f$z = 0\f$).
  void SetDefaultProjection();
  /// Rotate the viewing plane (angle in radian).
  void Rotate(const double angle);

  /// Set the number of contour levels.
  void SetNumberOfContours(const unsigned int n);
  /// Set the number of points used for drawing 1D functions.
  void SetNumberOfSamples1d(const unsigned int n);
  /// Set the number of points used for drawing 2D functions.
  void SetNumberOfSamples2d(const unsigned int nx, const unsigned int ny);

  /** Make a contour plot of the electric potential or field.
    * \param option quantity to be plotted
    * - potential: "v", "voltage", "p", "potential"
    * - magnitude of the electric field: "e", "field"
    * - x-component of the electric field: "ex"
    * - y-component of the electric field: "ey"
    * - z-component of the electric field: "ez"
    **/
  void PlotContour(const std::string& option = "v");
  /** Make a 2D plot of the electric potential or field.
    * \param option quantity to be plotted (see PlotContour)
    * \param drawopt option string passed to TF2::Draw
    **/
  void Plot(const std::string& option = "v",
            const std::string& drawopt = "arr");
  /** Make a 1D plot of the electric potential or field along a line.
    * \param x0,y0,z0 starting point
    * \param x1,y1,z1 end point
    * \param option quantity to be plotted (see PlotContour)
    **/
  void PlotProfile(const double x0, const double y0, const double z0,
                   const double x1, const double y1, const double z1,
                   const std::string& option = "v");

  /** Make a contour plot of the weighting potential or field.
    * \param label identifier of the electrode
    * \param option quantity to be plotted (see PlotContour)
    **/
  void PlotContourWeightingField(const std::string& label,
                                 const std::string& option);
  /** Make a 2D plot of the weighting potential or field.
    * \param label identifier of the electrode
    * \param option quantity to be plotted (see PlotContour)
    * \param drawopt option string passed to TF2::Draw
    **/
  void PlotWeightingField(const std::string& label, const std::string& option,
                          const std::string& drawopt);

  /** Make a 1D plot of the weighting potential or field along a line.
    * \param label identifier of the electrode
    * \param x0,y0,z0 starting point
    * \param x1,y1,z1 end point
    * \param option quantity to be plotted (see PlotContour)
    **/
  void PlotProfileWeightingField(const std::string& label, const double x0,
                                 const double y0, const double z0,
                                 const double x1, const double y1,
                                 const double z1,
                                 const std::string& option = "v");
  /// Determine the range of the potential/field automatically (true)
  /// or set it explicitly (false). 
  void EnableAutoRange(const bool on = true) { m_useAutoRange = on; }

  /** Make use of the status flag returned by the sensor/component.
    * \param v0 Value to be used for regions with status != 0.
    */
  void EnableAcknowledgeStatus(const double v0 = 0.) {
    m_useStatus = true;
    m_vBkg = v0;
  }
  /// Ignore the status flag returned by the sensor/component.
  void DisableAcknowledgeStatus() { m_useStatus = false; }

 private:
  enum class PlotType { Potential = 0, Magnitude, Ex, Ey, Ez, Unknown };

  bool m_useAutoRange = true;
  bool m_useStatus = false;
  double m_vBkg = 0.;

  // Sensor
  Sensor* m_sensor = nullptr;
  ComponentBase* m_component = nullptr;

  // Projection for viewing
  double m_project[3][3];
  double m_plane[4];
  char m_xLabel[50], m_yLabel[50], m_description[50];

  // Plot area
  bool m_hasUserArea = false;
  double m_xmin = -1., m_ymin = -1.;
  double m_xmax = 1., m_ymax = 1.;

  // Function range
  double m_vmin = 0., m_vmax = 100.;
  double m_emin = 0., m_emax = 10000.;
  double m_wmin = 0., m_wmax = 100.;

  // Number of contours
  unsigned int m_nContours = 20;
  // Number of points used to draw the functions
  unsigned int m_nSamples1d = 1000;
  unsigned int m_nSamples2dX = 200;
  unsigned int m_nSamples2dY = 200;

  void Labels();
  bool Range();
  void Draw2d(const std::string& option, const bool contour,
              const bool wfield, const std::string& electrode,
              const std::string& drawopt);
  void DrawProfile(const double x0, const double y0, const double z0,
                   const double x1, const double y1, const double z1,
                   const std::string& option, const bool wfield,
                   const std::string& electrode);
  PlotType GetPlotType(const std::string& option, std::string& title) const;
  double Field(const double x, const double y, const double z,
               const PlotType plotType) const;
  double Wfield(const double x, const double y, const double z,
                const PlotType plotType, const std::string& electrode) const;

};
}
#endif
