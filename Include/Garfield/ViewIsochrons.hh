#ifndef G_VIEW_ISOCHRONS
#define G_VIEW_ISOCHRONS

#include "ViewBase.hh"

namespace Garfield {

class Sensor;
class ComponentBase;

/// Draw equal time contour lines.

class ViewIsochrons : public ViewBase {
 public:
  /// Constructor.
  ViewIsochrons();
  /// Destructor.
  ~ViewIsochrons() = default;

  /// Set the sensor.
  void SetSensor(Sensor* s);
  /// Set the component.
  void SetComponent(ComponentBase* c);

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

  /** Draw equal time contour lines.
    * \param tstep Time interval. 
    * \param points List of starting points.
    * \param rev If true, the drift time is measured from the end 
    *            points of the drift lines. 
    * \param colour Draw contour lines using colours.
    * \param markers Draw markers (as opposed to lines).
    * \param plotDriftLines Draw drift lines together with the isochrons.
    */
  void PlotIsochrons(const double tstep,
                     const std::vector<std::array<double, 3> >& points,
                     const bool rev = false,
                     const bool colour = false, const bool markers = false,
                     const bool plotDriftLines = true);

  /// Request electron drift lines with negative (default) or 
  /// positive charge.
  void DriftElectrons(const bool positive = false) { 
    m_particle = Particle::Electron;
    m_positive = positive;
  }
  /// Request ion drift lines with positive (default) or negative charge.
  void DriftIons(const bool negative = false) { 
    m_particle = Particle::Ion;
    m_positive = !negative;
  }
  /// Sort (or not) the points on a contour line (default: sorting is done).
  void EnableSorting(const bool on = true) { m_sortContours = on; }
  /// Check (or not) that drift-lines do not cross isochrons 
  /// (default: check is done).
  void CheckCrossings(const bool on = true) { m_checkCrossings = on; }
  /// Set the aspect ratio above which an isochron is considered
  /// linear (as opposed to circular).
  void SetAspectRatioSwitch(const double ar);
  /// Fractional distance between two points for closing a 
  /// circular isochron (default: 0.2).
  void SetLoopThreshold(const double thr);
  /// Fractional distance over which isochron segments are connected
  /// (default: 0.2).
  void SetConnectionThreshold(const double thr);

 private:
  Sensor* m_sensor = nullptr;
  ComponentBase* m_component = nullptr;

  // Projection for viewing
  double m_proj[3][3];
  double m_plane[4];
  char m_xLabel[50], m_yLabel[50], m_description[50];

  // Plot area
  bool m_hasUserArea = false;
  double m_xmin = -1., m_ymin = -1.;
  double m_xmax = 1., m_ymax = 1.;

  enum class Particle { Electron = 0, Ion };
  // Type of particle to be used for computing drift lines.
  Particle m_particle = Particle::Electron;
  bool m_positive = false;

  short m_markerStyle = 5;
  short m_lineStyle = 2;

  bool m_sortContours = true;
  double m_aspectRatio = 3.;
  double m_loopThreshold = 0.2;
  double m_connectionThreshold = 0.2;
  bool m_checkCrossings = true;

  void SetupCanvas();
  bool Range();
  void Labels();

  void ComputeDriftLines(const double tstep,
      const std::vector<std::array<double, 3> >& points,
      std::vector<std::vector<std::array<double, 3> > >& driftLines,
      std::vector<std::array<double, 3> >& startPoints, 
      std::vector<std::array<double, 3> >& endPoints, 
      std::vector<int>& statusCodes, const bool rev = false);
  void SortContour(
      std::vector<std::pair<std::array<double, 4>, unsigned int> >& contour,
      bool& circle);

};
}
#endif
