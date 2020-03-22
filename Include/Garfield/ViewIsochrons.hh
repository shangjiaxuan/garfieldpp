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

  void PlotIsochrons(const double tstep,
                     const std::vector<std::array<double, 3> >& points,
                     const bool rev = false,
                     const bool colour = false, const bool markers = false,
                     const bool plotDriftLines = true);

  void DriftElectrons(const bool positive = false) { 
    m_particle = Particle::Electron;
    m_positive = positive;
  }
  void DriftIons(const bool negative = false) { 
    m_particle = Particle::Ion;
    m_positive = !negative;
  }

  void EnableSorting(const bool on = true) { m_sortContours = on; }
  void CheckCrossings(const bool on = true) { m_checkCrossings = on; }
  void SetAspectRatioSwitch(const double ar);
  void SetLoopThreshold(const double thr);
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
  std::array<double, 3> PLACO3(const double x1, const double y1, const double z1);

};
}
#endif
