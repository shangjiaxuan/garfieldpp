#ifndef G_VIEW_CELL
#define G_VIEW_CELL

#include <memory>
#include <string>

#include <TGeoManager.h>

#include "ViewBase.hh"

namespace Garfield {

class ComponentAnalyticField;
class ComponentNeBem2d;

/// Visualize the "cell" defined in an analytic-field component.

class ViewCell : public ViewBase {
 public:
  /// Constructor
  ViewCell();
  /// Destructor
  ~ViewCell() = default;

  /// Set the component for which to draw the cell geometry.
  void SetComponent(ComponentAnalyticField* comp);
  void SetComponent(ComponentNeBem2d* comp);

  /// Set the plot range explicitly.
  void SetArea(const double xmin, const double ymin, const double zmin,
               const double xmax, const double ymax, const double zmax);
  ///  Take the plot range from the bounding box of the component class.
  void SetArea() { m_hasUserArea = false; }

  /// Make a two-dimensional drawing of the cell geometry.
  void Plot2d();
  /// Make a three-dimensional drawing of the cell geometry (using TGeo).
  void Plot3d();

  /// Visualize wirers using markers or as a circle with the actual wire radius.
  /// The default is markers.
  void EnableWireMarkers(const bool on = true) { m_useWireMarker = on; }
  void DisableWireMarkers() { EnableWireMarkers(false); }

 private:
  // Options
  bool m_useWireMarker = true;

  std::string m_label = "Cell Layout";

  // Box dimensions
  bool m_hasUserArea = false;
  double m_xMin = -1., m_yMin = -1., m_zMin = -1.;
  double m_xMax = 1., m_yMax = 1., m_zMax = 1.;

  ComponentAnalyticField* m_component = nullptr;
  ComponentNeBem2d* m_nebem = nullptr;

  // 3D geometry.
  std::unique_ptr<TGeoManager> m_geo;

  bool Plot(const bool use3d);
  // Draw a wire in 2D.
  void PlotWire(const double x, const double y, const double d, const int type);
  // Draw a wire in 3D.
  void PlotWire(const double x, const double y, const double d, const int type,
                const double lz);
  // Draw a tube in 2D.
  void PlotTube(const double x0, const double y0, const double r, const int n);
  // Draw a tube in 3D.
  void PlotTube(const double x0, const double y0,
                const double r1, const double r2, const int n,
                const double lz); 
  // Draw a plane in 2D.
  void PlotPlane(const double x0, const double y0, 
                 const double x1, const double y1);
  // Draw a plane in 3D.
  void PlotPlane(const double dx, const double dy, const double dz,
                 const double x0, const double y0);
  // Draw a neBEM 2D layout.
  bool PlotNeBem(const bool use3d);
  // Setup the TGeoManager. 
  void SetupGeo(const double dx, const double dy, const double dz);
};
}
#endif
