#ifndef G_VIEW_GEOMETRY
#define G_VIEW_GEOMETRY

#include <memory>
#include <string>
#include <vector>

#include <TGeoManager.h>

#include "ViewBase.hh"

namespace Garfield {

class GeometrySimple;

/// Visualize a geometry defined using the "native" shapes.

class ViewGeometry : public ViewBase {
 public:
  /// Constructor.
  ViewGeometry();
  /// Destructor.
  ~ViewGeometry();

  /// Set the geometry to be drawn.
  void SetGeometry(GeometrySimple* geo);
  /// Draw the geometry.
  void Plot();

 private:
  GeometrySimple* m_geometry = nullptr;

  std::vector<TGeoVolume*> m_volumes;
  std::vector<TGeoMedium*> m_media;

  std::unique_ptr<TGeoManager> m_geoManager;

  void Reset();
};
}
#endif
