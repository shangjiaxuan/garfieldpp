#ifndef G_VIEW_FE_MESH
#define G_VIEW_FE_MESH

#include <memory>
#include <string>
#ifndef __CINT__
#include <map>
#endif

#include <TArrayD.h>
#include <TGaxis.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>

#include "ComponentCST.hh"
#include "ComponentFieldMap.hh"
#include "ViewBase.hh"
#include "ViewDrift.hh"

namespace Garfield {

/// Draw the mesh of a field-map component.

class ViewFEMesh : public ViewBase {
 public:
  /// Constructor.
  ViewFEMesh();
  /// Destructor.
  ~ViewFEMesh() = default;

  /// Set the component from which to retrieve the mesh and field.
  void SetComponent(ComponentFieldMap* comp);

  /// Set area to be plotted to the default.
  void SetArea();
  /// Set area to be plotted explicitly.
  void SetArea(const double xmin, const double ymin, const double zmin,
               const double xmax, const double ymax, const double zmax);
  /// Reset the projection plane.
  void SetDefaultProjection();
  /// Set the projection plane.
  void SetPlane(double fx, double fy, double fz, double x0, double y0,
                double z0);
  /// Set the projection plane specifying hint for in-plane x axis.
  void SetPlane(double fx, double fy, double fz, double x0, double y0,
                double z0, double hx, double hy, double hz);

  // Axes
  void SetXaxis(TGaxis* ax);
  void SetYaxis(TGaxis* ay);
  void SetXaxisTitle(const char* xtitle);
  void SetYaxisTitle(const char* ytitle);
  void EnableAxes() { m_drawAxes = true; }
  void DisableAxes() { m_drawAxes = false; }

  /// Plot method to be called by user
  bool Plot();

  /// Element fill switch; 2D only, set false for wireframe mesh
  void SetFillMesh(const bool f) { m_fillMesh = f; }

  /// Display intersection of projection plane with viewing area
  void SetDrawViewRegion(bool do_draw) { m_drawViewRegion = do_draw; }
  bool GetDrawViewRegion(void) const { return m_drawViewRegion; }

  /// Associate a color with each element material map ID;
  /// Uses ROOT color numberings
  void SetColor(int matID, int colorID) { m_colorMap[matID] = colorID; }
  void SetFillColor(int matID, int colorID) {
    m_colorMap_fill[matID] = colorID;
  }

  /// Set the optional associated ViewDrift
  void SetViewDrift(ViewDrift* vd) { m_viewDrift = vd; }

  /// Show filled mesh elements
  void SetFillMeshWithBorders() {
    m_plotMeshBorders = true;
    m_fillMesh = true;
  }

  /// Create a default set of custom-made axes.
  void CreateDefaultAxes();

  /// Disable a material so that its mesh cells are not drawn
  void DisableMaterial(int materialID) {
    m_disabledMaterial[materialID] = true;
  }

 private:
  std::string m_label = "Mesh";

  // Options
  bool m_fillMesh = false;

  // Viewing plane (plane normal is stored in m_proj[2])
  double m_proj[3][3];
  double m_dist;

  // Box dimensions
  bool m_userBox = false;
  double m_xMinBox = -1., m_xMaxBox = 1.;
  double m_yMinBox = -1., m_yMaxBox = 1.;
  double m_zMinBox = -1., m_zMaxBox = 1.;

  // Intersection of viewing plane with plotted area in planar coordinates
  bool m_drawViewRegion = false;
  std::vector<TPolyLine> m_viewRegionLines;

  // Viewing plane dimensions
  double m_xMinPlot = -1., m_xMaxPlot = 1.;
  double m_yMinPlot = -1., m_yMaxPlot = 1.;

  // The field map object
  ComponentFieldMap* m_component = nullptr;

  // Optional associated ViewDrift object
  ViewDrift* m_viewDrift = nullptr;
  bool m_plotMeshBorders = false;

  // Axes
  TGaxis* m_xaxis = nullptr;
  TGaxis* m_yaxis = nullptr;
  std::unique_ptr<TH2D> m_axes;
  bool m_drawAxes = false;

  // The mesh, stored as a vector of TPolyLine(3D) objects
  std::vector<TPolyLine> m_mesh;
  std::vector<TPolyLine> m_driftLines;

// The color map
#ifndef __CINT__
  std::map<int, int> m_colorMap;
  std::map<int, int> m_colorMap_fill;

  // Disabled materials -> not shown in the mesh view
  std::map<int, bool> m_disabledMaterial;
#endif
  // Element plotting methods
  void DrawElements();
  void DrawCST(ComponentCST* componentCST);

  /// Return true if the specified point is in the view region.
  bool InView(const double x, const double y) const;

  bool LinesCrossed(double x1, double y1, double x2, double y2, double u1,
                    double v1, double u2, double v2, double& xc,
                    double& yc) const;
  std::string CreateAxisTitle(const double* norm) const;
  bool IntersectPlaneArea(void);
  bool PlaneVector(double& x, double& y, double& z) const;
  bool OnLine(double x1, double y1, double x2, double y2, double u,
              double v) const;
  void RemoveCrossings(std::vector<double>& x, std::vector<double>& y);
  bool PlaneCut(double x1, double y1, double z1, double x2, double y2,
                double z2, TMatrixD& xMat);
  bool PlaneCoords(double x, double y, double z, const TMatrixD& projMat,
                   TMatrixD& xMat);
  void ClipToView(std::vector<double>& px, std::vector<double>& py,
                  std::vector<double>& cx, std::vector<double>& cy);
  bool IsInPolygon(double x, double y, std::vector<double>& px,
                   std::vector<double>& py, bool& edge) const;

  // Plot method to be called by Plot() for CST cubic elements
  // available are "xy", "yz" and "xz"
};
}  // namespace Garfield
#endif
