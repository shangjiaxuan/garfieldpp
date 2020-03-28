#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Garfield/ComponentCST.hh"
#include "Garfield/ComponentFieldMap.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"
#include "Garfield/ViewFEMesh.hh"

namespace Garfield {

ViewFEMesh::ViewFEMesh() : ViewBase("ViewFEMesh") {
  SetDefaultProjection();

  // Create a blank histogram for the axes.
  m_axes.reset(new TH2D());
  m_axes->SetStats(false);
  m_axes->GetXaxis()->SetTitle("x");
  m_axes->GetYaxis()->SetTitle("y");
}

void ViewFEMesh::SetDefaultProjection() {
  // Default projection: x-y at z=0
  m_proj[0][0] = 1;
  m_proj[0][1] = 0;
  m_proj[0][2] = 0;
  m_proj[1][0] = 0;
  m_proj[1][1] = 1;
  m_proj[1][2] = 0;
  m_proj[2][0] = 0;
  m_proj[2][1] = 0;
  m_proj[2][2] = 1;
  // Plane distance to (0,0,0)
  m_dist = 0;
}

void ViewFEMesh::SetComponent(ComponentFieldMap* comp) {
  if (!comp) {
    std::cerr << m_className << "::SetComponent: Null pointer.\n";
    return;
  }

  m_component = comp;
}

void ViewFEMesh::SetArea(double xmin, double ymin, double zmin, double xmax,
                         double ymax, double zmax) {
  // Check range, assign if non-null
  if (xmin == xmax || ymin == ymax || zmin == zmax) {
    std::cerr << m_className << "::SetArea: Null area range not permitted.\n";
    return;
  }
  m_xMin = std::min(xmin, xmax);
  m_yMin = std::min(ymin, ymax);
  m_zMin = std::min(zmin, zmax);
  m_xMax = std::max(xmin, xmax);
  m_yMax = std::max(ymin, ymax);
  m_zMax = std::max(zmin, zmax);

  m_hasUserArea = true;
  IntersectPlaneArea();
}

void ViewFEMesh::SetArea() { m_hasUserArea = false; }

// The plotting functionality here is ported from Garfield
//  with some inclusion of code from ViewCell.cc
bool ViewFEMesh::Plot() {
  if (!m_component) {
    std::cerr << m_className << "::Plot: Component is not defined.\n";
    return false;
  }

  double pmin = 0., pmax = 0.;
  if (!m_component->GetVoltageRange(pmin, pmax)) {
    std::cerr << m_className << "::Plot: Component is not ready.\n";
    return false;
  }

  // Get the bounding box.
  if (!m_hasUserArea) {
    std::cerr << m_className << "::Plot:\n"
              << "    Bounding box cannot be determined. Call SetArea first.\n";
    return false;
  }

  if (m_viewRegionLines.empty()) {
    std::cerr << m_className << "::Plot:\n"
              << "    Empty view. Make sure the viewing plane (SetPlane)\n"
              << "    intersects with the bounding box.\n";
    return false;
  }

  // Set up a canvas if one does not already exist.
  if (!m_canvas) {
    m_canvas = new TCanvas();
    m_canvas->SetTitle(m_label.c_str());
    if (m_hasExternalCanvas) m_hasExternalCanvas = false;
  }
  m_canvas->Range(m_xPlaneMin, m_yPlaneMin, m_xPlaneMax, m_yPlaneMax);

  // Plot the elements
  ComponentCST* componentCST = dynamic_cast<ComponentCST*>(m_component);
  if (componentCST) {
    std::cout << m_className << "::Plot: CST component. Calling DrawCST.\n";
    DrawCST(componentCST);
  } else {
    DrawElements();
  }
  m_canvas->Update();

  return true;
}

// Set the projection plane: modified from ViewField.cc
// to match functionality of Garfield
void ViewFEMesh::SetPlane(const double fx, const double fy, const double fz,
                          const double x0, const double y0, const double z0) {
  if (fy * fy + fz * fz > 0) {
    SetPlane(fx, fy, fz, x0, y0, z0, 1, 0, 0);
  } else {
    SetPlane(fx, fy, fz, x0, y0, z0, 0, 1, 0);
  }
}

// Set the projection plane specifying hint for in-plane x axis.
void ViewFEMesh::SetPlane(const double fx, const double fy, const double fz,
                          const double x0, const double y0, const double z0,
                          const double hx, const double hy, const double hz) {
  // Calculate 2 in-plane vectors for the normal vector
  double fnorm = sqrt(fx * fx + fy * fy + fz * fz);
  if (fnorm < Small) {
    std::cout << m_className << "::SetPlane:\n"
              << "    Normal vector has zero norm. No new projection set.\n";
    return;
  }
  double dist = (fx * x0 + fy * y0 + fz * z0) / fnorm;
  // Store the plane description
  m_proj[2][0] = fx / fnorm;
  m_proj[2][1] = fy / fnorm;
  m_proj[2][2] = fz / fnorm;
  m_dist = dist;

  double xx = hx, xy = hy, xz = hz;
  PlaneVector(xx, xy, xz);
  double vecx_norm = std::sqrt(xx * xx + xy * xy + xz * xz);
  if (vecx_norm < 1.0e-10) {  
    // Wrong in-plane x hint (close to norm).
    if (fy * fy + fz * fz > 0) {
      // Taking global x as in-plane x hint.
      xx = 1;
      xy = 0;
      xz = 0;  
    } else {
      // Taking global y as in-plane x hint.
      xx = 0;
      xy = 1;
      xz = 0;  
    }
    PlaneVector(xx, xy, xz);
    vecx_norm = std::sqrt(xx * xx + xy * xy + xz * xz);
  }
  m_proj[0][0] = xx / vecx_norm;
  m_proj[0][1] = xy / vecx_norm;
  m_proj[0][2] = xz / vecx_norm;
  // in-plane y === m_proj[1] = cross product [z,x];
  m_proj[1][0] = m_proj[2][1] * m_proj[0][2] - m_proj[2][2] * m_proj[0][1];
  m_proj[1][1] = m_proj[2][2] * m_proj[0][0] - m_proj[2][0] * m_proj[0][2];
  m_proj[1][2] = m_proj[2][0] * m_proj[0][1] - m_proj[2][1] * m_proj[0][0];

  IntersectPlaneArea();
}

// Set the x-axis.
void ViewFEMesh::SetXaxis(TGaxis* ax) { m_xaxis = ax; }

// Set the y-axis.
void ViewFEMesh::SetYaxis(TGaxis* ay) { m_yaxis = ay; }

// Set the x-axis title.
void ViewFEMesh::SetXaxisTitle(const char* xtitle) {
  if (m_axes) m_axes->GetXaxis()->SetTitle(xtitle);
}

// Set the y-axis title.
void ViewFEMesh::SetYaxisTitle(const char* ytitle) {
  if (m_axes) m_axes->GetYaxis()->SetTitle(ytitle);
}

// Create default axes
void ViewFEMesh::CreateDefaultAxes() {
  // Create a new x and y axis.
  const double dx = std::abs(m_xPlaneMax - m_xPlaneMin) * 0.1;
  const double dy = std::abs(m_yPlaneMax - m_yPlaneMin) * 0.1;
  const double x0 = m_xPlaneMin + dx;
  const double y0 = m_yPlaneMin + dy;
  const double x1 = m_xPlaneMax - dx;
  const double y1 = m_yPlaneMax - dy;
  m_xaxis = new TGaxis(x0, y0, x1, y0, x0, x1, 2405, "x");
  m_yaxis = new TGaxis(x0, y0, x0, y1, y0, y1, 2405, "y");

  // Label sizes
  m_xaxis->SetLabelSize(0.025);
  m_yaxis->SetLabelSize(0.025);

  // Titles
  m_xaxis->SetTitleSize(0.03);
  std::string name = CreateAxisTitle(m_proj[0]);
  m_xaxis->SetTitle(name.c_str());
  m_yaxis->SetTitleSize(0.03);
  name = CreateAxisTitle(m_proj[1]);
  m_yaxis->SetTitle(name.c_str());
}

std::string ViewFEMesh::CreateAxisTitle(const double* norm) const {
  std::string name;
  for (int i = 0; i < 3; ++i) {
    std::stringstream num;
    num << std::setprecision(6);
    double value = std::abs(norm[i]);
    if (value > 1.0e-10) num << value;
    std::string temp = num.str();
    if (temp.empty()) continue;
    if (norm[i] > 0) {
      name += (name.empty() ? "" : " + ");
    } else {
      name += (name.empty() ? "#minus" : " #minus ");
    } 
    // No need to write "1.00000*y [cm], use "y [cm]"
    name += (temp == "1" ? "" : temp + " ");  
    switch (i) {
      case 0: 
        name += "x";
        break;
      case 1: 
        name += "y";
        break;
      case 2: 
        name += "z";
        break;
    }
  }
  name += " [cm]";
  return name;
}

// Use ROOT plotting functions to draw the mesh elements on the canvas.
// General methodology ported from Garfield
void ViewFEMesh::DrawElements() {
  // Get the map boundaries from the component.
  double mapxmax = m_component->m_mapmax[0];
  double mapxmin = m_component->m_mapmin[0];
  double mapymax = m_component->m_mapmax[1];
  double mapymin = m_component->m_mapmin[1];
  double mapzmax = m_component->m_mapmax[2];
  double mapzmin = m_component->m_mapmin[2];

  // Get the periodicities.
  double sx = mapxmax - mapxmin;
  double sy = mapymax - mapymin;
  double sz = mapzmax - mapzmin;
  const bool perX =
      m_component->m_periodic[0] || m_component->m_mirrorPeriodic[0];
  const bool perY =
      m_component->m_periodic[1] || m_component->m_mirrorPeriodic[1];
  const bool perZ =
      m_component->m_periodic[2] || m_component->m_mirrorPeriodic[2];

  // Clear the meshes and drift line lists.
  m_mesh.clear();
  m_driftLines.clear();

  // Prepare the final projection matrix (the transpose of the 2D array
  // "project").
  TArrayD dataProj(9);
  dataProj[0] = m_proj[0][0];
  dataProj[1] = m_proj[1][0];
  dataProj[2] = m_proj[2][0];
  dataProj[3] = m_proj[0][1];
  dataProj[4] = m_proj[1][1];
  dataProj[5] = m_proj[2][1];
  dataProj[6] = m_proj[0][2];
  dataProj[7] = m_proj[1][2];
  dataProj[8] = m_proj[2][2];
  TMatrixD projMat(3, 3, dataProj.GetArray());

  // Calculate the determinant of the projection matrix.
  double projDet =
      projMat(0, 0) *
          (projMat(1, 1) * projMat(2, 2) - projMat(1, 2) * projMat(2, 1)) -
      projMat(0, 1) *
          (projMat(1, 0) * projMat(2, 2) - projMat(1, 2) * projMat(2, 0)) +
      projMat(0, 2) *
          (projMat(1, 0) * projMat(2, 1) - projMat(1, 1) * projMat(2, 0));

  // Calculate the inverse of the projection matrix for
  // calculating coordinates in the viewing plane.
  if (projDet != 0) {
    projMat.Invert();
  } else {
    std::cerr << m_className << "::DrawElements:\n";
    std::cerr << "    Projection matrix is not invertible.\n";
    std::cerr << "    Finite element mesh will not be drawn.\n";
  }

  // Get the plane information.
  double fx = m_proj[2][0];
  double fy = m_proj[2][1];
  double fz = m_proj[2][2];
  double dist = m_dist;

  // Construct two empty single-column matrices for use as coordinate vectors.
  TMatrixD xMat(3, 1);

  // Determine the number of periods present in the cell.
  const int nMinX = perX ? int(m_xMin / sx) - 1 : 0;
  const int nMaxX = perX ? int(m_xMax / sx) + 1 : 0;
  const int nMinY = perY ? int(m_yMin / sy) - 1 : 0;
  const int nMaxY = perY ? int(m_yMax / sy) + 1 : 0;
  const int nMinZ = perZ ? int(m_zMin / sz) - 1 : 0;
  const int nMaxZ = perZ ? int(m_zMax / sz) + 1 : 0;

  // Loop over all elements.
  for (const auto& element : m_component->elements) {
    const auto mat = element.matmap;
    // Do not plot the drift medium.
    if (m_component->materials[mat].driftmedium && !(m_plotMeshBorders)) {
      continue;
    }
    // Do not create Polygons for disabled materials
    if (m_disabledMaterial[mat]) continue;
    // -- Tetrahedral elements

    // Coordinates of vertices
    double vx1, vy1, vz1;
    double vx2, vy2, vz2;
    double vx3, vy3, vz3;
    double vx4, vy4, vz4;

    // Get the color for this element (default to 1).
    int colorID = m_colorMap.count(mat);
    if (colorID != 0)
      colorID = m_colorMap[mat];
    else
      colorID = 1;

    // Get the fill color for this element (default colorID).
    int colorID_fill = m_colorMap_fill.count(mat);
    if (colorID_fill != 0)
      colorID_fill = m_colorMap_fill[mat];
    else
      colorID_fill = colorID;

    const auto& n0 = m_component->nodes[element.emap[0]];
    const auto& n1 = m_component->nodes[element.emap[1]];
    const auto& n2 = m_component->nodes[element.emap[2]];
    const auto& n3 = m_component->nodes[element.emap[3]];
    // Loop over the periodicities in x.
    for (int nx = nMinX; nx <= nMaxX; nx++) {
      // Determine the x-coordinates of the tetrahedral vertices.
      if (m_component->m_mirrorPeriodic[0] && nx != 2 * (nx / 2)) {
        vx1 = mapxmin + (mapxmax - n0.x) + sx * nx;
        vx2 = mapxmin + (mapxmax - n1.x) + sx * nx;
        vx3 = mapxmin + (mapxmax - n2.x) + sx * nx;
        vx4 = mapxmin + (mapxmax - n3.x) + sx * nx;
      } else {
        vx1 = n0.x + sx * nx;
        vx2 = n1.x + sx * nx;
        vx3 = n2.x + sx * nx;
        vx4 = n3.x + sx * nx;
      }

      // Loop over the periodicities in y.
      for (int ny = nMinY; ny <= nMaxY; ny++) {
        // Determine the y-coordinates of the tetrahedral vertices.
        if (m_component->m_mirrorPeriodic[1] && ny != 2 * (ny / 2)) {
          vy1 = mapymin + (mapymax - n0.y) + sy * ny;
          vy2 = mapymin + (mapymax - n1.y) + sy * ny;
          vy3 = mapymin + (mapymax - n2.y) + sy * ny;
          vy4 = mapymin + (mapymax - n3.y) + sy * ny;
        } else {
          vy1 = n0.y + sy * ny;
          vy2 = n1.y + sy * ny;
          vy3 = n2.y + sy * ny;
          vy4 = n3.y + sy * ny;
        }

        // Loop over the periodicities in z.
        for (int nz = nMinZ; nz <= nMaxZ; nz++) {
          // Determine the z-coordinates of the tetrahedral vertices.
          if (m_component->m_mirrorPeriodic[2] && nz != 2 * (nz / 2)) {
            vz1 = mapzmin + (mapzmax - n0.z) + sz * nz;
            vz2 = mapzmin + (mapzmax - n1.z) + sz * nz;
            vz3 = mapzmin + (mapzmax - n2.z) + sz * nz;
            vz4 = mapzmin + (mapzmax - n3.z) + sz * nz;
          } else {
            vz1 = n0.z + sz * nz;
            vz2 = n1.z + sz * nz;
            vz3 = n2.z + sz * nz;
            vz4 = n3.z + sz * nz;
          }

          // Store the x and y coordinates of the relevant mesh vertices.
          std::vector<double> vX;
          std::vector<double> vY;

          // Value used to determine whether a vertex is in the plane.
          const double pcf = std::max(
              {std::abs(vx1), std::abs(vy1), std::abs(vz1), std::abs(fx),
               std::abs(fy), std::abs(fz), std::abs(dist)});
          const double tol = 1.e-4 * pcf;
          // First isolate the vertices that are in the viewing plane.
          bool in1 = (std::abs(fx * vx1 + fy * vy1 + fz * vz1 - dist) < tol);
          bool in2 = (std::abs(fx * vx2 + fy * vy2 + fz * vz2 - dist) < tol);
          bool in3 = (std::abs(fx * vx3 + fy * vy3 + fz * vz3 - dist) < tol);
          bool in4 = (std::abs(fx * vx4 + fy * vy4 + fz * vz4 - dist) < tol);

          // Calculate the planar coordinates for those edges that are in the
          // plane.
          if (in1) {
            PlaneCoords(vx1, vy1, vz1, projMat, xMat);
            vX.push_back(xMat(0, 0));
            vY.push_back(xMat(1, 0));
          }
          if (in2) {
            PlaneCoords(vx2, vy2, vz2, projMat, xMat);
            vX.push_back(xMat(0, 0));
            vY.push_back(xMat(1, 0));
          }
          if (in3) {
            PlaneCoords(vx3, vy3, vz3, projMat, xMat);
            vX.push_back(xMat(0, 0));
            vY.push_back(xMat(1, 0));
          }
          if (in4) {
            PlaneCoords(vx4, vy4, vz4, projMat, xMat);
            vX.push_back(xMat(0, 0));
            vY.push_back(xMat(1, 0));
          }

          // Cut the sides that are not in the plane.
          if (!(in1 || in2)) {
            if (PlaneCut(vx1, vy1, vz1, vx2, vy2, vz2, xMat)) {
              vX.push_back(xMat(0, 0));
              vY.push_back(xMat(1, 0));
            }
          }
          if (!(in1 || in3)) {
            if (PlaneCut(vx1, vy1, vz1, vx3, vy3, vz3, xMat)) {
              vX.push_back(xMat(0, 0));
              vY.push_back(xMat(1, 0));
            }
          }
          if (!(in1 || in4)) {
            if (PlaneCut(vx1, vy1, vz1, vx4, vy4, vz4, xMat)) {
              vX.push_back(xMat(0, 0));
              vY.push_back(xMat(1, 0));
            }
          }
          if (!(in2 || in3)) {
            if (PlaneCut(vx2, vy2, vz2, vx3, vy3, vz3, xMat)) {
              vX.push_back(xMat(0, 0));
              vY.push_back(xMat(1, 0));
            }
          }
          if (!(in2 || in4)) {
            if (PlaneCut(vx2, vy2, vz2, vx4, vy4, vz4, xMat)) {
              vX.push_back(xMat(0, 0));
              vY.push_back(xMat(1, 0));
            }
          }
          if (!(in3 || in4)) {
            if (PlaneCut(vx3, vy3, vz3, vx4, vy4, vz4, xMat)) {
              vX.push_back(xMat(0, 0));
              vY.push_back(xMat(1, 0));
            }
          }
          if (vX.size() < 3) continue;
          // Create a convex TPolyLine object connecting the points.

          // Eliminate crossings of the polygon lines
          // (known as "butterflies" in Garfield).
          RemoveCrossings(vX, vY);

          // Create vectors to store the clipped polygon.
          std::vector<double> cX;
          std::vector<double> cY;

          // Clip the polygon to the view area.
          ClipToView(vX, vY, cX, cY);

          // If we have more than 2 vertices, add the polygon.
          if (cX.size() <= 2) continue;

          // Again eliminate crossings of the polygon lines.
          RemoveCrossings(cX, cY);

          // Create the TPolyLine.
          TPolyLine poly;
          poly.SetLineColor(colorID);
          poly.SetFillColor(colorID_fill);
          poly.SetLineWidth(3);
          // Add all of the points.
          const auto nPoints = cX.size();
          for (size_t pt = 0; pt < nPoints; ++pt) {
            poly.SetPoint(pt, cX[pt], cY[pt]);
          }

          // Add the polygon to the mesh.
          m_mesh.push_back(std::move(poly));
        }  // end z-periodicity loop
      }    // end y-periodicity loop
    }      // end x-periodicity loop
  }        // end loop over elements

  // If we have an associated ViewDrift, plot projections of the drift lines.
  if (m_viewDrift) {
    for (const auto& dline : m_viewDrift->m_driftLines) {
      // Create a TPolyLine that is a 2D projection of the original.
      TPolyLine poly;
      if (dline.second == ViewDrift::Particle::Electron) {
        poly.SetLineColor(kOrange - 3);
      } else {
        poly.SetLineColor(kRed + 1);
      }
      int polyPts = 0;
      for (const auto& point : dline.first) {
        // Project this point onto the plane.
        PlaneCoords(point[0], point[1], point[2], projMat, xMat);
        // Add this point if it is within the view.
        if (InView(xMat(0, 0), xMat(1, 0))) {
          poly.SetPoint(polyPts, xMat(0, 0), xMat(1, 0));
          polyPts++;
        }
      }  // end loop over points

      // Add the drift line to the list.
      m_driftLines.push_back(std::move(poly));

    }  // end loop over drift lines
  }    // end if(m_viewDrift != 0)

  // Call the ROOT draw methods to plot the elements.
  m_canvas->cd();

  // Draw default axes by using a blank 2D histogram.
  if (!m_xaxis && !m_yaxis && m_drawAxes) {
    std::string name = CreateAxisTitle(m_proj[0]);
    m_axes->GetXaxis()->SetTitle(name.c_str());
    m_axes->GetXaxis()->SetLimits(m_xPlaneMin, m_xPlaneMax);
    name = CreateAxisTitle(m_proj[1]);
    m_axes->GetYaxis()->SetTitle(name.c_str());
    m_axes->GetYaxis()->SetTitleOffset(1.6);
    m_axes->GetYaxis()->SetLimits(m_yPlaneMin, m_yPlaneMax);
    m_axes->Draw();
  }

  // Draw custom axes.
  if (m_xaxis && m_drawAxes) m_xaxis->Draw();
  if (m_yaxis && m_drawAxes) m_yaxis->Draw();

  // Draw the mesh on the canvas.
  for (auto& m : m_mesh) {
    if (m_plotMeshBorders || !m_fillMesh) m.Draw("same");
    if (m_fillMesh) m.Draw("f:same");
  }
  // Draw the drift lines on the view.
  for (auto& dline : m_driftLines) dline.Draw("same");

  if (m_drawViewRegion)
    for (auto& m : m_viewRegionLines) m.Draw("same");

  // Draw axes again so they are on top
  gPad->RedrawAxis("g");
}

void ViewFEMesh::DrawCST(ComponentCST* componentCST) {
  /*The method is based on ViewFEMesh::Draw, thus the first part is copied from
   * there.
   * At the moment only x-y, x-z, and y-z are available due to the simple
   * implementation.
   * The advantage of this method is that there is no element loop and thus it
   * is much
   * faster.
   */
  // Get the map boundaries from the component
  double mapxmax = m_component->m_mapmax[0];
  double mapxmin = m_component->m_mapmin[0];
  double mapymax = m_component->m_mapmax[1];
  double mapymin = m_component->m_mapmin[1];
  double mapzmax = m_component->m_mapmax[2];
  double mapzmin = m_component->m_mapmin[2];

  // Get the periodicities.
  double sx = mapxmax - mapxmin;
  double sy = mapymax - mapymin;
  double sz = mapzmax - mapzmin;
  const bool perX =
      m_component->m_periodic[0] || m_component->m_mirrorPeriodic[0];
  const bool perY =
      m_component->m_periodic[1] || m_component->m_mirrorPeriodic[1];
  const bool perZ =
      m_component->m_periodic[2] || m_component->m_mirrorPeriodic[2];

  // Clear the meshes and drift line lists
  m_mesh.clear();
  m_driftLines.clear();

  // Prepare the final projection matrix (the transpose of the 2D array
  // "project")
  TArrayD dataProj(9);
  dataProj[0] = m_proj[0][0];
  dataProj[1] = m_proj[1][0];
  dataProj[2] = m_proj[2][0];
  dataProj[3] = m_proj[0][1];
  dataProj[4] = m_proj[1][1];
  dataProj[5] = m_proj[2][1];
  dataProj[6] = m_proj[0][2];
  dataProj[7] = m_proj[1][2];
  dataProj[8] = m_proj[2][2];
  TMatrixD projMat(3, 3, dataProj.GetArray());

  // Calculate the determinant of the projection matrix
  double projDet =
      projMat(0, 0) *
          (projMat(1, 1) * projMat(2, 2) - projMat(1, 2) * projMat(2, 1)) -
      projMat(0, 1) *
          (projMat(1, 0) * projMat(2, 2) - projMat(1, 2) * projMat(2, 0)) +
      projMat(0, 2) *
          (projMat(1, 0) * projMat(2, 1) - projMat(1, 1) * projMat(2, 0));

  // Calculate the inverse of the projection matrix for
  // calculating coordinates in the viewing plane
  if (projDet != 0) {
    projMat.Invert();
  } else {
    std::cerr << m_className << "::DrawCST:\n";
    std::cerr << "    Projection matrix is not invertible.\n";
    std::cerr << "    Finite element mesh will not be drawn.\n";
  }

  // Construct two empty single-column matrices for use as coordinate vectors
  TMatrixD xMat(3, 1);

  // Determine the number of periods present in the cell.
  const int nMinX = perX ? int(m_xMin / sx) - 1 : 0;
  const int nMaxX = perX ? int(m_xMax / sx) + 1 : 0;
  const int nMinY = perY ? int(m_yMin / sy) - 1 : 0;
  const int nMaxY = perY ? int(m_yMax / sy) + 1 : 0;
  const int nMinZ = perZ ? int(m_zMin / sz) - 1 : 0;
  const int nMaxZ = perZ ? int(m_zMax / sz) + 1 : 0;

  int elem = 0;
  std::vector<PolygonInfo> elements;
  int nMinU = 0, nMaxU = 0, nMinV = 0, nMaxV = 0;
  double mapumin = 0., mapumax = 0., mapvmin = 0., mapvmax = 0.;
  double su = 0., sv = 0.;
  bool mirroru = false, mirrorv = false;
  double uMin, vMin, uMax, vMax;
  unsigned int n_x, n_y, n_z;
  componentCST->GetNumberOfMeshLines(n_x, n_y, n_z);
  double e_xmin, e_xmax, e_ymin, e_ymax, e_zmin, e_zmax;
  // xy view
  if (m_proj[2][0] == 0 && m_proj[2][1] == 0 && m_proj[2][2] == 1) {
    std::cout << m_className << "::DrawCST: Creating x-y mesh view.\n";
    ViewFEMesh::SetXaxisTitle("x [cm]");
    ViewFEMesh::SetYaxisTitle("y [cm]");
    // calculate the z position
    unsigned int i, j, z;
    if (!componentCST->Coordinate2Index(0, 0, m_dist * m_proj[2][2], i, j, z)) {
      std::cerr << "    Could not determine the position of the plane in "
                << "z direction.\n";
      return;
    }
    std::cout << "    The plane position in z direction is: "
              << m_dist * m_proj[2][2] << "\n";
    nMinU = nMinX;
    nMaxU = nMaxX;
    nMinV = nMinY;
    nMaxV = nMaxY;
    uMin = m_xMin;
    uMax = m_xMax;
    vMin = m_yMin;
    vMax = m_yMax;

    mapumin = mapxmin;
    mapumax = mapxmax;
    mapvmin = mapymin;
    mapvmax = mapymax;
    su = sx;
    sv = sy;
    mirroru = perX;
    mirrorv = perY;
    for (unsigned int y = 0; y < (n_y - 1); y++) {
      for (unsigned int x = 0; x < (n_x - 1); x++) {
        elem = componentCST->Index2Element(x, y, z);
        componentCST->GetElementBoundaries(elem, e_xmin, e_xmax, e_ymin, e_ymax,
                                           e_zmin, e_zmax);
        PolygonInfo tmp_info;
        tmp_info.element = elem;
        tmp_info.p1[0] = e_xmin;
        tmp_info.p2[0] = e_xmax;
        tmp_info.p3[0] = e_xmax;
        tmp_info.p4[0] = e_xmin;
        tmp_info.p1[1] = e_ymin;
        tmp_info.p2[1] = e_ymin;
        tmp_info.p3[1] = e_ymax;
        tmp_info.p4[1] = e_ymax;
        tmp_info.material = componentCST->GetElementMaterial(elem);
        elements.push_back(std::move(tmp_info));
      }
    }
    // xz-view
  } else if (m_proj[2][0] == 0 && m_proj[2][1] == -1 && m_proj[2][2] == 0) {
    std::cout << m_className << "::DrawCST: Creating x-z mesh view.\n";
    ViewFEMesh::SetXaxisTitle("x [cm]");
    ViewFEMesh::SetYaxisTitle("z [cm]");
    // calculate the y position
    unsigned int i = 0, j = 0, y = 0;
    if (!componentCST->Coordinate2Index(0, m_dist * m_proj[2][1], 0, i, y, j)) {
      std::cerr << "    Could not determine the position of the plane in "
                << "y direction.\n";
      return;
    }
    std::cout << "    The plane position in y direction is: "
              << m_dist * m_proj[2][1] << "\n";

    nMinU = nMinX;
    nMaxU = nMaxX;
    nMinV = nMinZ;
    nMaxV = nMaxZ;
    uMin = m_xMin;
    uMax = m_xMax;
    vMin = m_zMin;
    vMax = m_zMax;

    mapumin = mapxmin;
    mapumax = mapxmax;
    mapvmin = mapzmin;
    mapvmax = mapzmax;
    su = sx;
    sv = sz;
    mirroru = perX;
    mirrorv = perZ;
    for (unsigned int z = 0; z < (n_z - 1); z++) {
      for (unsigned int x = 0; x < (n_x - 1); x++) {
        elem = componentCST->Index2Element(x, y, z);
        componentCST->GetElementBoundaries(elem, e_xmin, e_xmax, e_ymin, e_ymax,
                                           e_zmin, e_zmax);
        PolygonInfo tmp_info;
        tmp_info.element = elem;
        tmp_info.p1[0] = e_xmin;
        tmp_info.p2[0] = e_xmax;
        tmp_info.p3[0] = e_xmax;
        tmp_info.p4[0] = e_xmin;
        tmp_info.p1[1] = e_zmin;
        tmp_info.p2[1] = e_zmin;
        tmp_info.p3[1] = e_zmax;
        tmp_info.p4[1] = e_zmax;
        tmp_info.material = componentCST->GetElementMaterial(elem);
        elements.push_back(std::move(tmp_info));
      }
    }

    // yz-view
  } else if (m_proj[2][0] == -1 && m_proj[2][1] == 0 && m_proj[2][2] == 0) {
    std::cout << m_className << "::DrawCST: Creating z-y mesh view.\n";
    ViewFEMesh::SetXaxisTitle("z [cm]");
    ViewFEMesh::SetYaxisTitle("y [cm]");
    // calculate the x position
    unsigned int i, j, x;
    if (!componentCST->Coordinate2Index(m_dist * m_proj[2][0], 0, 0, x, i,
                                        j)) {
      std::cerr << "    Could not determine the position of the plane in "
                << "x direction.\n";
      return;
    }
    std::cout << "    The plane position in x direction is: "
              << m_dist * m_proj[2][0] << "\n";
    nMinU = nMinZ;
    nMaxU = nMaxZ;
    nMinV = nMinY;
    nMaxV = nMaxY;
    uMin = m_yMin;
    uMax = m_yMax;
    vMin = m_zMin;
    vMax = m_zMax;

    mapumin = mapzmin;
    mapumax = mapzmax;
    mapvmin = mapymin;
    mapvmax = mapymax;
    su = sz;
    sv = sy;
    mirroru = perZ;
    mirrorv = perY;
    for (unsigned int z = 0; z < (n_z - 1); z++) {
      for (unsigned int y = 0; y < (n_y - 1); y++) {
        elem = componentCST->Index2Element(x, y, z);
        componentCST->GetElementBoundaries(elem, e_xmin, e_xmax, e_ymin, e_ymax,
                                           e_zmin, e_zmax);
        PolygonInfo tmp_info;
        tmp_info.element = elem;
        tmp_info.p1[0] = e_zmin;
        tmp_info.p2[0] = e_zmax;
        tmp_info.p3[0] = e_zmax;
        tmp_info.p4[0] = e_zmin;
        tmp_info.p1[1] = e_ymin;
        tmp_info.p2[1] = e_ymin;
        tmp_info.p3[1] = e_ymax;
        tmp_info.p4[1] = e_ymax;
        tmp_info.material = componentCST->GetElementMaterial(elem);
        // Add the polygon to the mesh
        elements.push_back(std::move(tmp_info));
      }
    }
  } else {
    std::cerr << m_className << "::DrawCST:\n";
    std::cerr << "    The given plane name is not known.\n";
    std::cerr << "    Please choose one of the following: xy, xz, yz.\n";
    return;
  }
  std::cout << m_className << "::DrawCST:\n";
  std::cout << "    Number of elements in the projection of the unit cell:"
            << elements.size() << std::endl;
  std::vector<PolygonInfo>::iterator it;
  std::vector<PolygonInfo>::iterator itend = elements.end();

  for (int nu = nMinU; nu <= nMaxU; nu++) {
    for (int nv = nMinV; nv <= nMaxV; nv++) {
      it = elements.begin();
      while (it != itend) {
        if (m_disabledMaterial[(*it).material]) {
          // do not create Polygons for disabled materials
          it++;
          continue;
        }
        int colorID = m_colorMap.count((*it).material);
        if (colorID != 0)
          colorID = m_colorMap[(*it).material];
        else
          colorID = 1;

        // Get the fill color for this element (default colorID)
        int colorID_fill = m_colorMap_fill.count((*it).material);
        if (colorID_fill != 0)
          colorID_fill = m_colorMap_fill[(*it).material];
        else
          colorID_fill = colorID;

        TPolyLine poly;
        poly.SetLineColor(colorID);
        poly.SetFillColor(colorID_fill);
        if (m_plotMeshBorders)
          poly.SetLineWidth(3);
        else
          poly.SetLineWidth(1);
        // Add 4 points of the square
        Double_t tmp_u[4], tmp_v[4];
        if (mirroru && nu != 2 * (nu / 2)) {
          // nu is odd
          tmp_u[0] = mapumin + (mapumax - (*it).p1[0]) + su * nu;
          tmp_u[1] = mapumin + (mapumax - (*it).p2[0]) + su * nu;
          tmp_u[2] = mapumin + (mapumax - (*it).p3[0]) + su * nu;
          tmp_u[3] = mapumin + (mapumax - (*it).p4[0]) + su * nu;
        } else {
          // nu is even
          tmp_u[0] = (*it).p1[0] + su * nu;
          tmp_u[1] = (*it).p2[0] + su * nu;
          tmp_u[2] = (*it).p3[0] + su * nu;
          tmp_u[3] = (*it).p4[0] + su * nu;
        }
        if (mirrorv && nv != 2 * (nv / 2)) {
          tmp_v[0] = mapvmin + (mapvmax - (*it).p1[1]) + sv * nv;
          tmp_v[1] = mapvmin + (mapvmax - (*it).p2[1]) + sv * nv;
          tmp_v[2] = mapvmin + (mapvmax - (*it).p3[1]) + sv * nv;
          tmp_v[3] = mapvmin + (mapvmax - (*it).p4[1]) + sv * nv;
        } else {
          tmp_v[0] = (*it).p1[1] + sv * nv;
          tmp_v[1] = (*it).p2[1] + sv * nv;
          tmp_v[2] = (*it).p3[1] + sv * nv;
          tmp_v[3] = (*it).p4[1] + sv * nv;
        }
        if (tmp_u[0] < uMin || tmp_u[1] > uMax || tmp_v[0] < vMin ||
            tmp_v[2] > vMax) {
          it++;
          continue;
        }
        poly.SetPoint(0, tmp_u[0], tmp_v[0]);
        poly.SetPoint(1, tmp_u[1], tmp_v[1]);
        poly.SetPoint(2, tmp_u[2], tmp_v[2]);
        poly.SetPoint(3, tmp_u[3], tmp_v[3]);
        // Add the polygon to the mesh
        m_mesh.push_back(std::move(poly));
        it++;
      }  // end element loop
    }    // end v-periodicity loop
  }      // end u-periodicity loop
  std::cout << m_className << "::PlotCST:\n"
            << "    Number of polygons to be drawn:" << m_mesh.size() << "\n";
  // If we have an associated ViewDrift, plot projections of the drift lines
  if (m_viewDrift) {
    for (const auto& dline : m_viewDrift->m_driftLines) {
      // Create a TPolyLine that is a 2D projection of the original
      TPolyLine poly;
      if (dline.second == ViewDrift::Particle::Electron) {
        poly.SetLineColor(kOrange - 3);
      } else {
        poly.SetLineColor(kRed + 1);
      }
      int polyPts = 0;
      for (const auto& point : dline.first) {
        // Project this point onto the plane.
        PlaneCoords(point[0], point[1], point[2], projMat, xMat);
        // Add this point if it is within the view
        if (xMat(0, 0) >= uMin && xMat(0, 0) <= uMax && xMat(1, 0) >= vMin &&
            xMat(1, 0) <= vMax) {
          poly.SetPoint(polyPts, xMat(0, 0), xMat(1, 0));
          polyPts++;
        }
      }  // end loop over points
      // Add the drift line to the list
      m_driftLines.push_back(std::move(poly));
    }  // end loop over drift lines
  }    // end if(m_viewDrift != 0)

  // Call the ROOT draw methods to plot the elements
  m_canvas->cd();
  // Draw default axes by using a blank 2D histogram.
  if (!m_xaxis && !m_yaxis && m_drawAxes) {
    m_axes->GetXaxis()->SetLimits(uMin, uMax);
    m_axes->GetYaxis()->SetLimits(vMin, vMax);
    m_axes->Draw();
  }
  // Draw custom axes.
  if (m_xaxis && m_drawAxes) m_xaxis->Draw("");
  if (m_yaxis && m_drawAxes) m_yaxis->Draw("");
  // Draw the mesh on the canvas
  for (auto& m : m_mesh) {
    if (m_plotMeshBorders || !m_fillMesh) m.Draw("same");
    if (m_fillMesh) m.Draw("f:sames");
  }

  // Draw the drift lines on the view
  for (auto& dline : m_driftLines) {
    dline.Draw("sames");
  }
  // Draw axes again so they are on top
  gPad->RedrawAxis("g");
}

// Removes duplicate points and line crossings by correctly ordering
//  the points in the provided vectors.
//
//  NOTE: This is a 2D version of the BUTFLD method in Garfield.  It
//   follows the same general algorithm.
//
// TODO: there is an algorithm which always sorts points correctly using cross
// product, see IntersectPlaneArea.
void ViewFEMesh::RemoveCrossings(std::vector<double>& x,
                                 std::vector<double>& y) {
  // Determine element dimensions
  double xmin = x[0], xmax = x[0];
  double ymin = y[0], ymax = y[0];
  for (int i = 1; i < (int)x.size(); i++) {
    if (x[i] < xmin) xmin = x[i];
    if (x[i] > xmax) xmax = x[i];
    if (y[i] < ymin) ymin = y[i];
    if (y[i] > ymax) ymax = y[i];
  }

  // First remove duplicate points
  double xtol = 1e-10 * std::abs(xmax - xmin);
  double ytol = 1e-10 * std::abs(ymax - ymin);
  for (int i = 0; i < (int)x.size(); i++) {
    for (int j = i + 1; j < (int)x.size(); j++) {
      if (std::abs(x[i] - x[j]) < xtol && std::abs(y[i] - y[j]) < ytol) {
        x.erase(x.begin() + j);
        y.erase(y.begin() + j);
        j--;
      }
    }
  }

  // No crossings with 3 points or less
  if (x.size() <= 3) return;

  // Save the polygon size so it is easily accessible
  int NN = x.size();

  // Keep track of the number of attempts
  int attempts = 0;

  // Exchange points until crossings are eliminated or we have attempted NN
  // times
  bool crossings = true;
  while (crossings && (attempts < NN)) {
    // Assume we are done after this attempt.
    crossings = false;

    for (int i = 1; i <= NN; i++) {
      for (int j = i + 2; j <= NN; j++) {
        // End the j-loop if we have surpassed N and wrapped around to i
        if ((j + 1) > NN && 1 + (j % NN) >= i) break;
        // Otherwise, detect crossings and attempt to eliminate them.
        // Determine if we have a crossing.
        double xc = 0., yc = 0.;
        if (!LinesCrossed(x[(i - 1) % NN], y[(i - 1) % NN], x[i % NN],
                          y[i % NN], x[(j - 1) % NN], y[(j - 1) % NN],
                          x[j % NN], y[j % NN], xc, yc)) {
          continue;
        }
        // Swap each point from i towards j with each corresponding point
        // from j towards i.
        for (int k = 1; k <= (j - i) / 2; k++) {
          double xs = x[(i + k - 1) % NN];
          double ys = y[(i + k - 1) % NN];
          x[(i + k - 1) % NN] = x[(j - k) % NN];
          y[(i + k - 1) % NN] = y[(j - k) % NN];
          x[(j - k) % NN] = xs;
          y[(j - k) % NN] = ys;

          // Force another attempt
          crossings = true;
        }
      }  // end loop over j
    }    // end loop over i

    // Increment the number of attempts
    attempts++;

  }  // end while(crossings)

  if (attempts > NN) {
    std::cerr << m_className << "::RemoveCrossings:\n    Warning: "
              << "Maximum attempts reached. Crossings not removed.\n";
  }
}

/// Return true if the specified point is in the view region.
bool ViewFEMesh::InView(const double x, const double y) const {
  // Set up the view vertices.
  int vN = m_viewRegionLines[0].GetN() - 1;
  std::vector<double> vx(vN), vy(vN);
  vx.assign(m_viewRegionLines[0].GetX(), m_viewRegionLines[0].GetX() + vN);
  vy.assign(m_viewRegionLines[0].GetY(), m_viewRegionLines[0].GetY() + vN);
  // Test whether this vertex is inside the view.
  bool edge = false;
  return IsInPolygon(x, y, vx, vy, edge);
}

//
// Determines whether the line connecting points (x1,y1) and (x2,y2)
// and the line connecting points (u1,v1) and (u2,v2) cross somewhere
// between the 4 points.  Sets the crossing point in (xc, yc).
//
// Ported from Garfield function CROSSD
//
bool ViewFEMesh::LinesCrossed(double x1, double y1, double x2, double y2,
                              double u1, double v1, double u2, double v2,
                              double& xc, double& yc) const {
  // Set the tolerances.
  double xtol = 1.0e-10 * std::max({std::abs(x1), std::abs(x2), std::abs(u1),
                                    std::abs(u2)});
  double ytol = 1.0e-10 * std::max({std::abs(y1), std::abs(y2), std::abs(v1),
                                    std::abs(v2)});
  if (xtol <= 0) xtol = 1.0e-10;
  if (ytol <= 0) ytol = 1.0e-10;

  // Compute the distances and determinant (dx,dy) x (du,dv).
  double dy = y2 - y1;
  double dv = v2 - v1;
  double dx = x1 - x2;
  double du = u1 - u2;
  double det = dy * du - dx * dv;

  // Check for crossing because one of the endpoints is on both lines.
  if (OnLine(x1, y1, x2, y2, u1, v1)) {
    xc = u1;
    yc = v1;
    return true;
  } else if (OnLine(x1, y1, x2, y2, u2, v2)) {
    xc = u2;
    yc = v2;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x1, y1)) {
    xc = x1;
    yc = y1;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x2, y2)) {
    xc = x2;
    yc = y2;
    return true;
  }
  // Check if the lines are parallel (zero determinant).
  if (std::abs(det) < xtol * ytol) return false;
  // No special case: compute point of intersection.

  // Solve crossing equations.
  xc = (du * (x1 * y2 - x2 * y1) - dx * (u1 * v2 - u2 * v1)) / det;
  yc = ((-1 * dv) * (x1 * y2 - x2 * y1) + dy * (u1 * v2 - u2 * v1)) / det;

  // Determine if this point is on both lines.
  if (OnLine(x1, y1, x2, y2, xc, yc) && OnLine(u1, v1, u2, v2, xc, yc))
    return true;

  // The lines do not cross if we have reached this point.
  return false;
}

// Determines whether the point (u,v) lies on the line connecting
// points (x1,y1) and (x2,y2).
//
// Ported from Garfield function ONLIND
//
bool ViewFEMesh::OnLine(double x1, double y1, double x2, double y2, double u,
                        double v) const {
  // Set the tolerances
  double xtol = 1.e-10 * std::max({std::abs(x1), std::abs(x2), std::abs(u)});
  double ytol = 1.e-10 * std::max({std::abs(y1), std::abs(y2), std::abs(v)});
  if (xtol <= 0) xtol = 1.0e-10;
  if (ytol <= 0) ytol = 1.0e-10;

  // To store the coordinates of the comparison point
  double xc = 0, yc = 0;

  // Check if (u,v) is the point (x1,y1) or (x2,y2)
  if ((std::abs(x1 - u) <= xtol && std::abs(y1 - v) <= ytol) ||
      (std::abs(x2 - u) <= xtol && std::abs(y2 - v) <= ytol)) {
    return true;
  }
  // Check if the line is actually a point
  if (std::abs(x1 - x2) <= xtol && std::abs(y1 - y2) <= ytol) {
    return false;
  }
  // Choose (x1,y1) as starting point if closer to (u,v)
  if (std::abs(u - x1) + std::abs(v - y1) <
      std::abs(u - x2) + std::abs(v - y2)) {
    // Compute the component of the line from (x1,y1) to (u,v)
    // along the line from (x1,y1) to (x2,y2)
    double dpar = ((u - x1) * (x2 - x1) + (v - y1) * (y2 - y1)) /
                  ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

    // Determine the point on the line to which to compare (u,v)
    if (dpar < 0.0) {
      xc = x1;
      yc = y1;
    } else if (dpar > 1.0) {
      xc = x2;
      yc = y2;
    } else {
      xc = x1 + dpar * (x2 - x1);
      yc = y1 + dpar * (y2 - y1);
    }
  } else { 
    // Choose (x2,y2) as starting point if closer to (u,v)
    // Compute the component of the line from (x2,y2) to (u,v)
    //  along the line from (x2,y2) to (x1,y1)
    double dpar = ((u - x2) * (x1 - x2) + (v - y2) * (y1 - y2)) /
                  ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

    // Determine the point on the line to which to compare (u,v)
    if (dpar < 0.0) {
      xc = x2;
      yc = y2;
    } else if (dpar > 1.0) {
      xc = x1;
      yc = y1;
    } else {
      xc = x2 + dpar * (x1 - x2);
      yc = y2 + dpar * (y1 - y2);
    }
  }

  // Compare the calculated point to (u,v)
  if (std::abs(u - xc) < xtol && std::abs(v - yc) < ytol) return true;

  return false;
}

// Ported from Garfield: determines the point of intersection, in planar
// coordinates, of a plane with the line connecting multiple points
// x1,y1,z1;x2,y2,z2: the world coordinates of the two points
// projMat;planeMat: the projection and plane matrices
// xMat: the resulting planar coordinates of the intersection point
bool ViewFEMesh::PlaneCut(double x1, double y1, double z1, double x2, double y2,
                          double z2, TMatrixD& xMat) {
  // Set up the matrix for cutting edges not in the plane
  TArrayD dataCut(9);
  TMatrixD cutMat(3, 3);
  dataCut[0] = m_proj[0][0];
  dataCut[1] = m_proj[1][0];
  dataCut[2] = x1 - x2;
  dataCut[3] = m_proj[0][1];
  dataCut[4] = m_proj[1][1];
  dataCut[5] = y1 - y2;
  dataCut[6] = m_proj[0][2];
  dataCut[7] = m_proj[1][2];
  dataCut[8] = z1 - z2;
  cutMat.SetMatrixArray(dataCut.GetArray());

  // Calculate the determinant of the cut matrix
  double cutDet =
      cutMat(0, 0) *
          (cutMat(1, 1) * cutMat(2, 2) - cutMat(1, 2) * cutMat(2, 1)) -
      cutMat(0, 1) *
          (cutMat(1, 0) * cutMat(2, 2) - cutMat(1, 2) * cutMat(2, 0)) +
      cutMat(0, 2) *
          (cutMat(1, 0) * cutMat(2, 1) - cutMat(1, 1) * cutMat(2, 0));

  // Do not proceed if the matrix is singular
  if (std::abs(cutDet) < 1e-20) return false;

  // Set up a coordinate vector (RHS of equation)
  TArrayD dataCoords(3);
  TMatrixD coordMat(3, 1);
  dataCoords[0] = x1 - m_dist * m_proj[2][0];
  dataCoords[1] = y1 - m_dist * m_proj[2][1];
  dataCoords[2] = z1 - m_dist * m_proj[2][2];
  coordMat.SetMatrixArray(dataCoords.GetArray());

  // Invert the cut matrix and multiply to get the solution
  cutMat.SetTol(1e-20);
  cutMat.Invert();
  // Do not proceed if the matrix is singular
  if (!cutMat.IsValid()) return false;
  xMat = cutMat * coordMat;

  // Return success if the plane point is between the two vertices
  if (xMat(2, 0) < 0 || xMat(2, 0) > 1) return false;
  return true;
}

// Calculates m_viewRegionLines and canvas dimensions based on projection plane
// and view area
bool ViewFEMesh::IntersectPlaneArea(void) {
  std::vector<TMatrixD> intersect_points;
  m_viewRegionLines.clear();
  // Loop over box edges
  for (int x0 = 0; x0 < 2; ++x0) {
    for (int y0 = 0; y0 < 2; ++y0) {
      for (int z0 = 0; z0 < 2; ++z0) {
        for (int x1 = x0; x1 < 2; ++x1) {
          for (int y1 = y0; y1 < 2; ++y1) {
            for (int z1 = z0; z1 < 2; ++z1) {
              if (x1 - x0 + y1 - y0 + z1 - z0 != 1) continue;
              double X0 = (x0 ? m_xMin : m_xMax);
              double Y0 = (y0 ? m_yMin : m_yMax);
              double Z0 = (z0 ? m_zMin : m_zMax);
              double X1 = (x1 ? m_xMin : m_xMax);
              double Y1 = (y1 ? m_yMin : m_yMax);
              double Z1 = (z1 ? m_zMin : m_zMax);
              TMatrixD xMat(3, 1);
              if (!PlaneCut(X0, Y0, Z0, X1, Y1, Z1, xMat)) continue;
              if (m_debug) {
                std::cout << m_className << "::IntersectPlaneArea:\n"
                          << "    Intersection of plane at (" << xMat(0, 0) 
                          << ", " << xMat(1, 0) << ", " << xMat(2, 0) 
                          << ") with edge\n    (" 
                          << X0 << ", " << Y0 << ", " << Z0 << ")-(" 
                          << X1 << ", " << Y1 << ", " << Z1 << ")\n";
              }
              // Do not add same points (the case when plane contains an edge)
              bool skip = false;
              for (auto& p : intersect_points) {
                const double dx = xMat(0, 0) - p(0, 0);
                const double dy = xMat(1, 0) - p(1, 0);
                if (std::sqrt(dx * dx + dy * dy) < 1e-10) {
                  skip = true;
                  break;
                }
              }
              if (!skip) intersect_points.push_back(xMat);
            }
          }
        }
      }
    }
  }
  if (intersect_points.size() < 3) {
    std::cerr << m_className << "::IntersectPlaneArea:\n";
    std::cerr << "    WARNING: Empty intersection of view plane with area.\n";
    return false;
  }
  TMatrixD offset = intersect_points[0];
  m_xPlaneMin = m_xPlaneMax = intersect_points[0](0, 0);
  m_yPlaneMin = m_yPlaneMax = intersect_points[0](1, 0);
  // Remove crossings in resulting polyline by sorting points rotation-wise.
  for (auto& p : intersect_points) p -= offset;
  std::sort(intersect_points.begin(), intersect_points.end(),
            [](const TMatrixD& a, const TMatrixD& b) -> bool {
              double cross_z = a(0, 0) * b(1, 0) - a(1, 0) * b(0, 0);
              return cross_z < 0;
            });
  TPolyLine poly;
  poly.SetLineColor(plottingEngine.GetRootColorLine2());
  poly.SetLineWidth(3);
  std::size_t pn = 0;
  for (auto& p : intersect_points) {
    p += offset;
    poly.SetPoint(pn, p(0, 0), p(1, 0));
    m_xPlaneMin = std::min(p(0, 0), m_xPlaneMin);
    m_yPlaneMin = std::min(p(1, 0), m_yPlaneMin);
    m_xPlaneMax = std::max(p(0, 0), m_xPlaneMax);
    m_yPlaneMax = std::max(p(1, 0), m_yPlaneMax);
    ++pn;
  }
  poly.SetPoint(pn, offset(0, 0), offset(1, 0));
  m_viewRegionLines.push_back(poly);
  return true;
}

// In x,y,z: vector coordinates
// Out x,y,z: vector parallel to the viewing plane (project[3][3])
bool ViewFEMesh::PlaneVector(double& x, double& y, double& z) const {
  double dist = x * m_proj[2][0] + y * m_proj[2][1] + z * m_proj[2][2];
  x = x - dist * m_proj[2][0];
  y = y - dist * m_proj[2][1];
  z = z - dist * m_proj[2][2];
  return true;
}

// Ported from Garfield: calculates the planar coordinates
// x,y,z: original world coordinates
// projMat: the projection matrix
// xMat: the resulting planar coordinates in single-column (vector) form
bool ViewFEMesh::PlaneCoords(double x, double y, double z,
                             const TMatrixD& projMat, TMatrixD& xMat) {
  // Set up the coordinate vector
  TArrayD dataCoords(3);
  TMatrixD coordMat(3, 1);
  dataCoords[0] = x;
  dataCoords[1] = y;
  dataCoords[2] = z;
  coordMat.SetMatrixArray(dataCoords.GetArray());
  xMat = projMat * coordMat;

  return true;
}

// Ported from Garfield (function INTERD):
// Returns true if the point (x,y) is inside of the specified polygon.
// x: the x-coordinate
// y: the y-coordinate
// px: the x-vertices of the polygon
// py: the y-vertices of the polygon
// edge: a variable set to true if the point is located on the polygon edge
bool ViewFEMesh::IsInPolygon(double x, double y, std::vector<double>& px,
                             std::vector<double>& py, bool& edge) const {
  // Get the number and coordinates of the polygon vertices.
  int pN = (int)px.size();

  // Handle the special case of less than 2 vertices.
  if (pN < 2) return false;
  // Handle the special case of exactly 2 vertices (a line).
  if (pN == 2) return OnLine(px[0], py[0], px[1], py[1], x, y);

  // Set the minimum and maximum coordinates of all polygon vertices.
  double px_min = px[0], py_min = py[0];
  double px_max = px[0], py_max = py[0];
  for (int i = 0; i < pN; i++) {
    px_min = std::min(px_min, px[i]);
    py_min = std::min(py_min, py[i]);
    px_max = std::max(px_max, px[i]);
    py_max = std::max(py_max, py[i]);
  }

  // Set the tolerances
  double xtol = 1.0e-10 * std::max(std::abs(px_min), std::abs(px_max));
  double ytol = 1.0e-10 * std::max(std::abs(py_min), std::abs(py_max));
  if (xtol <= 0) xtol = 1.0e-10;
  if (ytol <= 0) ytol = 1.0e-10;

  // If we have essentially one x value, check to see if y is in range.
  if (std::abs(px_max - px_min) < xtol) {
    edge = (y > (py_min - ytol) && y < (py_max + ytol) &&
            std::abs(px_max + px_min - 2 * x) < xtol);
    return false;
  }
  // If we have essentially one y value, check to see if x is in range.
  if (std::abs(py_max - py_min) < ytol) {
    edge = (x > (px_min - xtol) && x < (px_max + xtol) &&
            std::abs(py_max + py_min - 2 * y) < ytol);
    return false;
  }

  // Set "infinity" points.
  double xinf = px_min - std::abs(px_max - px_min);
  double yinf = py_min - std::abs(py_max - py_min);

  // Loop until successful or maximum iterations (100) reached.
  int niter = 0;
  bool done = false;
  int ncross = 0;

  while (!done && niter < 100) {
    // Assume we will finish on this loop.
    done = true;

    // Loop over all edges, counting the number of edges crossed by a line
    // extending from (x, y) to (xinf, yinf).
    ncross = 0;
    for (int i = 0; (done && i < pN); i++) {
      // Determine whether the point lies on the edge.
      if (OnLine(px[i % pN], py[i % pN], px[(i + 1) % pN], py[(i + 1) % pN], x,
                 y)) {
        edge = true;
        return false;
      }

      // Determine whether this edge is crossed; if so increment the counter.
      double xc = 0., yc = 0.;
      if (LinesCrossed(x, y, xinf, yinf, px[i % pN], py[i % pN],
                       px[(i + 1) % pN], py[(i + 1) % pN], xc, yc))
        ncross++;

      // Ensure this vertex is not crossed by the line from (x,y)
      //  to (xinf,yinf); if so recompute (xinf,yinf) and start over.
      if (OnLine(x, y, xinf, yinf, px[i], py[i])) {
        // Recompute (xinf,yinf).
        xinf = px_min - RndmUniform() * std::abs(px_max - xinf);
        yinf = py_min - RndmUniform() * std::abs(py_max - yinf);

        // Start over.
        done = false;
        niter++;
      }
    }
  }

  // If we failed to finish iterating, return false.
  if (niter >= 100) {
    std::cerr << m_className << "::IsInPolygon: Unable to determine whether ("
              << x << ", " << y << ") is inside a polygon. Returning false.\n";
    return false;
  }

  // Point is inside for an odd, nonzero number of crossings.
  return (ncross != 2 * (ncross / 2));
}

// Ported from Garfield (method GRCONV):
// Clip the specified polygon to the view region; return the clipped polygon.
// px: the x-vertices of the polygon
// py: the y-vertices of the polygon
// cx: to contain the x-vertices of the clipped polygon
// cy: to contain the y-vertices of the clipped polygon
void ViewFEMesh::ClipToView(std::vector<double>& px, std::vector<double>& py,
                            std::vector<double>& cx, std::vector<double>& cy) {
  // Get the number and coordinates of the polygon vertices.
  int pN = (int)px.size();

  // Clear the vectors to contain the final polygon.
  cx.clear();
  cy.clear();

  // Set up the view vertices.
  int vN = m_viewRegionLines[0].GetN() - 1;
  std::vector<double> vx(vN), vy(vN);
  vx.assign(m_viewRegionLines[0].GetX(), m_viewRegionLines[0].GetX() + vN);
  vy.assign(m_viewRegionLines[0].GetY(), m_viewRegionLines[0].GetY() + vN);

  // Do nothing if we have less than 2 points.
  if (pN < 2) return;

  // Loop over the polygon vertices.
  for (int i = 0; i < pN; i++) {
    // Flag for skipping check for edge intersection.
    bool skip = false;

    // Loop over the view vertices.
    for (int j = 0; j < vN; j++) {
      // Determine whether this vertex lies on a view edge:
      // if so add the vertex to the final polygon.
      if (OnLine(vx[j % vN], vy[j % vN], vx[(j + 1) % vN], vy[(j + 1) % vN],
                 px[i], py[i])) {
        // Add the vertex.
        cx.push_back(px[i]);
        cy.push_back(py[i]);

        // Skip edge intersection check in this case.
        skip = true;
      }

      // Determine whether a corner of the view area lies on this edge:
      // if so add the corner to the final polygon.
      if (OnLine(px[i % pN], py[i % pN], px[(i + 1) % pN], py[(i + 1) % pN],
                 vx[j], vy[j])) {
        // Add the vertex.
        cx.push_back(vx[j]);
        cy.push_back(vy[j]);

        // Skip edge intersection check in this case.
        skip = true;
      }
    }

    // If we have not skipped the edge intersection check, look for an
    // intersection between this edge and the view edges.
    if (skip) continue;
    // Loop over the view vertices.
    for (int j = 0; j < vN; j++) {
      // Check for a crossing with this edge;
      // if one exists, add the crossing point.
      double xc = 0., yc = 0.;
      if (LinesCrossed(vx[j % vN], vy[j % vN], vx[(j + 1) % vN],
                       vy[(j + 1) % vN], px[i % pN], py[i % pN],
                       px[(i + 1) % pN], py[(i + 1) % pN], xc, yc)) {
        // Add a vertex.
        cx.push_back(xc);
        cy.push_back(yc);
      }
    }
  }

  // Find all view field vertices inside the polygon.
  for (int j = 0; j < vN; j++) {
    // Test whether this vertex is inside the polygon.
    // If so, add it to the final polygon.
    bool edge = false;
    if (IsInPolygon(vx[j], vy[j], px, py, edge)) {
      // Add the view vertex.
      cx.push_back(vx[j]);
      cy.push_back(vy[j]);
    }
  }

  // Find all polygon vertices inside the box.
  for (int i = 0; i < pN; i++) {
    // Test whether this vertex is inside the view.
    // If so, add it to the final polygon.
    bool edge = false;
    if (IsInPolygon(px[i], py[i], vx, vy, edge)) {
      // Add the polygon vertex.
      cx.push_back(px[i]);
      cy.push_back(py[i]);
    }
  }
}
}  // namespace Garfield
