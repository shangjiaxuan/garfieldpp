#include <cmath>
#include <iostream>

#include <TGeoBBox.h>
#include <TGeoCone.h>
#include <TGeoBoolNode.h>
#include <TGeoCompositeShape.h>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Solid.hh"
#include "Garfield/ViewGeometry.hh"

namespace Garfield {

ViewGeometry::ViewGeometry() : ViewBase("ViewGeometry") { 
  plottingEngine.SetDefaultStyle(); 
}

ViewGeometry::~ViewGeometry() {
  Reset();
}

void ViewGeometry::SetGeometry(GeometrySimple* geo) {
  if (!geo) {
    std::cerr << m_className << "::SetGeometry: Null pointer.\n";
    return;
  }

  m_geometry = geo;
}

void ViewGeometry::Plot() {
  if (!m_geometry) {
    std::cerr << m_className << "::Plot: Geometry is not defined.\n";
    return;
  }

  if (!m_canvas) {
    m_canvas = new TCanvas();
    if (m_hasExternalCanvas) m_hasExternalCanvas = false;
  }
  m_canvas->cd();

  const unsigned int nSolids = m_geometry->GetNumberOfSolids();
  if (nSolids == 0) {
    std::cerr << m_className << "::Plot: Geometry is empty.\n";
    return;
  }

  // Get the bounding box.
  double xMin = 0., yMin = 0., zMin = 0.;
  double xMax = 0., yMax = 0., zMax = 0.;
  if (!m_geometry->GetBoundingBox(xMin, yMin, zMin, xMax, yMax, zMax)) {
    std::cerr << m_className << "::Plot: Cannot retrieve bounding box.\n";
    return;
  }
  m_geoManager.reset(new TGeoManager("ViewGeometryGeoManager", ""));
  TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0., 0., 0.);
  TGeoMedium* medVacuum = new TGeoMedium("Vacuum", 1, matVacuum);
  m_media.push_back(medVacuum);
  // Use silicon as "default" material.
  TGeoMaterial* matDefault = new TGeoMaterial("Default", 28.085, 14., 2.329);
  TGeoMedium* medDefault = new TGeoMedium("Default", 1, matDefault);
  TGeoVolume* world = m_geoManager->MakeBox(
      "World", medVacuum, std::max(fabs(xMin), fabs(xMax)),
      std::max(fabs(yMin), fabs(yMax)), std::max(fabs(zMin), fabs(zMax)));
  m_geoManager->SetTopVolume(world);
  m_volumes.push_back(world);

  for (unsigned int i = 0; i < nSolids; ++i) {
    Solid* solid = m_geometry->GetSolid(i);
    if (!solid) {
      std::cerr << m_className << "::Plot:\n"
                << "    Could not get solid " << i << " from geometry.\n";
      continue;
    }
    // Get the center coordinates.
    double x0 = 0., y0 = 0., z0 = 0.;
    if (!solid->GetCentre(x0, y0, z0)) {
      std::cerr << m_className << "::Plot: Could not determine solid centre.\n";
      continue;
    }
    // Get the rotation.
    double ctheta = 1., stheta = 0.;
    double cphi = 1., sphi = 0.;
    if (!solid->GetOrientation(ctheta, stheta, cphi, sphi)) {
      std::cerr << m_className << "::Plot:\n"
                << "    Could not determine solid orientation.\n";
      continue;
    }
    double matrix[9] = {cphi * ctheta, -sphi, cphi * stheta,
                        sphi * ctheta, cphi,  sphi * stheta,
                        -stheta,       0,     ctheta};
    TGeoVolume* volume = nullptr;
    if (solid->IsTube()) {
      const double rmin = solid->GetInnerRadius();
      const double rmax = solid->GetOuterRadius();
      const double lz = solid->GetHalfLengthZ();
      volume = m_geoManager->MakeTube("Tube", medDefault, rmin, rmax, lz);
    } else if (solid->IsBox()) {
      const double dx = solid->GetHalfLengthX();
      const double dy = solid->GetHalfLengthY();
      const double dz = solid->GetHalfLengthZ();
      volume = m_geoManager->MakeBox("Box", medDefault, dx, dy, dz);
    } else if (solid->IsSphere()) {
      const double rmin = solid->GetInnerRadius();
      const double rmax = solid->GetOuterRadius();
      volume = m_geoManager->MakeSphere("Sphere", medDefault, rmin, rmax);
    } else if (solid->IsHole()) {
      const double r1 = solid->GetLowerRadius();
      const double r2 = solid->GetUpperRadius();
      const double rm = 0.5 * (r1 + r2);
      const double dr = (r2 - r1);
      const double dx = solid->GetHalfLengthX();
      const double dy = solid->GetHalfLengthY();
      const double dz = solid->GetHalfLengthZ();
      TGeoBBox* box = new TGeoBBox("HoleBox", dx, dy, dz);
      TGeoCone* cone = new TGeoCone("HoleCone", 2 * dz, 0, rm - dr, 0, rm + dr);
      TGeoCompositeShape* hole = new TGeoCompositeShape("Hole", 
        new TGeoSubtraction(box, cone));
      hole->RegisterYourself();
      volume = new TGeoVolume("Hole", hole, medDefault); 
    } else {
      std::cerr << m_className << "::Plot: Unknown type of solid.\n";
      continue;
    }
    Medium* medium = m_geometry->GetMedium(x0, y0, z0);
    if (!medium) {
      volume->SetLineColor(kGreen + 2);
      volume->SetTransparency(50);
    } else if (medium->IsGas()) {
      volume->SetLineColor(kBlue + medium->GetId());
      volume->SetTransparency(50);
    } else if (medium->IsSemiconductor()) {
      volume->SetLineColor(kRed + medium->GetId());
      volume->SetTransparency(50);
    } else {
      volume->SetLineColor(kViolet + medium->GetId());
      volume->SetTransparency(0);
    }
    TGeoRotation r;
    r.SetMatrix(matrix);
    TGeoTranslation t(x0, y0, z0);
    TGeoCombiTrans* transform = new TGeoCombiTrans(t, r);
    m_volumes.push_back(volume);
    m_geoManager->GetTopVolume()->AddNode(volume, 1, transform);
  }
  m_geoManager->CloseGeometry();
  m_geoManager->GetTopNode()->Draw("ogl");
}

void ViewGeometry::Reset() {
  for (auto it = m_volumes.begin(), end = m_volumes.end(); it != end; ++it) {
    if (*it) {
      TGeoShape* shape = (*it)->GetShape();
      if (shape) delete shape;
      delete *it;
    }
  }
  m_volumes.clear();
  for (auto it = m_media.begin(), end = m_media.end(); it != end; ++it) {
    if (*it) {
      TGeoMaterial* material = (*it)->GetMaterial();
      if (material) delete material;
      delete *it;
    }
  }
  m_media.clear();

  m_geoManager.reset(nullptr);
}
}
