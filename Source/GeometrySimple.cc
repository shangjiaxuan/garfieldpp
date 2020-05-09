#include <algorithm>
#include <iostream>

#include "Garfield/GeometrySimple.hh"

namespace Garfield {

GeometrySimple::GeometrySimple() : Geometry("GeometrySimple") {}

void GeometrySimple::AddSolid(Solid* solid, Medium* medium) {
  // Make sure the solid and the medium are defined.
  if (!solid || !medium) {
    std::cerr << m_className << "::AddSolid: Null pointer.\n";
    return;
  }

  // Update the bounding box ranges
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  if (!solid->GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)) {
    std::cerr << m_className << "::AddSolid: Solid has no bounding box.\n";
    return;
  }

  if (m_hasBoundingBox) {
    m_xMinBoundingBox = std::min(m_xMinBoundingBox, xmin);
    m_yMinBoundingBox = std::min(m_yMinBoundingBox, ymin);
    m_zMinBoundingBox = std::min(m_zMinBoundingBox, zmin);
    m_xMaxBoundingBox = std::max(m_xMaxBoundingBox, xmax);
    m_yMaxBoundingBox = std::max(m_yMaxBoundingBox, ymax);
    m_zMaxBoundingBox = std::max(m_zMaxBoundingBox, zmax);
  } else {
    m_xMinBoundingBox = xmin;
    m_yMinBoundingBox = ymin;
    m_zMinBoundingBox = zmin;
    m_xMaxBoundingBox = xmax;
    m_yMaxBoundingBox = ymax;
    m_zMaxBoundingBox = zmax;
    m_hasBoundingBox = true;
  }

  // Add the new solid to the list.
  m_solids.emplace_back(std::make_pair(solid, medium));
}

Solid* GeometrySimple::GetSolid(const double x, const double y,
                                const double z) const {
  for (const auto& solid : m_solids) {
    if (solid.first->IsInside(x, y, z)) return solid.first;
  }
  return nullptr;
}

Medium* GeometrySimple::GetMedium(const double x, const double y,
                                  const double z) const {
  for (const auto& solid : m_solids) {
    if (solid.first->IsInside(x, y, z)) {
      return solid.second;
    }
  }
  return m_medium;
}

Solid* GeometrySimple::GetSolid(const unsigned int i) const {
  if (i >= m_solids.size()) {
    std::cerr << m_className << "::GetSolid:\n"
              << "    Requested solid " << i << " does not exist.\n";
    return nullptr;
  }
  return m_solids[i].first;
}

Solid* GeometrySimple::GetSolid(const unsigned int i, Medium*& medium) const {
  if (i >= m_solids.size()) {
    std::cerr << m_className << "::GetSolid:\n"
              << "    Requested solid " << i << " does not exist.\n";
    return nullptr;
  }
  medium = m_solids[i].second;
  return m_solids[i].first;
}

void GeometrySimple::Clear() {
  m_solids.clear();
  m_medium = nullptr;
}

void GeometrySimple::PrintSolids() {
  std::cout << m_className << "::PrintSolids:\n";
  const unsigned int nSolids = m_solids.size();
  if (nSolids == 1) {
    std::cout << "    1 solid\n";
  } else {
    std::cout << "    " << nSolids << " solids\n";
  }
  if (m_solids.empty()) return;
  std::cout << "      Index      Type    Medium\n";
  for (unsigned int i = 0; i < nSolids; ++i) {
    std::cout << "        " << i << "         ";
    if (m_solids[i].first->IsBox()) {
      std::cout << "box      ";
    } else if (m_solids[i].first->IsTube()) {
      std::cout << "tube     ";
    } else if (m_solids[i].first->IsSphere()) {
      std::cout << "sphere   ";
    } else if (m_solids[i].first->IsHole()) {
      std::cout << "hole     ";
    } else if (m_solids[i].first->IsRidge()) {
      std::cout << "ridge    ";
    } else {
      std::cout << "unknown  ";
    }
    if (m_solids[i].second) {
      std::cout << m_solids[i].second->GetName() << "\n";
    } else {
      std::cout << " ---\n";
    }
  }
}

bool GeometrySimple::IsInside(const double x, const double y,
                              const double z) const {
  if (!IsInBoundingBox(x, y, z)) return false;

  for (const auto& solid : m_solids) {
    if (solid.first->IsInside(x, y, z)) return true;
  }
  return false;
}

bool GeometrySimple::IsInBoundingBox(const double x, const double y,
                                     const double z) const {
  if (!m_hasBoundingBox) {
    if (m_debug) {
      std::cerr << m_className << "::IsInBoundingBox:\n"
                << "    Bounding box is not defined.\n";
    }
    return true;
  }

  if (x >= m_xMinBoundingBox && x <= m_xMaxBoundingBox &&
      y >= m_yMinBoundingBox && y <= m_yMaxBoundingBox &&
      z >= m_zMinBoundingBox && z <= m_zMaxBoundingBox)
    return true;
  return false;
}
}
