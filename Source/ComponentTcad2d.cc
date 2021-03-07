#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Utilities.hh"

namespace Garfield {

ComponentTcad2d::ComponentTcad2d() : ComponentTcadBase("Tcad2d") {}

void ComponentTcad2d::ElectricField(const double xin, const double yin,
                                    const double zin, double& ex, double& ey,
                                    double& ez, double& p, Medium*& m,
                                    int& status) {
  // Assume this will work.
  status = 0;
  // Initialise.
  ex = ey = ez = p = 0.;
  m = nullptr;

  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::ElectricField:\n"
              << "    Field map is not available for interpolation.\n";
    status = -10;
    return;
  }

  if (m_hasRangeZ && (zin < m_bbMin[2] || zin > m_bbMax[2])) {
    status = -6;
    return;
  }
  // In case of periodicity, reduce to the cell volume.
  std::array<double, 2> x = {xin, yin};
  std::array<bool, 2> mirr = {false, false};
  MapCoordinates(x, mirr);
  // Check if the point is inside the bounding box.
  if (!InBoundingBox(x)) {
    status = -6;
    return;
  }

  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x[0], x[1], w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    status = -6;
    return;
  }

  const Element& element = m_elements[i];
  const size_t nVertices = element.type + 1;
  for (size_t j = 0; j < nVertices; ++j) {
    const size_t index = element.vertex[j];
    ex += w[j] * m_efield[index][0];
    ey += w[j] * m_efield[index][1];
    p += w[j] * m_potential[index];
  }
  if (mirr[0]) ex = -ex;
  if (mirr[1]) ey = -ey;
  m = m_regions[element.region].medium;
  if (!m_regions[element.region].drift || !m) status = -5;
  m_lastElement = i;
}

bool ComponentTcad2d::Interpolate(const double xin, const double yin,
    const double z,
    const std::vector<std::array<double, 2> >& field,
    double& fx, double& fy, double& fz) {

  fx = fy = fz = 0.;
  if (field.empty()) return false;
  if (m_hasRangeZ && (z < m_bbMin[2] || z > m_bbMax[2])) return false;
  std::array<double, 2> x = {xin, yin};
  std::array<bool, 2> mirr = {false, false};
  // In case of periodicity, reduce to the cell volume.
  MapCoordinates(x, mirr);
  // Make sure the point is inside the bounding box.
  if (!InBoundingBox(x)) return false;

  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x[0], x[1], w);
  // Stop if the point is outside the mesh.
  if (i >= m_elements.size()) return false;

  const Element& element = m_elements[i];
  const size_t nVertices = element.type + 1;
  for (size_t j = 0; j < nVertices; ++j) {
    const auto index = element.vertex[j];
    fx += w[j] * field[index][0];
    fy += w[j] * field[index][1];
  }
  if (mirr[0]) fx = -fx;
  if (mirr[1]) fy = -fy;
  m_lastElement = i;
  return true;
}

bool ComponentTcad2d::Interpolate(const double xin, const double yin,
    const double z,
    const std::vector<double>& field, double& f) {

  f = 0.;
  if (field.empty()) return false;
  if (m_hasRangeZ && (z < m_bbMin[2] || z > m_bbMax[2])) return false;
  std::array<double, 2> x = {xin, yin};
  std::array<bool, 2> mirr = {false, false};
  // In case of periodicity, reduce to the cell volume.
  MapCoordinates(x, mirr);
  if (!InBoundingBox(x)) return false;

  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x[0], x[1], w);
  // Stop if the point is outside the mesh.
  if (i >= m_elements.size()) return false;

  const Element& element = m_elements[i];
  const size_t nVertices = element.type + 1;
  for (size_t j = 0; j < nVertices; ++j) {
    f += w[j] * field[element.vertex[j]];
  }
  m_lastElement = i;
  return true;
}

Medium* ComponentTcad2d::GetMedium(const double xin, const double yin,
                                   const double zin) {
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::GetMedium:\n"
              << "    Field map not available for interpolation.\n";
    return nullptr;
  }

  if (m_hasRangeZ && (zin < m_bbMin[2] || zin > m_bbMax[2])) return nullptr;
  std::array<double, 2> x = {xin, yin};
  std::array<bool, 2> mirr = {false, false};
  MapCoordinates(x, mirr);
  // Check if the point is inside the bounding box.
  if (!InBoundingBox(x)) return nullptr;

  // Determine the shape functions.
  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x[0], x[1], w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    return nullptr;
  }
  m_lastElement = i;
  return m_regions[m_elements[i].region].medium;
}

bool ComponentTcad2d::Initialise(const std::string& gridfilename,
                                 const std::string& datafilename) {
  m_ready = false;
  Cleanup();

  // Import mesh data from .grd file.
  if (!LoadGrid(gridfilename)) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Importing mesh data failed.\n";
    Cleanup();
    return false;
  }

  // Import electric field, potential and other data from .dat file.
  if (!LoadData(datafilename)) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Importing electric field and potential failed.\n";
    Cleanup();
    return false;
  }

  // Find min./max. coordinates and potentials.
  m_bbMax[0] = m_vertices[m_elements[0].vertex[0]][0];
  m_bbMax[1] = m_vertices[m_elements[0].vertex[0]][1];
  m_bbMin = m_bbMax;
  for (auto& element : m_elements) {
    const auto& v0 = m_vertices[element.vertex[0]];
    double xmin = v0[0];
    double xmax = v0[0];
    double ymin = v0[1];
    double ymax = v0[1];
    const size_t nVertices = element.type + 1;
    for (size_t j = 0; j < nVertices; ++j) {
      const auto& v = m_vertices[element.vertex[j]];
      xmin = std::min(xmin, v[0]);
      xmax = std::max(xmax, v[0]);
      ymin = std::min(ymin, v[1]);
      ymax = std::max(ymax, v[1]);
    }
    constexpr double tol = 1.e-6;
    element.bbMin[0] = xmin - tol;
    element.bbMax[0] = xmax + tol;
    element.bbMin[1] = ymin - tol;
    element.bbMax[1] = ymax + tol;
    m_bbMin[0] = std::min(m_bbMin[0], xmin);
    m_bbMax[0] = std::max(m_bbMax[0], xmax);
    m_bbMin[1] = std::min(m_bbMin[1], ymin);
    m_bbMax[1] = std::max(m_bbMax[1], ymax);
  }
  m_pMin = *std::min_element(m_potential.begin(), m_potential.end());
  m_pMax = *std::max_element(m_potential.begin(), m_potential.end());
  std::cout << m_className << "::Initialise:\n"
            << "    Available data:\n";
  if (!m_potential.empty()) std::cout << "      Electrostatic potential\n";
  if (!m_efield.empty()) std::cout << "      Electric field\n";
  if (!m_eMobility.empty()) std::cout << "      Electron mobility\n";
  if (!m_hMobility.empty()) std::cout << "      Hole mobility\n";
  if (!m_eVelocity.empty()) std::cout << "      Electron velocity\n";
  if (!m_hVelocity.empty()) std::cout << "      Hole velocity\n";
  if (!m_eLifetime.empty()) std::cout << "      Electron lifetime\n";
  if (!m_hLifetime.empty()) std::cout << "      Hole lifetime\n";
  if (!m_donors.empty()) {
    std::cout << "      " << m_donors.size() << " donor-type traps\n";
  }
  if (!m_acceptors.empty()) {
    std::cout << "      " << m_acceptors.size() << " acceptor-type traps\n";
  }
  std::cout << "    Bounding box:\n"
            << "      " << m_bbMin[0] << " < x [cm] < " << m_bbMax[0] << "\n"
            << "      " << m_bbMin[1] << " < y [cm] < " << m_bbMax[1] << "\n"
            << "    Voltage range:\n"
            << "      " << m_pMin << " < V < " << m_pMax << "\n";

  bool ok = true;

  // Count the number of elements belonging to a region.
  const auto nRegions = m_regions.size();
  std::vector<int> nElementsRegion(nRegions, 0);

  // Count the different element shapes.
  unsigned int nPoints = 0;
  unsigned int nLines = 0;
  unsigned int nTriangles = 0;
  unsigned int nRectangles = 0;
  unsigned int nOtherShapes = 0;

  // Check if there are elements which are not part of any region.
  unsigned int nLoose = 0;
  std::vector<int> looseElements;

  // Check if there are degenerate elements.
  unsigned int nDegenerate = 0;
  std::vector<int> degenerateElements;

  const auto nElements = m_elements.size();
  for (size_t i = 0; i < nElements; ++i) {
    const Element& element = m_elements[i];
    if (element.type == 0) {
      ++nPoints;
    } else if (element.type == 1) {
      ++nLines;
      if (element.vertex[0] == element.vertex[1]) {
        degenerateElements.push_back(i);
        ++nDegenerate;
      }
    } else if (element.type == 2) {
      ++nTriangles;
      if (element.vertex[0] == element.vertex[1] ||
          element.vertex[1] == element.vertex[2] ||
          element.vertex[2] == element.vertex[0]) {
        degenerateElements.push_back(i);
        ++nDegenerate;
      }
    } else if (element.type == 3) {
      ++nRectangles;
      if (element.vertex[0] == element.vertex[1] ||
          element.vertex[0] == element.vertex[2] ||
          element.vertex[0] == element.vertex[3] ||
          element.vertex[1] == element.vertex[2] ||
          element.vertex[1] == element.vertex[3] ||
          element.vertex[2] == element.vertex[3]) {
        degenerateElements.push_back(i);
        ++nDegenerate;
      }
    } else {
      // Other shapes should not occur, since they were excluded in LoadGrid.
      ++nOtherShapes;
    }
    if (element.region < nRegions) {
      ++nElementsRegion[element.region];
    } else {
      looseElements.push_back(i);
      ++nLoose;
    }
  }

  if (nDegenerate > 0) {
    std::cerr << m_className << "::Initialise:\n"
              << "    The following elements are degenerate:\n";
    for (unsigned int i = 0; i < nDegenerate; ++i) {
      std::cerr << "      " << degenerateElements[i] << "\n";
    }
    ok = false;
  }

  if (nLoose > 0) {
    std::cerr << m_className << "::Initialise:\n"
              << "    The following elements are not part of any region:\n";
    for (unsigned int i = 0; i < nLoose; ++i) {
      std::cerr << "      " << looseElements[i] << "\n";
    }
    ok = false;
  }

  std::cout << m_className << "::Initialise:\n"
            << "    Number of regions: " << nRegions << "\n";
  for (size_t i = 0; i < nRegions; ++i) {
    std::cout << "      " << i << ": " << m_regions[i].name << ", "
              << nElementsRegion[i] << " elements\n";
  }

  std::cout << "    Number of elements: " << nElements << "\n";
  if (nPoints > 0) {
    std::cout << "      " << nPoints << " points\n";
  }
  if (nLines > 0) {
    std::cout << "      " << nLines << " lines\n";
  }
  if (nTriangles > 0) {
    std::cout << "      " << nTriangles << " triangles\n";
  }
  if (nRectangles > 0) {
    std::cout << "      " << nRectangles << " rectangles\n";
  }
  if (nOtherShapes > 0) {
    std::cerr << "      " << nOtherShapes << " elements of unknown type.\n"
              << "      Program bug!\n";
    m_ready = false;
    Cleanup();
    return false;
  }

  std::cout << "    Number of vertices: " << m_vertices.size() << "\n";
  if (!ok) {
    m_ready = false;
    Cleanup();
    return false;
  }

  // Set up the quad tree.
  const double hx = 0.5 * (m_bbMax[0] - m_bbMin[0]);
  const double hy = 0.5 * (m_bbMax[1] - m_bbMin[1]);
  m_tree.reset(new QuadTree(m_bbMin[0] + hx, m_bbMin[1] + hy, hx, hy));
  // Insert the mesh nodes in the tree.
  const auto nVertices = m_vertices.size();
  for (size_t i = 0; i < nVertices; ++i) {
    m_tree->InsertMeshNode(m_vertices[i][0], m_vertices[i][1], i);
  }

  // Insert the mesh elements in the tree.
  for (size_t i = 0; i < nElements; ++i) {
    const Element& e = m_elements[i];
    const double bb[4] = {e.bbMin[0], e.bbMin[1], e.bbMax[0], e.bbMax[1]};
    m_tree->InsertMeshElement(bb, i);
  }

  m_ready = true;
  UpdatePeriodicity();
  std::cout << m_className << "::Initialise: Initialisation finished.\n";
  return true;
}

bool ComponentTcad2d::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                     double& xmax, double& ymax, double& zmax) {
  if (!m_ready) return false;
  if (m_periodic[0] || m_mirrorPeriodic[0]) {
    xmin = -std::numeric_limits<double>::infinity();
    xmax = +std::numeric_limits<double>::infinity();
  } else {
    xmin = m_bbMin[0];
    xmax = m_bbMax[0];
  }

  if (m_periodic[1] || m_mirrorPeriodic[1]) {
    ymin = -std::numeric_limits<double>::infinity();
    ymax = +std::numeric_limits<double>::infinity();
  } else {
    ymin = m_bbMin[1];
    ymax = m_bbMax[1];
  }

  if (m_hasRangeZ) {
    zmin = m_bbMin[2];
    zmax = m_bbMax[2];
  }
  return true;
}

void ComponentTcad2d::SetRangeZ(const double zmin, const double zmax) {
  if (fabs(zmax - zmin) <= 0.) {
    std::cerr << m_className << "::SetRangeZ: Zero range is not permitted.\n";
    return;
  }
  m_bbMin[2] = std::min(zmin, zmax);
  m_bbMax[2] = std::max(zmin, zmax);
  m_hasRangeZ = true;
}

bool ComponentTcad2d::GetElement(const size_t i, double& vol,
                                 double& dmin, double& dmax, int& type,
                                 std::vector<size_t>& nodes, int& reg) const {
  nodes.clear();
  if (i >= m_elements.size()) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }

  const Element& element = m_elements[i];
  if (element.type == 0) {
    dmin = dmax = vol = 0;
  } else if (element.type == 1) {
    const auto& v0 = m_vertices[element.vertex[0]];
    const auto& v1 = m_vertices[element.vertex[1]];
    const double d = std::hypot(v1[0] - v0[0], v1[1] - v0[1]);
    dmin = dmax = vol = d;
  } else if (m_elements[i].type == 2) {
    const auto& v0 = m_vertices[element.vertex[0]];
    const auto& v1 = m_vertices[element.vertex[1]];
    const auto& v2 = m_vertices[element.vertex[2]];
    vol = 0.5 * fabs((v2[0] - v0[0]) * (v1[1] - v0[1]) - 
                     (v2[1] - v0[1]) * (v1[0] - v0[0]));
    const double a = std::hypot(v1[0] - v0[0], v1[1] - v0[1]);
    const double b = std::hypot(v2[0] - v0[0], v2[1] - v0[1]);
    const double c = std::hypot(v1[0] - v2[0], v1[1] - v2[1]);
    dmin = std::min({a, b, c});
    dmax = std::max({a, b, c});
  } else if (m_elements[i].type == 3) {
    const auto& v0 = m_vertices[element.vertex[0]];
    const auto& v1 = m_vertices[element.vertex[1]];
    const auto& v3 = m_vertices[element.vertex[3]];
    const double a = std::hypot(v1[0] - v0[0], v1[1] - v0[1]);
    const double b = std::hypot(v3[0] - v0[0], v3[1] - v0[1]);
    vol = a * b;
    dmin = std::min(a, b);
    dmax = sqrt(a * a + b * b);
  } else {
    std::cerr << m_className << "::GetElement:\n"
              << "    Unexpected element type (" << type << ")\n";
    return false;
  }
  const size_t nVertices = element.type + 1;
  for (size_t i = 0; i < nVertices; ++i) {
    nodes.push_back(element.vertex[0]);
  }
  reg = element.region;
  return true;
}

bool ComponentTcad2d::GetNode(const size_t i, double& x, double& y,
                              double& v, double& ex, double& ey) const {
  if (i >= m_vertices.size()) {
    std::cerr << m_className << "::GetNode: Index out of range.\n";
    return false;
  }

  x = m_vertices[i][0];
  y = m_vertices[i][1];
  if (!m_potential.empty()) v = m_potential[i];
  if (!m_efield.empty()) {
    ex = m_efield[i][0];
    ey = m_efield[i][1];
  }
  return true;
}

size_t ComponentTcad2d::FindElement(const double x, const double y, 
    std::array<double, nMaxVertices>& w) const {

  w.fill(0.);
 
  if (m_lastElement < m_elements.size()) {
    // Check if the point is still located in the previously found element.
    const Element& last = m_elements[m_lastElement];
    if (InElement(x, y, last, w)) return m_lastElement;
  }

  // The point is not in the previous element.
  std::vector<int> elementsToSearch;
  if (m_tree) elementsToSearch = m_tree->GetElementsInBlock(x, y);
  const size_t nElementsToSearch = m_tree ? elementsToSearch.size() : m_elements.size(); 
  for (size_t i = 0; i < nElementsToSearch; ++i) {
    const size_t idx = m_tree ? elementsToSearch[i] : i;
    if (InElement(x, y, m_elements[idx], w)) return idx;
  }
  // Point is outside the mesh.
  if (m_debug) {
    std::cerr << m_className << "::FindElement:\n"
              << "    Point (" << x << ", " << y << ") is outside the mesh.\n";
  }
  return m_elements.size();
}

bool ComponentTcad2d::InElement(const double x, const double y,
                                const Element& element,
                                std::array<double, nMaxVertices>& w) const {
  if (x < element.bbMin[0] || x > element.bbMax[0] || 
      y < element.bbMin[1] || y > element.bbMax[1]) {
    return false;
  }
  switch (element.type) {
    case 0:
      return AtPoint(x, y, element, w);
    case 1:
      return OnLine(x, y, element, w);
    case 2:
      return InTriangle(x, y, element, w);
    case 3:
      return InRectangle(x, y, element, w);
    default:
      std::cerr << m_className << "::InElement:\n"
                << "    Unknown element type. Program bug!\n";
      break;
  }
  return false;
}

bool ComponentTcad2d::InRectangle(const double x, const double y,
                                  const Element& element,
                                  std::array<double, nMaxVertices>& w) const {
  const auto& v0 = m_vertices[element.vertex[0]];
  const auto& v1 = m_vertices[element.vertex[1]];
  const auto& v3 = m_vertices[element.vertex[3]];
  if (y < v0[1] || x > v3[0] || y > v1[1]) return false;

  // Map (x, y) to local variables (u, v) in [-1, 1].
  const double u = (x - 0.5 * (v0[0] + v3[0])) / (v3[0] - v0[0]);
  const double v = (y - 0.5 * (v0[1] + v1[1])) / (v1[1] - v0[1]);
  // Compute weighting factors for each corner.
  w[0] = (0.5 - u) * (0.5 - v);
  w[1] = (0.5 - u) * (0.5 + v);
  w[2] = (0.5 + u) * (0.5 + v);
  w[3] = (0.5 + u) * (0.5 - v);
  return true;
}

bool ComponentTcad2d::InTriangle(const double x, const double y,
                                 const Element& element,
                                 std::array<double, nMaxVertices>& w) const {
  const auto& v0 = m_vertices[element.vertex[0]];
  const auto& v1 = m_vertices[element.vertex[1]];
  const auto& v2 = m_vertices[element.vertex[2]];
  if (x > v1[0] && x > v2[0]) return false;
  if (y < v0[1] && y < v1[1] && y < v2[1]) return false;
  if (y > v0[1] && y > v1[1] && y > v2[1]) return false;

  // Map (x, y) onto local variables (b, c) such that
  // P = A + b * (B - A) + c * (C - A)
  // A point P is inside the triangle ABC if b, c > 0 and b + c < 1;
  // b, c are also weighting factors for points B, C
  const double sx = v1[0] - v0[0];
  const double sy = v1[1] - v0[1];
  const double tx = v2[0] - v0[0];
  const double ty = v2[1] - v0[1];
  const double d = 1. / (sx * ty - sy * tx);
  w[1] = ((x - v0[0]) * ty - (y - v0[1]) * tx) * d;
  if (w[1] < 0. || w[1] > 1.) return false;
  w[2] = ((v0[0] - x) * sy - (v0[1] - y) * sx) * d;
  if (w[2] < 0. || w[1] + w[2] > 1.) return false;

  // Weighting factor for point A
  w[0] = 1. - w[1] - w[2];

  return true;
}

bool ComponentTcad2d::OnLine(const double x, const double y,
                             const Element& element,
                             std::array<double, nMaxVertices>& w) const {
  const auto& v0 = m_vertices[element.vertex[0]];
  const auto& v1 = m_vertices[element.vertex[1]];
  if (x > v1[0]) return false;
  if (y < v0[1] && y < v1[1]) return false;
  if (y > v0[1] && y > v1[1]) return false;
  const double tx = (x - v0[0]) / (v1[0] - v0[0]);
  if (tx < 0. || tx > 1.) return false;
  const double ty = (y - v0[1]) / (v1[1] - v0[1]);
  if (ty < 0. || ty > 1.) return false;
  if (tx == ty) {
    // Compute weighting factors for endpoints A, B
    w[0] = tx;
    w[1] = 1. - w[0];
    return true;
  }
  return false;
}

bool ComponentTcad2d::AtPoint(const double x, const double y,
                              const Element& element,
                              std::array<double, nMaxVertices>& w) const {
  const auto& v0 = m_vertices[element.vertex[0]];
  if (x != v0[0] || y != v0[1]) return false;
  w[0] = 1;
  return true;
}

}
