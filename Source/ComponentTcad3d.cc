#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "Garfield/ComponentTcad3d.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Utilities.hh"

namespace Garfield {

ComponentTcad3d::ComponentTcad3d() : Component("Tcad3d") {
  m_regions.reserve(10);
  m_vertices.reserve(10000);
  m_elements.reserve(10000);
}

void ComponentTcad3d::ElectricField(const double xin, const double yin,
                                    const double zin, double& ex, double& ey,
                                    double& ez, double& p, Medium*& m,
                                    int& status) {
  ex = ey = ez = p = 0.;
  m = nullptr;
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::ElectricField:\n"
              << "    Field map is not available for interpolation.\n";
    status = -10;
    return;
  }

  double x = xin, y = yin, z = zin;
  // In case of periodicity, reduce to the cell volume.
  bool xmirr = false, ymirr = false, zmirr = false;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr);

  // Check if the point is inside the bounding box.
  if (x < m_xMinBB || x > m_xMaxBB || y < m_yMinBB || y > m_yMaxBB ||
      z < m_zMinBB || z > m_zMaxBB) {
    status = -6;
    return;
  }

  // Assume this will work.
  status = 0;
  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x, y, z, w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    status = -6;
    return;
  }
  const Element& element = m_elements[i];
  const unsigned int nVertices = element.type == 2 ? 3 : 4;
  for (unsigned int j = 0; j < nVertices; ++j) {
    const Vertex& vj = m_vertices[element.vertex[j]];
    ex += w[j] * vj.ex;
    ey += w[j] * vj.ey;
    ez += w[j] * vj.ez;
    p += w[j] * vj.p;
  }
  if (xmirr) ex = -ex;
  if (ymirr) ey = -ey;
  if (zmirr) ez = -ez;
  const auto& region = m_regions[element.region];
  m = region.medium;
  if (!region.drift || !m) status = -5;
}

void ComponentTcad3d::ElectricField(const double x, const double y,
                                    const double z, double& ex, double& ey,
                                    double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

void ComponentTcad3d::WeightingField(const double xin, const double yin,
                                     const double zin, double& wx, double& wy,
                                     double& wz, const std::string& /*label*/) {
  wx = wy = wz = 0.;
  if (m_wf.empty()) {
    std::cerr << m_className << "::WeightingField: Not available.\n";
    return;
  }
  double x = xin, y = yin, z = zin;
  // In case of periodicity, reduce to the cell volume.
  bool xmirr = false, ymirr = false, zmirr = false;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr);

  // Check if the point is inside the bounding box.
  if (x < m_xMinBB || x > m_xMaxBB || y < m_yMinBB || y > m_yMaxBB ||
      z < m_zMinBB || z > m_zMaxBB) return;

  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x, y, z, w);
  if (i >= m_elements.size()) return;

  const Element& element = m_elements[i];
  const unsigned int nVertices = element.type == 2 ? 3 : 4;
  for (unsigned int j = 0; j < nVertices; ++j) {
    const auto& f = m_wf[element.vertex[j]];
    wx += w[j] * f[0];
    wy += w[j] * f[1];
    wz += w[j] * f[2];
  }
  if (xmirr) wx = -wx;
  if (ymirr) wy = -wy;
  if (zmirr) wz = -wz;
}

double ComponentTcad3d::WeightingPotential(const double xin, const double yin,
                                           const double zin, 
                                           const std::string& /*label*/) {

  if (m_wp.empty()) {
    std::cerr << m_className << "::WeightingPotential: Not available.\n";
    return 0.;
  }
  // In case of periodicity, reduce to the cell volume.
  double x = xin, y = yin, z = zin;
  bool xmirr = false, ymirr = false, zmirr = false;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr);

  // Check if the point is inside the bounding box.
  if (x < m_xMinBB || x > m_xMaxBB || y < m_yMinBB || y > m_yMaxBB ||
      z < m_zMinBB || z > m_zMaxBB) return 0.;

  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x, y, z, w);
  if (i >= m_elements.size()) return 0.;

  double v = 0.;
  const Element& element = m_elements[i];
  const unsigned int nVertices = element.type == 2 ? 3 : 4;
  for (unsigned int j = 0; j < nVertices; ++j) {
    v += w[j] * m_wp[element.vertex[j]];
  }
  return v;
}

Medium* ComponentTcad3d::GetMedium(const double xin, const double yin,
                                   const double zin) {
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::GetMedium:\n"
              << "    Field map not available for interpolation.\n";
    return nullptr;
  }

  double x = xin, y = yin, z = zin;
  bool xmirr = false, ymirr = false, zmirr = false;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr);

  // Check if the point is inside the bounding box.
  if (x < m_xMinBB || x > m_xMaxBB || y < m_yMinBB || y > m_yMaxBB ||
      z < m_zMinBB || z > m_zMaxBB) {
    return nullptr;
  }

  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x, y, z, w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    return nullptr;
  }
  const Element& element = m_elements[i];
  return m_regions[element.region].medium;
}

bool ComponentTcad3d::Initialise(const std::string& gridfilename,
                                 const std::string& datafilename) {
  m_ready = false;
  // Import mesh data from .grd file.
  if (!LoadGrid(gridfilename)) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Importing mesh data failed.\n";
    return false;
  }

  // Import electric field and potential from .dat file.
  if (!LoadData(datafilename)) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Importing electric field and potential failed.\n";
    return false;
  }

  // Find min./max. coordinates and potentials.
  m_xMaxBB = m_vertices[m_elements[0].vertex[0]].x;
  m_yMaxBB = m_vertices[m_elements[0].vertex[0]].y;
  m_zMaxBB = m_vertices[m_elements[0].vertex[0]].z;
  m_xMinBB = m_xMaxBB;
  m_yMinBB = m_yMaxBB;
  m_zMinBB = m_zMaxBB;
  m_pMax = m_pMin = m_vertices[m_elements[0].vertex[0]].p;
  const size_t nElements = m_elements.size();
  for (size_t i = 0; i < nElements; ++i) {
    Element& element = m_elements[i];
    double xmin = m_vertices[element.vertex[0]].x;
    double ymin = m_vertices[element.vertex[0]].y;
    double zmin = m_vertices[element.vertex[0]].z;
    double xmax = xmin;
    double ymax = ymin;
    double zmax = zmin;
    const size_t nVertices = element.type == 2 ? 3 : 4;
    for (size_t j = 0; j < nVertices; ++j) {
      const Vertex& vj = m_vertices[element.vertex[j]];
      if (vj.x < xmin) xmin = vj.x;
      if (vj.x > xmax) xmax = vj.x;
      if (vj.y < ymin) ymin = vj.y;
      if (vj.y > ymax) ymax = vj.y;
      if (vj.z < zmin) zmin = vj.z;
      if (vj.z > zmax) zmax = vj.z;
      if (vj.p < m_pMin) m_pMin = vj.p;
      if (vj.p > m_pMax) m_pMax = vj.p;
    }
    constexpr double tol = 1.e-6;
    element.xmin = xmin - tol;
    element.xmax = xmax + tol;
    element.ymin = ymin - tol;
    element.ymax = ymax + tol;
    element.zmin = zmin - tol;
    element.zmax = zmax + tol;
    m_xMinBB = std::min(m_xMinBB, xmin);
    m_xMaxBB = std::max(m_xMaxBB, xmax);
    m_yMinBB = std::min(m_yMinBB, ymin);
    m_yMaxBB = std::max(m_yMaxBB, ymax);
    m_zMinBB = std::min(m_zMinBB, zmin);
    m_zMaxBB = std::max(m_zMaxBB, zmax);
  }

  std::cout << m_className << "::Initialise:\n"
            << "    Bounding box:\n"
            << "      " << m_xMinBB << " < x [cm] < " << m_xMaxBB << "\n"
            << "      " << m_yMinBB << " < y [cm] < " << m_yMaxBB << "\n"
            << "      " << m_zMinBB << " < z [cm] < " << m_zMaxBB << "\n"
            << "    Voltage range:\n"
            << "      " << m_pMin << " < V < " << m_pMax << "\n";

  bool ok = true;

  // Count the number of elements belonging to a region.
  const size_t nRegions = m_regions.size();
  std::vector<size_t> nElementsRegion(nRegions, 0);

  // Count the different element shapes.
  size_t nTriangles = 0;
  size_t nTetrahedra = 0;
  size_t nOtherShapes = 0;

  // Check if there are elements which are not part of any region.
  std::vector<size_t> looseElements;
  // Check if there are degenerate elements.
  std::vector<size_t> degenerateElements;

  for (size_t i = 0; i < nElements; ++i) {
    const Element& element = m_elements[i];
    if (element.type == 2) {
      ++nTriangles;
      if (element.vertex[0] == element.vertex[1] ||
          element.vertex[1] == element.vertex[2] ||
          element.vertex[2] == element.vertex[0]) {
        degenerateElements.push_back(i);
      }
    } else if (element.type == 5) {
      if (element.vertex[0] == element.vertex[1] ||
          element.vertex[0] == element.vertex[2] ||
          element.vertex[0] == element.vertex[3] ||
          element.vertex[1] == element.vertex[2] ||
          element.vertex[1] == element.vertex[3] ||
          element.vertex[2] == element.vertex[3]) {
        degenerateElements.push_back(i);
      }
      ++nTetrahedra;
    } else {
      // Other shapes should not occur, since they were excluded in LoadGrid.
      ++nOtherShapes;
    }
    if (element.region >= 0 && element.region < (int)nRegions) {
      ++nElementsRegion[element.region];
    } else {
      looseElements.push_back(i);
    }
  }

  if (!degenerateElements.empty()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    The following elements are degenerate:\n";
    for (size_t i : degenerateElements) std::cerr << "      " << i << "\n";
    ok = false;
  }

  if (!looseElements.empty()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    The following elements are not part of any region:\n";
    for (size_t i : looseElements) std::cerr << "      " << i << "\n";
    ok = false;
  }

  std::cout << m_className << "::Initialise:\n"
            << "    Number of regions: " << nRegions << "\n";
  for (size_t i = 0; i < nRegions; ++i) {
    std::cout << "      " << i << ": " << m_regions[i].name << ", "
              << nElementsRegion[i] << " elements\n";
  }

  std::cout << "    Number of elements: " << nElements << "\n";
  if (nTriangles > 0) {
    std::cout << "      " << nTriangles << " triangles\n";
  }
  if (nTetrahedra > 0) {
    std::cout << "      " << nTetrahedra << " tetrahedra\n";
  }
  if (nOtherShapes > 0) {
    std::cerr << "      " << nOtherShapes << " elements of unknown type\n"
              << "      Program bug!\n";
    m_ready = false;
    Cleanup();
    return false;
  }
  if (m_debug) {
    // For each element, print the indices of the constituting vertices.
    for (size_t i = 0; i < nElements; ++i) {
      const Element& element = m_elements[i];
      if (element.type == 2) {
        std::cout << "      " << i << ": " << element.vertex[0] << "  "
                  << element.vertex[1] << "  " << element.vertex[2]
                  << " (triangle, region " << element.region << ")\n";
      } else if (element.type == 5) {
        std::cout << "      " << i << ": " << element.vertex[0] << "  "
                  << element.vertex[1] << "  " << element.vertex[2] << "  "
                  << element.vertex[3] << " (tetrahedron, region "
                  << element.region << ")\n";
      }
    }
  }

  const size_t nVertices = m_vertices.size();
  std::cout << "    Number of vertices: " << nVertices << "\n";
  if (m_debug) {
    for (size_t i = 0; i < nVertices; ++i) {
      const Vertex& vi = m_vertices[i];
      std::cout << "      " << i << ": (x, y, z) = (" << vi.x << ", " << vi.y
                << ", " << vi.z << "), V = " << vi.p << "\n";
    }
  }

  if (!ok) {
    m_ready = false;
    Cleanup();
    return false;
  }

  // Set up the octree.
  const float hx = 0.5 * (m_xMaxBB - m_xMinBB);
  const float hy = 0.5 * (m_yMaxBB - m_yMinBB);
  const float hz = 0.5 * (m_zMaxBB - m_zMinBB);
  m_tree.reset(new TetrahedralTree(Vec3(m_xMinBB + hx, m_yMinBB + hy, m_zMinBB + hz),
                                   Vec3(hx, hy, hz)));

  // Insert the mesh nodes in the tree.
  for (size_t i = 0; i < nVertices; ++i) {
    const Vertex& vtx = m_vertices[i];
    m_tree->InsertMeshNode(Vec3(vtx.x, vtx.y, vtx.z), i);
  }

  // Insert the mesh elements in the tree.
  for (size_t i = 0; i < nElements; ++i) {
    const Element& e = m_elements[i];
    const double bb[6] = {e.xmin, e.ymin, e.zmin, e.xmax, e.ymax, e.zmax};
    m_tree->InsertMeshElement(bb, i);
  }

  m_ready = true;
  UpdatePeriodicity();
  std::cout << m_className << "::Initialise: Initialisation finished.\n";
  return true;
}

bool ComponentTcad3d::SetWeightingField(const std::string& datfile1,
                                        const std::string& datfile2,
                                        const double dv) {

  if (!m_ready) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Mesh is not available. Call Initialise first.\n";
    return false;
  }
  if (dv < Small) {
     std::cerr << m_className << "::SetWeightingField:\n"
               << "    Voltage difference must be > 0.\n";
     return false;
  }
  const double s = 1. / dv;

  m_wf.clear();
  m_wp.clear();
  // Load first the field/potential at nominal bias.
  std::vector<std::array<double, 3> > wf1;
  std::vector<double> wp1;
  if (!LoadWeightingField(datfile1, wf1, wp1)) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Could not import data from " << datfile1 << ".\n";
    return false;
  }
  // Then load the field/potential for the configuration with the potential 
  // at the electrode to be read out increased by small voltage dv. 
  std::vector<std::array<double, 3> > wf2;
  std::vector<double> wp2;
  if (!LoadWeightingField(datfile2, wf2, wp2)) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Could not import data from " << datfile2 << ".\n";
    return false;
  }
  const size_t nVertices = m_vertices.size();
  bool foundField = true;
  if (wf1.size() != nVertices || wf2.size() != nVertices) {
    foundField = false;
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Could not load electric field values.\n";
  }
  bool foundPotential = true;
  if (wp1.size() != nVertices || wp2.size() != nVertices) {
    foundPotential = false;
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Could not load electrostatic potentials.\n";
  }
  if (!foundField && !foundPotential) return false;
  if (foundField) {
    m_wf.assign(nVertices, {0., 0., 0.});
    for (unsigned int i = 0; i < nVertices; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
        m_wf[i][j] = (wf2[i][j] - wf1[i][j]) * s;
      } 
    }
  }
  if (foundPotential) {
    m_wp.assign(nVertices, 0.);
    for (unsigned int i = 0; i < nVertices; ++i) {
      m_wp[i] = (wp2[i] - wp1[i]) * s; 
    }
  }
  return true;
}

bool ComponentTcad3d::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                     double& xmax, double& ymax, double& zmax) {
  if (!m_ready) return false;
  xmin = m_xMinBB;
  ymin = m_yMinBB;
  zmin = m_zMinBB;
  xmax = m_xMaxBB;
  ymax = m_yMaxBB;
  zmax = m_zMaxBB;
  if (m_periodic[0] || m_mirrorPeriodic[0]) {
    xmin = -std::numeric_limits<double>::infinity();
    xmax = +std::numeric_limits<double>::infinity();
  }
  if (m_periodic[1] || m_mirrorPeriodic[1]) {
    ymin = -std::numeric_limits<double>::infinity();
    ymax = +std::numeric_limits<double>::infinity();
  }
  if (m_periodic[2] || m_mirrorPeriodic[2]) {
    zmin = -std::numeric_limits<double>::infinity();
    zmax = +std::numeric_limits<double>::infinity();
  }
  return true;
}

bool ComponentTcad3d::GetVoltageRange(double& vmin, double& vmax) {
  if (!m_ready) return false;
  vmin = m_pMin;
  vmax = m_pMax;
  return true;
}

void ComponentTcad3d::PrintRegions() {
  // Do not proceed if not properly initialised.
  if (!m_ready) {
    std::cerr << m_className << "::PrintRegions:\n"
              << "    Field map not yet initialised.\n";
    return;
  }

  if (m_regions.empty()) {
    std::cerr << m_className << "::PrintRegions:\n"
              << "    No regions are currently defined.\n";
    return;
  }

  const size_t nRegions = m_regions.size();
  std::cout << m_className << "::PrintRegions:\n"
            << "    Currently " << nRegions << " regions are defined.\n"
            << "      Index  Name      Medium\n";
  for (size_t i = 0; i < nRegions; ++i) {
    std::cout << "      " << i << "  " << m_regions[i].name;
    if (!m_regions[i].medium) {
      std::cout << "      none  ";
    } else {
      std::cout << "      " << m_regions[i].medium->GetName();
    }
    if (m_regions[i].drift) {
      std::cout << " (active region)\n";
    } else {
      std::cout << "\n";
    }
  }
}

void ComponentTcad3d::GetRegion(const size_t i, std::string& name,
                                bool& active) const {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::GetRegion: Index out of range.\n";
    return;
  }
  name = m_regions[i].name;
  active = m_regions[i].drift;
}

void ComponentTcad3d::SetDriftRegion(const size_t i) {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::SetDriftRegion: Index out of range.\n";
    return;
  }
  m_regions[i].drift = true;
}

void ComponentTcad3d::UnsetDriftRegion(const size_t i) {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::UnsetDriftRegion: Index out of range.\n";
    return;
  }
  m_regions[i].drift = false;
}

void ComponentTcad3d::SetMedium(const size_t i, Medium* medium) {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::SetMedium: Index out of range.\n";
    return;
  }

  if (!medium) {
    std::cerr << m_className << "::SetMedium: Null pointer.\n";
    return;
  }
  m_regions[i].medium = medium;
}

bool ComponentTcad3d::GetMedium(const size_t i, Medium*& m) const {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::GetMedium: Index out of range.\n";
    return false;
  }

  m = m_regions[i].medium;
  if (!m) return false;
  return true;
}

size_t ComponentTcad3d::FindElement(
    const double x, const double y, const double z,
    std::array<double, nMaxVertices>& w) const {

  w.fill(0.);

  std::vector<int> elementsToSearch;
  if (m_tree) elementsToSearch = m_tree->GetElementsInBlock(Vec3(x, y, z));
  const size_t nElementsToSearch = m_tree ? elementsToSearch.size() : m_elements.size(); 
  // Loop over the elements.
  for (size_t i = 0; i < nElementsToSearch; ++i) {
    const size_t idx = m_tree ? elementsToSearch[i] : i;
    const Element& element = m_elements[idx];
    if (x < element.xmin || x > element.xmax || y < element.ymin ||
        y > element.ymax || z < element.zmin || z > element.zmax)
      continue;
    if (CheckElement(x, y, z, element, w)) return idx;
  }

  if (m_debug) {
    std::cerr << m_className << "::FindElement:\n"
              << "    Point (" << x << ", " << y << ", " << z
              << ") is outside the mesh.\n";
  }
  return m_elements.size();
}

bool ComponentTcad3d::GetElement(const size_t i, double& vol,
                                 double& dmin, double& dmax, int& type) const {
  if (i >= m_elements.size()) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }

  const Element& element = m_elements[i];
  if (element.type == 2) {
    // Triangle
    const Vertex& v0 = m_vertices[element.vertex[0]];
    const Vertex& v1 = m_vertices[element.vertex[1]];
    const Vertex& v2 = m_vertices[element.vertex[2]];
    const double vx =
        (v1.y - v0.y) * (v2.z - v0.z) - (v1.z - v0.z) * (v2.y - v0.y);
    const double vy =
        (v1.z - v0.z) * (v2.x - v0.x) - (v1.x - v0.x) * (v2.z - v0.z);
    const double vz =
        (v1.x - v0.x) * (v2.y - v0.y) - (v1.y - v0.y) * (v2.x - v0.x);
    vol = sqrt(vx * vx + vy * vy + vz * vz);
    const double a =
        sqrt(pow(v1.x - v0.x, 2) + pow(v1.y - v0.y, 2) + pow(v1.z - v0.z, 2));
    const double b =
        sqrt(pow(v2.x - v0.x, 2) + pow(v2.y - v0.y, 2) + pow(v2.z - v0.z, 2));
    const double c =
        sqrt(pow(v1.x - v2.x, 2) + pow(v1.y - v2.y, 2) + pow(v1.z - v2.z, 2));
    dmin = std::min({a, b, c});
    dmax = std::max({a, b, c});
  } else if (element.type == 5) {
    // Tetrahedron
    const Vertex& v0 = m_vertices[element.vertex[0]];
    const Vertex& v1 = m_vertices[element.vertex[1]];
    const Vertex& v2 = m_vertices[element.vertex[2]];
    const Vertex& v3 = m_vertices[element.vertex[3]];
    vol = fabs((v3.x - v0.x) * ((v1.y - v0.y) * (v2.z - v0.z) -
                                (v2.y - v0.y) * (v1.z - v0.z)) +
               (v3.y - v0.y) * ((v1.z - v0.z) * (v2.x - v0.x) -
                                (v2.z - v0.z) * (v1.x - v0.x)) +
               (v3.z - v0.z) * ((v1.x - v0.x) * (v2.y - v0.y) -
                                (v3.x - v0.x) * (v1.y - v0.y))) /
          6.;
    // Loop over all pairs of m_vertices.
    for (size_t j = 0; j < nMaxVertices - 1; ++j) {
      const Vertex& vj = m_vertices[element.vertex[j]];
      for (size_t k = j + 1; k < nMaxVertices; ++k) {
        const Vertex& vk = m_vertices[element.vertex[k]];
        // Compute distance.
        const double dist = sqrt(pow(vj.x - vk.x, 2) + pow(vj.y - vk.y, 2) +
                                 pow(vj.z - vk.z, 2));
        if (k == 1) {
          dmin = dmax = dist;
        } else {
          if (dist < dmin) dmin = dist;
          if (dist > dmax) dmax = dist;
        }
      }
    }
  } else {
    std::cerr << m_className << "::GetElement:\n"
              << "    Unexpected element type (" << type << ").\n";
    return false;
  }
  return true;
}

bool ComponentTcad3d::GetElement(const size_t i, double& vol, 
                                 double& dmin, double& dmax, int& type,
                                 std::vector<size_t>& nodes,
                                 int& reg) const {
  if (!GetElement(i, vol, dmin, dmax, type)) return false;
  const Element& element = m_elements[i];
  const size_t nVertices = element.type + 1;
  for (size_t j = 0; j < nVertices; ++j) {
    nodes.push_back(element.vertex[j]);
  }
  reg = element.region;
  return true;
}

bool ComponentTcad3d::GetNode(const size_t i, double& x, double& y,
                              double& z, double& v, double& ex, double& ey,
                              double& ez) const {
  if (i >= m_vertices.size()) {
    std::cerr << m_className << "::GetNode: Index out of range.\n";
    return false;
  }

  const Vertex& vi = m_vertices[i];
  x = vi.x;
  y = vi.y;
  z = vi.z;
  v = vi.p;
  ex = vi.ex;
  ey = vi.ey;
  ez = vi.ez;
  return true;
}

bool ComponentTcad3d::LoadData(const std::string& datafilename) {
  std::ifstream datafile;
  datafile.open(datafilename.c_str(), std::ios::in);
  if (!datafile) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Could not open file " << datafilename << ".\n";
    return false;
  }

  const size_t nVertices = m_vertices.size();
  for (size_t i = 0; i < nVertices; ++i) {
    m_vertices[i].p = 0.;
    m_vertices[i].ex = 0.;
    m_vertices[i].ey = 0.;
    m_vertices[i].ez = 0.;
  }
  // Read the file line by line.
  std::string line;
  while (std::getline(datafile, line)) {
    // Strip white space from beginning of line.
    ltrim(line);
    // Find data section.
    if (line.substr(0, 8) != "function") continue;
    // Read type of data set.
    const auto pEq = line.find('=');
    if (pEq == std::string::npos) {
      // No "=" found.
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading file " << datafilename << ".\n"
                << "    Line:\n    " << line << "\n";
      datafile.close();
      Cleanup();
      return false;
    }
    line = line.substr(pEq + 1);
    std::string dataset;
    std::istringstream data;
    data.str(line);
    data >> dataset;
    data.clear();
    if (dataset == "ElectrostaticPotential") {
      if (!ReadDataset(datafile, dataset)) return false;
    } else if (dataset == "ElectricField") {
      if (!ReadDataset(datafile, dataset)) return false;
    }
  }
  if (datafile.fail() && !datafile.eof()) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Error reading file " << datafilename << "\n";
    datafile.close();
    Cleanup();
    return false;
  }
  datafile.close();
  return true;
}

bool ComponentTcad3d::ReadDataset(std::ifstream& datafile,
                                  const std::string& dataset) {
  if (!datafile.is_open()) return false;
  enum DataSet { ElectrostaticPotential, EField, Unknown };
  DataSet ds = Unknown;
  if (dataset == "ElectrostaticPotential") {
    ds = ElectrostaticPotential;
  } else if (dataset == "ElectricField") {
    ds = EField;
  } else {
    std::cerr << m_className << "::ReadDataset:\n"
              << "    Unexpected dataset " << dataset << ".\n";
    return false;
  }
  bool isVector = false;
  if (ds == EField) {
    isVector = true;
  }

  std::string line;
  std::getline(datafile, line);
  std::getline(datafile, line);
  std::getline(datafile, line);
  std::getline(datafile, line);
  // Get the region name (given in brackets).
  auto bra = line.find('[');
  auto ket = line.find(']');
  if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
    std::cerr << m_className << "::ReadDataset:\n"
              << "    Cannot extract region name.\n"
              << "    Line:\n    " << line << "\n";
    datafile.close();
    Cleanup();
    return false;
  }
  line = line.substr(bra + 1, ket - bra - 1);
  std::string name;
  std::istringstream data;
  data.str(line);
  data >> name;
  data.clear();
  // Check if the region name matches one from the mesh file.
  const size_t index = FindRegion(name);
  if (index >= m_regions.size()) {
    std::cerr << m_className << "::ReadDataset:\n"
              << "    Unknown region " << name << ".\n";
    datafile.close();
    Cleanup();
    return false;
  }
  // Get the number of values.
  std::getline(datafile, line);
  bra = line.find('(');
  ket = line.find(')');
  if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
    std::cerr << m_className << "::ReadDataset:\n"
              << "    Cannot extract number of values to be read.\n"
              << "    Line:\n    " << line << "\n";
    datafile.close();
    Cleanup();
    return false;
  }
  line = line.substr(bra + 1, ket - bra - 1);
  int nValues;
  data.str(line);
  data >> nValues;
  if (isVector) nValues /= 3;
  // Mark the vertices belonging to this region.
  const size_t nVertices = m_vertices.size();
  std::vector<bool> isInRegion(nVertices, false);
  const size_t nElements = m_elements.size();
  for (size_t j = 0; j < nElements; ++j) {
    if (m_elements[j].region != (int)index) continue;
    for (int k = 0; k <= m_elements[j].type; ++k) {
      isInRegion[m_elements[j].vertex[k]] = true;
    }
  }

  unsigned int ivertex = 0;
  for (int j = 0; j < nValues; ++j) {
    // Read the next value.
    double val1 = 0., val2 = 0., val3 = 0.;
    if (isVector) {
      datafile >> val1 >> val2 >> val3;
    } else {
      datafile >> val1;
    }
    // Find the next vertex belonging to the region.
    while (ivertex < nVertices) {
      if (isInRegion[ivertex]) break;
      ++ivertex;
    }
    // Check if there is a mismatch between the number of m_vertices
    // and the number of potential values.
    if (ivertex >= nVertices) {
      std::cerr << m_className << "::ReadDataset:\n"
                << "    Dataset " << dataset << " has more values than "
                << "there are vertices in region " << name << "\n";
      datafile.close();
      Cleanup();
      return false;
    }
    switch (ds) {
      case ElectrostaticPotential:
        m_vertices[ivertex].p = val1;
        break;
      case EField:
        m_vertices[ivertex].ex = val1;
        m_vertices[ivertex].ey = val2;
        m_vertices[ivertex].ez = val3;
        break;
      default:
        std::cerr << m_className << "::ReadDataset:\n"
                  << "    Unexpected dataset (" << ds << "). Program bug!\n";
        datafile.close();
        Cleanup();
        return false;
    }
    ++ivertex;
  }
  return true;
}

bool ComponentTcad3d::LoadWeightingField(const std::string& datafilename,
    std::vector<std::array<double, 3> >& wf, std::vector<double>& wp) {

  std::ifstream datafile;
  datafile.open(datafilename.c_str(), std::ios::in);
  if (!datafile) {
    std::cerr << m_className << "::LoadWeightingField:\n"
              << "    Could not open file " << datafilename << ".\n";
    return false;
  }

  const size_t nVertices = m_vertices.size();
  bool ok = true;
  // Read the file line by line.
  std::string line;
  while (std::getline(datafile, line)) {
    // Strip white space from beginning of line.
    ltrim(line);
    // Find data section.
    if (line.substr(0, 8) != "function") continue;
    // Read type of data set.
    const auto pEq = line.find('=');
    if (pEq == std::string::npos) {
      // No "=" found.
      std::cerr << m_className << "::LoadWeightingField:\n"
                << "    Error reading file " << datafilename << ".\n"
                << "    Line:\n    " << line << "\n";
      datafile.close();
      return false;
    }
    line = line.substr(pEq + 1);
    std::string dataset;
    std::istringstream data;
    data.str(line);
    data >> dataset;
    data.clear();
    if (dataset != "ElectrostaticPotential" && 
        dataset != "ElectricField") continue;
    bool field = false;
    if (dataset == "ElectricField") {
      wf.assign(nVertices, {0., 0., 0.});
      field = true;
    } else {
      wp.assign(nVertices, 0.);
    }
    std::getline(datafile, line);
    std::getline(datafile, line);
    std::getline(datafile, line);
    std::getline(datafile, line);
    // Get the region name (given in brackets).
    auto bra = line.find('[');
    auto ket = line.find(']');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      std::cerr << m_className << "::LoadWeightingField:\n"
                << "    Cannot extract region name.\n"
                << "    Line:\n    " << line << "\n";
      ok = false;
      break;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    std::string name;
    data.str(line);
    data >> name;
    data.clear();
    // Check if the region name matches one from the mesh file.
    const size_t index = FindRegion(name);
    if (index >= m_regions.size()) {
      std::cerr << m_className << "::LoadWeightingField:\n"
                << "    Unknown region " << name << ".\n";
      ok = false;
      break;
    }
    // Get the number of values.
    std::getline(datafile, line);
    bra = line.find('(');
    ket = line.find(')');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      std::cerr << m_className << "::LoadWeightingField:\n"
                << "    Cannot extract number of values to be read.\n"
                << "    Line:\n    " << line << "\n";
      ok = false;
      break;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    int nValues;
    data.str(line);
    data >> nValues;
    if (field) nValues /= 3;
    // Mark the vertices belonging to this region.
    std::vector<bool> isInRegion(nVertices, false);
    const size_t nElements = m_elements.size();
    for (size_t j = 0; j < nElements; ++j) {
      if (m_elements[j].region != (int)index) continue;
      for (int k = 0; k <= m_elements[j].type; ++k) {
        isInRegion[m_elements[j].vertex[k]] = true;
      }
    }
    unsigned int ivertex = 0;
    for (int j = 0; j < nValues; ++j) {
      // Read the next value.
      double val1 = 0., val2 = 0., val3 = 0.;
      if (field) {
        datafile >> val1 >> val2 >> val3;
      } else {
        datafile >> val1;
      }
      // Find the next vertex belonging to the region.
      while (ivertex < nVertices) {
        if (isInRegion[ivertex]) break;
        ++ivertex;
      }
      // Check if there is a mismatch between the number of vertices
      // and the number of potential values.
      if (ivertex >= nVertices) {
        std::cerr << m_className << "::LoadWeightingField:\n"
                  << "    Dataset " << dataset << " has more values than "
                  << "there are vertices in region " << name << "\n";
        ok = false;
        break;
      }
      if (field) {
        wf[ivertex] = {val1, val2, val3};
      } else {
        wp[ivertex] = val1;
      }
      ++ivertex;
    }
  }

  if (!ok || (datafile.fail() && !datafile.eof())) {
    std::cerr << m_className << "::LoadWeightingField:\n"
              << "    Error reading file " << datafilename << "\n";
    datafile.close();
    return false;
  }
  datafile.close();
  return true;
}

bool ComponentTcad3d::LoadGrid(const std::string& gridfilename) {
  // Open the file containing the mesh description.
  std::ifstream gridfile;
  gridfile.open(gridfilename.c_str(), std::ios::in);
  if (!gridfile) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not open file " << gridfilename << ".\n";
    return false;
  }

  // Delete existing mesh information.
  Cleanup();

  // Count line numbers.
  int iLine = 0;

  // Get the number of regions.
  size_t nRegions = 0;
  // Read the file line by line.
  std::string line;
  while (std::getline(gridfile, line)) {
    ++iLine;
    // Strip white space from the beginning of the line.
    ltrim(line);
    // Find entry 'nb_regions'.
    if (line.substr(0, 10) != "nb_regions") continue;
    const auto pEq = line.find('=');
    if (pEq == std::string::npos) {
      // No "=" sign found.
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of regions.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    line = line.substr(pEq + 1);
    std::istringstream data;
    data.str(line);
    data >> nRegions;
    break;
  }
  if (gridfile.eof()) {
    // Reached end of file.
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find entry 'nb_regions' in file\n"
              << "    " << gridfilename << ".\n";
    Cleanup();
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    // Error reading from the file.
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << gridfilename << " (line " << iLine
              << ").\n";
    Cleanup();
    gridfile.close();
    return false;
  }
  m_regions.resize(nRegions);
  for (size_t j = 0; j < nRegions; ++j) {
    m_regions[j].name = "";
    m_regions[j].drift = false;
    m_regions[j].medium = nullptr;
  }

  if (m_debug) {
    std::cout << m_className << "::LoadGrid:\n"
              << "    Found " << nRegions << " regions.\n";
  }

  // Get the region names.
  while (std::getline(gridfile, line)) {
    ++iLine;
    ltrim(line);
    // Find entry 'regions'.
    if (line.substr(0, 7) != "regions") continue;
    // Get region names (given in brackets).
    const auto bra = line.find('[');
    const auto ket = line.find(']');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      // No closed brackets [].
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read region names.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    std::istringstream data;
    data.str(line);
    for (size_t j = 0; j < nRegions; ++j) {
      data >> m_regions[j].name;
      data.clear();
      // Assume by default that all regions are active.
      m_regions[j].drift = true;
      m_regions[j].medium = nullptr;
    }
    break;
  }
  if (gridfile.eof()) {
    // Reached end of file.
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find entry 'regions' in file\n"
              << "    " << gridfilename << ".\n";
    Cleanup();
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    // Error reading from the file.
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << gridfilename << " (line " << iLine
              << ").\n";
    Cleanup();
    gridfile.close();
    return false;
  }

  // Get the vertices.
  size_t nVertices = 0;
  while (std::getline(gridfile, line)) {
    ++iLine;
    ltrim(line);
    // Find section 'Vertices'.
    if (line.substr(0, 8) != "Vertices") continue;
    // Get number of vertices (given in brackets).
    const auto bra = line.find('(');
    const auto ket = line.find(')');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      // No closed brackets [].
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of vertices.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    std::istringstream data;
    data.str(line);
    data >> nVertices;
    m_vertices.resize(nVertices);
    // Get the coordinates of this vertex.
    for (size_t j = 0; j < nVertices; ++j) {
      gridfile >> m_vertices[j].x >> m_vertices[j].y >> m_vertices[j].z;
      // Change units from micron to cm.
      m_vertices[j].x *= 1.e-4;
      m_vertices[j].y *= 1.e-4;
      m_vertices[j].z *= 1.e-4;
    }
    iLine += nVertices - 1;
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Vertices' in file\n"
              << "    " << gridfilename << ".\n";
    Cleanup();
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << gridfilename << " (line " << iLine
              << ").\n";
    Cleanup();
    gridfile.close();
    return false;
  }

  // Get the "edges" (lines connecting two vertices).
  size_t nEdges = 0;
  // Temporary arrays for storing edge points.
  std::vector<size_t> edgeP1;
  std::vector<size_t> edgeP2;
  while (std::getline(gridfile, line)) {
    ++iLine;
    ltrim(line);
    // Find section 'Edges'.
    if (line.substr(0, 5) != "Edges") continue;
    // Get the number of edges (given in brackets).
    const auto bra = line.find('(');
    const auto ket = line.find(')');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      // No closed brackets ()
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of edges.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    std::istringstream data;
    data.str(line);
    data >> nEdges;
    edgeP1.resize(nEdges);
    edgeP2.resize(nEdges);
    // Get the indices of the two endpoints.
    for (size_t j = 0; j < nEdges; ++j) {
      gridfile >> edgeP1[j] >> edgeP2[j];
    }
    iLine += nEdges - 1;
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Edges' in file\n"
              << "    " << gridfilename << ".\n";
    Cleanup();
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << gridfilename << " (line " << iLine
              << ").\n";
    Cleanup();
    gridfile.close();
    return false;
  }

  for (size_t i = 0; i < nEdges; ++i) {
    // Make sure the indices of the edge endpoints are not out of range.
    if (edgeP1[i] >= nVertices || edgeP2[i] >= nVertices) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Vertex index of edge " << i << " out of range.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    // Make sure the edge is non-degenerate.
    if (edgeP1[i] == edgeP2[i]) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Edge " << i << " is degenerate.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
  }

  // Get the "faces".
  size_t nFaces = 0;
  std::vector<Face> faces;
  while (std::getline(gridfile, line)) {
    ++iLine;
    ltrim(line);
    // Find section 'Faces'.
    if (line.substr(0, 5) != "Faces") continue;
    // Get the number of faces (given in brackets).
    const auto bra = line.find('(');
    const auto ket = line.find(')');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      // No closed brackets ()
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of faces.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    std::istringstream data;
    data.str(line);
    data >> nFaces;
    faces.resize(nFaces);
    // Get the indices of the edges constituting this face.
    for (size_t j = 0; j < nFaces; ++j) {
      gridfile >> faces[j].type;
      if (faces[j].type != 3 && faces[j].type != 4) {
        std::cerr << m_className << "::LoadGrid:\n"
                  << "    Face with index " << j
                  << " has invalid number of edges, " << faces[j].type << ".\n";
        Cleanup();
        gridfile.close();
        return false;
      }
      for (int k = 0; k < faces[j].type; ++k) {
        gridfile >> faces[j].edge[k];
      }
    }
    iLine += nFaces - 1;
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Faces' in file\n"
              << "    " << gridfilename << ".\n";
    Cleanup();
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << gridfilename << " (line " << iLine
              << ").\n";
    Cleanup();
    gridfile.close();
    return false;
  }

  // Get the elements.
  size_t nElements = 0;
  while (std::getline(gridfile, line)) {
    ++iLine;
    ltrim(line);
    // Find section 'Elements'.
    if (line.substr(0, 8) != "Elements") continue;
    // Get number of elements (given in brackets).
    const auto bra = line.find('(');
    const auto ket = line.find(')');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      // No closed brackets ().
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of elements.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    std::istringstream data;
    data.str(line);
    data >> nElements;
    data.clear();
    // Resize array of elements.
    m_elements.resize(nElements);
    // Get type and constituting edges of each element.
    for (size_t j = 0; j < nElements; ++j) {
      ++iLine;
      int type = 0;
      gridfile >> type;
      if (type == 2) {
        // Triangle
        int edge0, edge1, edge2;
        gridfile >> edge0 >> edge1 >> edge2;
        // Get the vertices.
        // Negative edge index means that the sequence of the two points
        // is supposed to be inverted.
        // The actual index is then given by "-index - 1".
        // For our purposes, the orientation does not matter.
        if (edge0 < 0) edge0 = -edge0 - 1;
        if (edge1 < 0) edge1 = -edge1 - 1;
        if (edge2 < 0) edge2 = -edge2 - 1;
        // Make sure the indices are not out of range.
        if (edge0 >= (int)nEdges || edge1 >= (int)nEdges || 
            edge2 >= (int)nEdges) {
          std::cerr << m_className << "::LoadGrid:\n    Error reading file "
                    << gridfilename << " (line " << iLine << ").\n"
                    << "    Edge index out of range.\n";
          Cleanup();
          gridfile.close();
          return false;
        }
        m_elements[j].vertex[0] = edgeP1[edge0];
        m_elements[j].vertex[1] = edgeP2[edge0];
        if (edgeP1[edge1] != m_elements[j].vertex[0] &&
            edgeP1[edge1] != m_elements[j].vertex[1]) {
          m_elements[j].vertex[2] = edgeP1[edge1];
        } else {
          m_elements[j].vertex[2] = edgeP2[edge1];
        }
      } else if (type == 5) {
        // Tetrahedron
        // Get the faces.
        // Negative face index means that the sequence of the edges
        // is supposed to be inverted.
        // For our purposes, the orientation does not matter.
        int face0, face1, face2, face3;
        gridfile >> face0 >> face1 >> face2 >> face3;
        if (face0 < 0) face0 = -face0 - 1;
        if (face1 < 0) face1 = -face1 - 1;
        if (face2 < 0) face2 = -face2 - 1;
        if (face3 < 0) face3 = -face3 - 1;
        // Make sure the face indices are not out of range.
        if (face0 >= (int)nFaces || face1 >= (int)nFaces || 
            face2 >= (int)nFaces || face3 >= (int)nFaces) {
          std::cerr << m_className << "::LoadGrid:\n    Error reading file "
                    << gridfilename << " (line " << iLine << ").\n"
                    << "    Face index out of range.\n";
          Cleanup();
          gridfile.close();
          return false;
        }
        // Get the edges of the first face.
        int edge0 = faces[face0].edge[0];
        int edge1 = faces[face0].edge[1];
        int edge2 = faces[face0].edge[2];
        if (edge0 < 0) edge0 = -edge0 - 1;
        if (edge1 < 0) edge1 = -edge1 - 1;
        if (edge2 < 0) edge2 = -edge2 - 1;
        // Make sure the edge indices are not out of range.
        if (edge0 >= (int)nEdges || edge1 >= (int)nEdges || 
            edge2 >= (int)nEdges) {
          std::cerr << m_className << "::LoadGrid:\n"
                    << "    Error reading file " << gridfilename << "\n"
                    << "    Edge index in element " << j << " out of range.\n";
          Cleanup();
          gridfile.close();
          return false;
        }
        // Get the first three vertices.
        m_elements[j].vertex[0] = edgeP1[edge0];
        m_elements[j].vertex[1] = edgeP2[edge0];
        if (edgeP1[edge1] != m_elements[j].vertex[0] &&
            edgeP1[edge1] != m_elements[j].vertex[1]) {
          m_elements[j].vertex[2] = edgeP1[edge1];
        } else {
          m_elements[j].vertex[2] = edgeP2[edge1];
        }
        // Get the fourth vertex from face 1.
        edge0 = faces[face1].edge[0];
        edge1 = faces[face1].edge[1];
        edge2 = faces[face1].edge[2];
        if (edge0 < 0) edge0 = -edge0 - 1;
        if (edge1 < 0) edge1 = -edge1 - 1;
        if (edge2 < 0) edge2 = -edge2 - 1;
        const auto v0 = m_elements[j].vertex[0];
        const auto v1 = m_elements[j].vertex[1];
        const auto v2 = m_elements[j].vertex[2];
        if (edgeP1[edge0] != v0 && edgeP1[edge0] != v1 && edgeP1[edge0] != v2) {
          m_elements[j].vertex[3] = edgeP1[edge0];
        } else if (edgeP2[edge0] != v0 && edgeP2[edge0] != v1 && 
                   edgeP2[edge0] != v2) {
          m_elements[j].vertex[3] = edgeP2[edge0];
        } else if (edgeP1[edge1] != v0 &&
                   edgeP1[edge1] != v1 &&
                   edgeP1[edge1] != v2) {
          m_elements[j].vertex[3] = edgeP1[edge1];
        } else if (edgeP2[edge1] != v0 &&
                   edgeP2[edge1] != v1 &&
                   edgeP2[edge1] != v2) {
          m_elements[j].vertex[3] = edgeP2[edge1];
        } else {
          std::cerr << m_className << "::LoadGrid:\n"
                    << "    Error reading file " << gridfilename << "\n"
                    << "    Face 1 of element " << j << " is degenerate.\n";
          Cleanup();
          gridfile.close();
          return false;
        }
      } else {
        // Other element types are not allowed.
        std::cerr << m_className << "::LoadGrid:\n"
                  << "    Error reading file " << gridfilename << " (line "
                  << iLine << ").\n";
        if (type == 0 || type == 1) {
          std::cerr << "    Invalid element type (" << type
                    << ") for 3d mesh.\n";
        } else {
          std::cerr << "    Element type " << type << " is not supported.\n"
                    << "    Remesh with option -t to create only"
                    << " triangles and tetrahedra.\n";
        }
        Cleanup();
        gridfile.close();
        return false;
      }
      m_elements[j].type = type;
      m_elements[j].region = -1;
    }
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Elements' in file\n"
              << "    " << gridfilename << ".\n";
    Cleanup();
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n    Error reading file "
              << gridfilename << " (line " << iLine << ").\n";
    Cleanup();
    gridfile.close();
    return false;
  }

  // Assign regions to elements.
  while (std::getline(gridfile, line)) {
    ltrim(line);
    // Find section 'Region'.
    if (line.substr(0, 6) != "Region") continue;
    // Get region name (given in brackets).
    auto bra = line.find('(');
    auto ket = line.find(')');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read region name.\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    std::istringstream data;
    data.str(line);
    std::string name;
    data >> name;
    data.clear();
    const size_t index = FindRegion(name);
    if (index >= m_regions.size()) {
      // Specified region name is not in the list.
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Error reading file " << gridfilename << ".\n"
                << "    Unknown region " << name << ".\n";
      continue;
    }
    std::getline(gridfile, line);
    std::getline(gridfile, line);
    bra = line.find('(');
    ket = line.find(')');
    if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
      // No closed brackets ().
      std::cerr << m_className << "::LoadGrid:\n    Error reading file "
                << gridfilename << ".\n    Could not read number of elements "
                << "in region " << name << ".\n";
      Cleanup();
      gridfile.close();
      return false;
    }
    line = line.substr(bra + 1, ket - bra - 1);
    int nElementsRegion;
    size_t iElement;
    data.str(line);
    data >> nElementsRegion;
    data.clear();
    for (int j = 0; j < nElementsRegion; ++j) {
      gridfile >> iElement;
      m_elements[iElement].region = index;
    }
  }

  gridfile.close();
  if (gridfile.fail() && !gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << gridfilename << ".\n";
    Cleanup();
    return false;
  }

  return true;
}

void ComponentTcad3d::Cleanup() {
  // Vertices
  m_vertices.clear();
  // Elements
  m_elements.clear();
  // Regions
  m_regions.clear();
  // Weighting fields and potentials.
  m_wf.clear();
  m_wp.clear();
}

bool ComponentTcad3d::CheckTetrahedron(
    const double x, const double y, const double z, const Element& element,
    std::array<double, nMaxVertices>& w) const {
  const Vertex& v0 = m_vertices[element.vertex[0]];
  const Vertex& v1 = m_vertices[element.vertex[1]];
  const Vertex& v2 = m_vertices[element.vertex[2]];
  const Vertex& v3 = m_vertices[element.vertex[3]];
  const double x10 = v1.x - v0.x;
  const double y10 = v1.y - v0.y;
  const double z10 = v1.z - v0.z;

  const double x20 = v2.x - v0.x;
  const double y20 = v2.y - v0.y;
  const double z20 = v2.z - v0.z;

  const double x30 = v3.x - v0.x;
  const double y30 = v3.y - v0.y;
  const double z30 = v3.z - v0.z;

  const double x21 = v2.x - v1.x;
  const double y21 = v2.y - v1.y;
  const double z21 = v2.z - v1.z;

  const double x31 = v3.x - v1.x;
  const double y31 = v3.y - v1.y;
  const double z31 = v3.z - v1.z;

  const double x32 = v3.x - v2.x;
  const double y32 = v3.y - v2.y;
  const double z32 = v3.z - v2.z;

  w[0] = (x - v1.x) * (y21 * z31 - y31 * z21) +
         (y - v1.y) * (z21 * x31 - z31 * x21) +
         (z - v1.z) * (x21 * y31 - x31 * y21);

  w[0] /= x10 * (y31 * z21 - y21 * z31) + y10 * (z31 * x21 - z21 * x31) +
          z10 * (x31 * y21 - x21 * y31);
  if (w[0] < 0.) return false;

  w[1] = (x - v2.x) * (-y20 * z32 + y32 * z20) +
         (y - v2.y) * (-z20 * x32 + z32 * x20) +
         (z - v2.z) * (-x20 * y32 + x32 * y20);

  w[1] /= x21 * (y20 * z32 - y32 * z20) + y21 * (z20 * x32 - z32 * x20) +
          z21 * (x20 * y32 - x32 * y20);
  if (w[1] < 0.) return false;

  w[2] = (x - v3.x) * (y30 * z31 - y31 * z30) +
         (y - v3.y) * (z30 * x31 - z31 * x30) +
         (z - v3.z) * (x30 * y31 - x31 * y30);

  w[2] /= x32 * (y31 * z30 - y30 * z31) + y32 * (z31 * x30 - z30 * x31) +
          z32 * (x31 * y30 - x30 * y31);
  if (w[2] < 0.) return false;

  w[3] = (x - v0.x) * (y20 * z10 - y10 * z20) +
         (y - v0.y) * (z20 * x10 - z10 * x20) +
         (z - v0.z) * (x20 * y10 - x10 * y20);

  w[3] /= x30 * (y20 * z10 - y10 * z20) + y30 * (z20 * x10 - z10 * x20) +
          z30 * (x20 * y10 - x10 * y20);
  if (w[3] < 0.) return false;

  if (m_debug) {
    // Reconstruct the point from the local coordinates.
    const double xr = w[0] * v0.x + w[1] * v1.x + w[2] * v2.x + w[3] * v3.x;
    const double yr = w[0] * v0.y + w[1] * v1.y + w[2] * v2.y + w[3] * v3.y;
    const double zr = w[0] * v0.z + w[1] * v1.z + w[2] * v2.z + w[3] * v3.z;
    std::cout << m_className << "::CheckTetrahedron:\n"
              << "    Original coordinates:      (" << x << ", " << y << ", "
              << z << ")\n"
              << "    Local coordinates:         (" << w[0] << ", " << w[1]
              << ", " << w[2] << ", " << w[3] << ")\n"
              << "    Reconstructed coordinates: (" << xr << ", " << yr << ", "
              << zr << ")\n"
              << "    Checksum: " << w[0] + w[1] + w[2] + w[3] - 1. << "\n";
  }

  return true;
}

bool ComponentTcad3d::CheckTriangle(
    const double x, const double y, const double z, const Element& element,
    std::array<double, nMaxVertices>& w) const {
  const Vertex& v0 = m_vertices[element.vertex[0]];
  const Vertex& v1 = m_vertices[element.vertex[1]];
  const Vertex& v2 = m_vertices[element.vertex[2]];

  const double v1x = v1.x - v0.x;
  const double v2x = v2.x - v0.x;
  const double v1y = v1.y - v0.y;
  const double v2y = v2.y - v0.y;
  const double v1z = v1.z - v0.z;
  const double v2z = v2.z - v0.z;

  // Check whether the point lies in the plane of the triangle.
  // Compute the coefficients of the plane equation.
  const double a = v1y * v2z - v2y * v1z;
  const double b = v1z * v2x - v2z * v1x;
  const double c = v1x * v2y - v2x * v1y;
  const double d = a * v0.x + b * v0.y + c * v0.z;
  // Check if the point satisfies the plane equation.
  if (a * x + b * y + c * z != d) return false;

  // Map (x, y) onto local variables (b, c) such that
  // P = A + b * (B - A) + c * (C - A)
  // A point P is inside the triangle ABC if b, c > 0 and b + c < 1;
  // b, c are also weighting factors for points B, C
  w[1] = ((x - v0.x) * v2y - (y - v0.y) * v2x) / (v1x * v2y - v1y * v2x);
  if (w[1] < 0. || w[1] > 1.) return false;
  w[2] = ((v0.x - x) * v1y - (v0.y - y) * v1x) / (v1x * v2y - v1y * v2x);
  if (w[2] < 0. || w[1] + w[2] > 1.) return false;

  // Weighting factor for point A
  w[0] = 1. - w[1] - w[2];

  return true;
}

void ComponentTcad3d::Reset() {
  Cleanup();
  m_ready = false;
}

void ComponentTcad3d::UpdatePeriodicity() {
  if (!m_ready) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Field map not available.\n";
    return;
  }

  for (size_t i = 0; i < 3; ++i) {
    // Check for conflicts.
    if (m_periodic[i] && m_mirrorPeriodic[i]) {
      std::cerr << m_className << "::UpdatePeriodicity:\n"
                << "    Both simple and mirror periodicity requested. Reset.\n";
      m_periodic[i] = m_mirrorPeriodic[i] = false;
    }
  }

  if (m_axiallyPeriodic[0] || m_axiallyPeriodic[1] || m_axiallyPeriodic[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Axial symmetry is not supported. Reset.\n";
    m_axiallyPeriodic.fill(false);
  }

  if (m_rotationSymmetric[0] || m_rotationSymmetric[1] ||
      m_rotationSymmetric[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Rotation symmetry is not supported. Reset.\n";
    m_rotationSymmetric.fill(false);
  }
}

void ComponentTcad3d::MapCoordinates(double& x, double& y, double& z,
                                     bool& xmirr, bool& ymirr,
                                     bool& zmirr) const {
  xmirr = false;
  const double cellsx = m_xMaxBB - m_xMinBB;
  if (m_periodic[0]) {
    x = m_xMinBB + fmod(x - m_xMinBB, cellsx);
    if (x < m_xMinBB) x += cellsx;
  } else if (m_mirrorPeriodic[0]) {
    double xNew = m_xMinBB + fmod(x - m_xMinBB, cellsx);
    if (xNew < m_xMinBB) xNew += cellsx;
    const int nx = int(floor(0.5 + (xNew - x) / cellsx));
    if (nx != 2 * (nx / 2)) {
      xNew = m_xMinBB + m_xMaxBB - xNew;
      xmirr = true;
    }
    x = xNew;
  }
  ymirr = false;
  const double cellsy = m_yMaxBB - m_yMinBB;
  if (m_periodic[1]) {
    y = m_yMinBB + fmod(y - m_yMinBB, cellsy);
    if (y < m_yMinBB) y += cellsy;
  } else if (m_mirrorPeriodic[1]) {
    double yNew = m_yMinBB + fmod(y - m_yMinBB, cellsy);
    if (yNew < m_yMinBB) yNew += cellsy;
    const int ny = int(floor(0.5 + (yNew - y) / cellsy));
    if (ny != 2 * (ny / 2)) {
      yNew = m_yMinBB + m_yMaxBB - yNew;
      ymirr = true;
    }
    y = yNew;
  }
  zmirr = false;
  const double cellsz = m_zMaxBB - m_zMinBB;
  if (m_periodic[2]) {
    z = m_zMinBB + fmod(z - m_zMinBB, cellsz);
    if (z < m_zMinBB) z += cellsz;
  } else if (m_mirrorPeriodic[2]) {
    double zNew = m_zMinBB + fmod(z - m_zMinBB, cellsz);
    if (zNew < m_zMinBB) zNew += cellsz;
    const int nz = int(floor(0.5 + (zNew - z) / cellsz));
    if (nz != 2 * (nz / 2)) {
      zNew = m_zMinBB + m_zMaxBB - zNew;
      zmirr = true;
    }
    z = zNew;
  }
}

size_t ComponentTcad3d::FindRegion(const std::string& name) const {
  const size_t nRegions = m_regions.size();
  for (size_t j = 0; j < nRegions; ++j) {
    if (name == m_regions[j].name) return j;
  }
  return nRegions;
}
}
