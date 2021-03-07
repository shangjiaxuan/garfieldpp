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

namespace {

bool ExtractFromSquareBrackets(std::string& line) {

  const auto bra = line.find('[');
  const auto ket = line.find(']');
  if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
    return false;
  }
  line = line.substr(bra + 1, ket - bra - 1);
  return true;
}

bool ExtractFromBrackets(std::string& line) {

  const auto bra = line.find('(');
  const auto ket = line.find(')');
  if (ket < bra || bra == std::string::npos || ket == std::string::npos) {
    return false;
  }
  line = line.substr(bra + 1, ket - bra - 1);
  return true;
}

}

namespace Garfield {

ComponentTcad3d::ComponentTcad3d() : ComponentTcadBase("Tcad3d") {}

void ComponentTcad3d::ElectricField(const double xin, const double yin,
                                    const double zin, double& ex, double& ey,
                                    double& ez, double& p, Medium*& m,
                                    int& status) {
  // Assume this will work.
  status = 0;
  ex = ey = ez = p = 0.;
  m = nullptr;
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::ElectricField:\n"
              << "    Field map is not available for interpolation.\n";
    status = -10;
    return;
  }
  std::array<double, 3> x = {xin, yin, zin};
  std::array<bool, 3> mirr = {false, false, false};
  // In case of periodicity, reduce to the cell volume.
  MapCoordinates(x, mirr);
  // Check if the point is inside the bounding box.
  if (!InBoundingBox(x)) {
    status = -6;
    return;
  }

  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x[0], x[1], x[2], w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    status = -6;
    return;
  }
  const Element& element = m_elements[i];
  const unsigned int nVertices = element.type == 2 ? 3 : 4;
  for (unsigned int j = 0; j < nVertices; ++j) {
    const auto index = element.vertex[j];
    ex += w[j] * m_efield[index][0];
    ey += w[j] * m_efield[index][1];
    ez += w[j] * m_efield[index][2];
    p += w[j] * m_potential[index];
  }
  if (mirr[0]) ex = -ex;
  if (mirr[1]) ey = -ey;
  if (mirr[2]) ez = -ez;
  m = m_regions[element.region].medium;
  if (!m_regions[element.region].drift || !m) status = -5;
}

void ComponentTcad3d::ElectricField(const double x, const double y,
                                    const double z, double& ex, double& ey,
                                    double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

bool ComponentTcad3d::Interpolate(
    const double xin, const double yin, const double zin,
    const std::vector<std::array<double, 3> >& field,
    double& fx, double& fy, double& fz) {

  if (field.empty()) return false;
  std::array<double, 3> x = {xin, yin, zin};
  std::array<bool, 3> mirr = {false, false, false};
  // In case of periodicity, reduce to the cell volume.
  MapCoordinates(x, mirr);
  // Make sure the point is inside the bounding box.
  if (!InBoundingBox(x)) return false;

  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x[0], x[1], x[2], w);
  // Stop if the point is outside the mesh.
  if (i >= m_elements.size()) return false;

  const Element& element = m_elements[i];
  const size_t nVertices = element.type == 2 ? 3 : 4;
  for (size_t j = 0; j < nVertices; ++j) {
    const auto index = element.vertex[j];
    fx += w[j] * field[index][0];
    fy += w[j] * field[index][1];
    fz += w[j] * field[index][2];
  }
  if (mirr[0]) fx = -fx;
  if (mirr[1]) fy = -fy;
  if (mirr[2]) fz = -fz;
  return true;
} 

bool ComponentTcad3d::Interpolate(
    const double xin, const double yin, const double zin,
    const std::vector<double>& field, double& f) {

  f = 0.;
  if (field.empty()) return false;
  std::array<double, 3> x = {xin, yin, zin};
  std::array<bool, 3> mirr = {false, false, false};
  // In case of periodicity, reduce to the cell volume.
  MapCoordinates(x, mirr);
  // Make sure the point is inside the bounding box.
  if (!InBoundingBox(x)) return false;

  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x[0], x[1], x[2], w);
  // Stop if the point is outside the mesh.
  if (i >= m_elements.size()) return false;

  const Element& element = m_elements[i];
  const size_t nVertices = element.type == 2 ? 3 : 4;
  for (size_t j = 0; j < nVertices; ++j) {
    f += w[j] * field[element.vertex[j]];
  }
  return true;
}

Medium* ComponentTcad3d::GetMedium(const double xin, const double yin,
                                   const double zin) {
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::GetMedium:\n"
              << "    Field map not available for interpolation.\n";
    return nullptr;
  }

  std::array<double, 3> x = {xin, yin, zin};
  std::array<bool, 3> mirr = {false, false, false};
  MapCoordinates(x, mirr);
  // Check if the point is inside the bounding box.
  if (!InBoundingBox(x)) return nullptr;

  // Determine the shape functions.
  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x[0], x[1], x[2], w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    return nullptr;
  }
  return m_regions[m_elements[i].region].medium;
}

bool ComponentTcad3d::Initialise(const std::string& gridfilename,
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

  // Import electric field and potential from .dat file.
  if (!LoadData(datafilename)) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Importing electric field and potential failed.\n";
    Cleanup();
    return false;
  }

  // Find min./max. coordinates and potentials.
  m_bbMax[0] = m_vertices[m_elements[0].vertex[0]][0];
  m_bbMax[1] = m_vertices[m_elements[0].vertex[0]][1];
  m_bbMax[2] = m_vertices[m_elements[0].vertex[0]][2];
  m_bbMin = m_bbMax;
  const size_t nElements = m_elements.size();
  for (size_t i = 0; i < nElements; ++i) {
    Element& element = m_elements[i];
    double xmin = m_vertices[element.vertex[0]][0];
    double ymin = m_vertices[element.vertex[0]][1];
    double zmin = m_vertices[element.vertex[0]][2];
    double xmax = xmin;
    double ymax = ymin;
    double zmax = zmin;
    const size_t nVertices = element.type == 2 ? 3 : 4;
    for (size_t j = 0; j < nVertices; ++j) {
      const auto& vj = m_vertices[element.vertex[j]];
      if (vj[0] < xmin) xmin = vj[0];
      if (vj[0] > xmax) xmax = vj[0];
      if (vj[1] < ymin) ymin = vj[1];
      if (vj[1] > ymax) ymax = vj[1];
      if (vj[2] < zmin) zmin = vj[2];
      if (vj[2] > zmax) zmax = vj[2];
    }
    constexpr double tol = 1.e-6;
    element.bbMin[0] = xmin - tol;
    element.bbMax[0] = xmax + tol;
    element.bbMin[1] = ymin - tol;
    element.bbMax[1] = ymax + tol;
    element.bbMin[2] = zmin - tol;
    element.bbMax[2] = zmax + tol;
    m_bbMin[0] = std::min(m_bbMin[0], xmin);
    m_bbMax[0] = std::max(m_bbMax[0], xmax);
    m_bbMin[1] = std::min(m_bbMin[1], ymin);
    m_bbMax[1] = std::max(m_bbMax[1], ymax);
    m_bbMin[2] = std::min(m_bbMin[2], zmin);
    m_bbMax[2] = std::max(m_bbMax[2], zmax);
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
  std::cout << m_className << "::Initialise:\n"
            << "    Bounding box:\n"
            << "      " << m_bbMin[0] << " < x [cm] < " << m_bbMax[0] << "\n"
            << "      " << m_bbMin[1] << " < y [cm] < " << m_bbMax[1] << "\n"
            << "      " << m_bbMin[2] << " < z [cm] < " << m_bbMax[2] << "\n"
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
    if (element.region < nRegions) {
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
  if (!ok) {
    m_ready = false;
    Cleanup();
    return false;
  }

  // Set up the octree.
  const float hx = 0.5 * (m_bbMax[0] - m_bbMin[0]);
  const float hy = 0.5 * (m_bbMax[1] - m_bbMin[1]);
  const float hz = 0.5 * (m_bbMax[2] - m_bbMin[2]);
  m_tree.reset(new TetrahedralTree(Vec3(m_bbMin[0] + hx, m_bbMin[1] + hy, m_bbMin[2] + hz),
                                   Vec3(hx, hy, hz)));

  // Insert the mesh nodes in the tree.
  for (size_t i = 0; i < nVertices; ++i) {
    const auto& vtx = m_vertices[i];
    m_tree->InsertMeshNode(Vec3(vtx[0], vtx[1], vtx[2]), i);
  }

  // Insert the mesh elements in the tree.
  for (size_t i = 0; i < nElements; ++i) {
    const Element& e = m_elements[i];
    const double bb[6] = {e.bbMin[0], e.bbMin[1], e.bbMin[2], 
                          e.bbMax[0], e.bbMax[1], e.bbMax[2]};
    m_tree->InsertMeshElement(bb, i);
  }

  m_ready = true;
  UpdatePeriodicity();
  std::cout << m_className << "::Initialise: Initialisation finished.\n";
  return true;
}

bool ComponentTcad3d::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                     double& xmax, double& ymax, double& zmax) {
  if (!m_ready) return false;
  xmin = m_bbMin[0];
  ymin = m_bbMin[1];
  zmin = m_bbMin[2];
  xmax = m_bbMax[0];
  ymax = m_bbMax[1];
  zmax = m_bbMax[2];
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

size_t ComponentTcad3d::FindElement(
    const double x, const double y, const double z,
    std::array<double, nMaxVertices>& w) const {

  w.fill(0.);

  std::vector<int> elementsToSearch;
  if (m_tree) {
    elementsToSearch = m_tree->GetElementsInBlock(Vec3(x, y, z));
  }
  const size_t nElementsToSearch = m_tree ? elementsToSearch.size() : m_elements.size(); 
  // Loop over the elements.
  for (size_t i = 0; i < nElementsToSearch; ++i) {
    const size_t idx = m_tree ? elementsToSearch[i] : i;
    if (InElement(x, y, z, m_elements[idx], w)) return idx;
  }

  if (m_debug) {
    std::cerr << m_className << "::FindElement:\n"
              << "    Point (" << x << ", " << y << ", " << z
              << ") is outside the mesh.\n";
  }
  return m_elements.size();
}

bool ComponentTcad3d::GetElement(const size_t i, double& vol,
                                 double& dmin, double& dmax, int& type,
                                 std::vector<size_t>& nodes, int& reg) const {
  nodes.clear();
  if (i >= m_elements.size()) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }

  const Element& element = m_elements[i];
  if (element.type == 2) {
    // Triangle
    const auto& v0 = m_vertices[element.vertex[0]];
    const auto& v1 = m_vertices[element.vertex[1]];
    const auto& v2 = m_vertices[element.vertex[2]];
    const double vx = (v1[1] - v0[1]) * (v2[2] - v0[2]) - 
                      (v1[2] - v0[2]) * (v2[1] - v0[1]);
    const double vy = (v1[2] - v0[2]) * (v2[0] - v0[0]) - 
                      (v1[0] - v0[0]) * (v2[2] - v0[2]);
    const double vz = (v1[0] - v0[0]) * (v2[1] - v0[1]) - 
                      (v1[1] - v0[1]) * (v2[0] - v0[0]);
    vol = sqrt(vx * vx + vy * vy + vz * vz);
    const double a = sqrt(pow(v1[0] - v0[0], 2) + pow(v1[1] - v0[1], 2) + 
                          pow(v1[2] - v0[2], 2));
    const double b = sqrt(pow(v2[0] - v0[0], 2) + pow(v2[1] - v0[1], 2) + 
                          pow(v2[2] - v0[2], 2));
    const double c = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2) + 
                          pow(v1[2] - v2[2], 2));
    dmin = std::min({a, b, c});
    dmax = std::max({a, b, c});
  } else if (element.type == 5) {
    // Tetrahedron
    const auto& v0 = m_vertices[element.vertex[0]];
    const auto& v1 = m_vertices[element.vertex[1]];
    const auto& v2 = m_vertices[element.vertex[2]];
    const auto& v3 = m_vertices[element.vertex[3]];
    vol = fabs((v3[0] - v0[0]) * ((v1[1] - v0[1]) * (v2[2] - v0[2]) -
                                  (v2[1] - v0[1]) * (v1[2] - v0[2])) +
               (v3[1] - v0[1]) * ((v1[2] - v0[2]) * (v2[0] - v0[0]) -
                                  (v2[2] - v0[2]) * (v1[0] - v0[0])) +
               (v3[2] - v0[2]) * ((v1[0] - v0[0]) * (v2[1] - v0[1]) -
                                  (v3[0] - v0[0]) * (v1[1] - v0[1]))) /
          6.;
    // Loop over all pairs of m_vertices.
    for (size_t j = 0; j < nMaxVertices - 1; ++j) {
      const auto& vj = m_vertices[element.vertex[j]];
      for (size_t k = j + 1; k < nMaxVertices; ++k) {
        const auto& vk = m_vertices[element.vertex[k]];
        // Compute distance.
        const double dist = sqrt(pow(vj[0] - vk[0], 2) + pow(vj[1] - vk[1], 2) +
                                 pow(vj[2] - vk[2], 2));
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

  x = m_vertices[i][0];
  y = m_vertices[i][1];
  z = m_vertices[i][2];
  if (!m_potential.empty()) v = m_potential[i];
  if (!m_efield.empty()) {
    ex = m_efield[i][0];
    ey = m_efield[i][1];
    ez = m_efield[i][2];
  }
  return true;
}

bool ComponentTcad3d::LoadGrid(const std::string& filename) {
  // Open the file containing the mesh description.
  std::ifstream gridfile;
  gridfile.open(filename.c_str(), std::ios::in);
  if (!gridfile) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not open file " << filename << ".\n";
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
              << "    " << filename << ".\n";
    return false;
  } else if (gridfile.fail()) {
    // Error reading from the file.
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
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
    if (!ExtractFromSquareBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read region names.\n";
      return false;
    }
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
              << "    " << filename << ".\n";
    return false;
  } else if (gridfile.fail()) {
    // Error reading from the file.
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
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
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of vertices.\n";
      return false;
    }
    std::istringstream data;
    data.str(line);
    data >> nVertices;
    m_vertices.resize(nVertices);
    // Get the coordinates of this vertex.
    for (size_t j = 0; j < nVertices; ++j) {
      gridfile >> m_vertices[j][0] >> m_vertices[j][1] >> m_vertices[j][2];
      // Change units from micron to cm.
      m_vertices[j][0] *= 1.e-4;
      m_vertices[j][1] *= 1.e-4;
      m_vertices[j][2] *= 1.e-4;
      ++iLine;
    }
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Vertices' in file\n"
              << "    " << filename << ".\n";
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
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
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of edges.\n";
      return false;
    }
    std::istringstream data;
    data.str(line);
    data >> nEdges;
    edgeP1.resize(nEdges);
    edgeP2.resize(nEdges);
    // Get the indices of the two endpoints.
    for (size_t j = 0; j < nEdges; ++j) {
      gridfile >> edgeP1[j] >> edgeP2[j];
      ++iLine;
    }
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Edges' in file\n"
              << "    " << filename << ".\n";
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
    return false;
  }

  for (size_t i = 0; i < nEdges; ++i) {
    // Make sure the indices of the edge endpoints are not out of range.
    if (edgeP1[i] >= nVertices || edgeP2[i] >= nVertices) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Vertex index of edge " << i << " out of range.\n";
      return false;
    }
    // Make sure the edge is non-degenerate.
    if (edgeP1[i] == edgeP2[i]) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Edge " << i << " is degenerate.\n";
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
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of faces.\n";
      return false;
    }
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
              << "    " << filename << ".\n";
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
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
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of elements.\n";
      return false;
    }
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
                    << filename << " (line " << iLine << ").\n"
                    << "    Edge index out of range.\n";
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
                    << filename << " (line " << iLine << ").\n"
                    << "    Face index out of range.\n";
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
                    << "    Error reading file " << filename << "\n"
                    << "    Edge index in element " << j << " out of range.\n";
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
                    << "    Error reading file " << filename << "\n"
                    << "    Face 1 of element " << j << " is degenerate.\n";
          return false;
        }
      } else {
        // Other element types are not allowed.
        std::cerr << m_className << "::LoadGrid:\n"
                  << "    Error reading file " << filename << " (line "
                  << iLine << ").\n";
        if (type == 0 || type == 1) {
          std::cerr << "    Invalid element type (" << type
                    << ") for 3d mesh.\n";
        } else {
          std::cerr << "    Element type " << type << " is not supported.\n"
                    << "    Remesh with option -t to create only"
                    << " triangles and tetrahedra.\n";
        }
        return false;
      }
      m_elements[j].type = type;
      m_elements[j].region = m_regions.size();
    }
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Elements' in file\n"
              << "    " << filename << ".\n";
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine 
              << ").\n";
    return false;
  }

  // Assign regions to elements.
  while (std::getline(gridfile, line)) {
    ltrim(line);
    // Find section 'Region'.
    if (line.substr(0, 6) != "Region") continue;
    // Get region name (given in brackets).
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read region name.\n";
      return false;
    }
    std::istringstream data;
    data.str(line);
    std::string name;
    data >> name;
    data.clear();
    const size_t index = FindRegion(name);
    if (index >= m_regions.size()) {
      // Specified region name is not in the list.
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Error reading file " << filename << ".\n"
                << "    Unknown region " << name << ".\n";
      continue;
    }
    std::getline(gridfile, line);
    std::getline(gridfile, line);
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Error reading file " << filename << ".\n"
                << "    Could not read number of elements in region " 
                << name << ".\n";
      return false;
    }
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

  if (gridfile.fail() && !gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << ".\n";
    return false;
  }
  return true;
}

bool ComponentTcad3d::InTetrahedron(
    const double x, const double y, const double z, const Element& element,
    std::array<double, nMaxVertices>& w) const {
  const auto& v0 = m_vertices[element.vertex[0]];
  const auto& v1 = m_vertices[element.vertex[1]];
  const auto& v2 = m_vertices[element.vertex[2]];
  const auto& v3 = m_vertices[element.vertex[3]];
  const double x10 = v1[0] - v0[0];
  const double y10 = v1[1] - v0[1];
  const double z10 = v1[2] - v0[2];

  const double x20 = v2[0] - v0[0];
  const double y20 = v2[1] - v0[1];
  const double z20 = v2[2] - v0[2];

  const double x30 = v3[0] - v0[0];
  const double y30 = v3[1] - v0[1];
  const double z30 = v3[2] - v0[2];

  const double x21 = v2[0] - v1[0];
  const double y21 = v2[1] - v1[1];
  const double z21 = v2[2] - v1[2];

  const double x31 = v3[0] - v1[0];
  const double y31 = v3[1] - v1[1];
  const double z31 = v3[2] - v1[2];

  const double x32 = v3[0] - v2[0];
  const double y32 = v3[1] - v2[1];
  const double z32 = v3[2] - v2[2];

  w[0] = (x - v1[0]) * (y21 * z31 - y31 * z21) +
         (y - v1[1]) * (z21 * x31 - z31 * x21) +
         (z - v1[2]) * (x21 * y31 - x31 * y21);

  w[0] /= x10 * (y31 * z21 - y21 * z31) + y10 * (z31 * x21 - z21 * x31) +
          z10 * (x31 * y21 - x21 * y31);
  if (w[0] < 0.) return false;

  w[1] = (x - v2[0]) * (-y20 * z32 + y32 * z20) +
         (y - v2[1]) * (-z20 * x32 + z32 * x20) +
         (z - v2[2]) * (-x20 * y32 + x32 * y20);

  w[1] /= x21 * (y20 * z32 - y32 * z20) + y21 * (z20 * x32 - z32 * x20) +
          z21 * (x20 * y32 - x32 * y20);
  if (w[1] < 0.) return false;

  w[2] = (x - v3[0]) * (y30 * z31 - y31 * z30) +
         (y - v3[1]) * (z30 * x31 - z31 * x30) +
         (z - v3[2]) * (x30 * y31 - x31 * y30);

  w[2] /= x32 * (y31 * z30 - y30 * z31) + y32 * (z31 * x30 - z30 * x31) +
          z32 * (x31 * y30 - x30 * y31);
  if (w[2] < 0.) return false;

  w[3] = (x - v0[0]) * (y20 * z10 - y10 * z20) +
         (y - v0[1]) * (z20 * x10 - z10 * x20) +
         (z - v0[2]) * (x20 * y10 - x10 * y20);

  w[3] /= x30 * (y20 * z10 - y10 * z20) + y30 * (z20 * x10 - z10 * x20) +
          z30 * (x20 * y10 - x10 * y20);
  if (w[3] < 0.) return false;

  if (m_debug) {
    // Reconstruct the point from the local coordinates.
    const double xr = w[0] * v0[0] + w[1] * v1[0] + w[2] * v2[0] + w[3] * v3[0];
    const double yr = w[0] * v0[1] + w[1] * v1[1] + w[2] * v2[1] + w[3] * v3[1];
    const double zr = w[0] * v0[2] + w[1] * v1[2] + w[2] * v2[2] + w[3] * v3[2];
    std::cout << m_className << "::InTetrahedron:\n"
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

bool ComponentTcad3d::InTriangle(
    const double x, const double y, const double z, const Element& element,
    std::array<double, nMaxVertices>& w) const {
  const auto& v0 = m_vertices[element.vertex[0]];
  const auto& v1 = m_vertices[element.vertex[1]];
  const auto& v2 = m_vertices[element.vertex[2]];

  const double v1x = v1[0] - v0[0];
  const double v2x = v2[0] - v0[0];
  const double v1y = v1[1] - v0[1];
  const double v2y = v2[1] - v0[1];
  const double v1z = v1[2] - v0[2];
  const double v2z = v2[2] - v0[2];

  // Check whether the point lies in the plane of the triangle.
  // Compute the coefficients of the plane equation.
  const double a = v1y * v2z - v2y * v1z;
  const double b = v1z * v2x - v2z * v1x;
  const double c = v1x * v2y - v2x * v1y;
  const double d = a * v0[0] + b * v0[1] + c * v0[2];
  // Check if the point satisfies the plane equation.
  if (a * x + b * y + c * z != d) return false;

  // Map (x, y) onto local variables (b, c) such that
  // P = A + b * (B - A) + c * (C - A)
  // A point P is inside the triangle ABC if b, c > 0 and b + c < 1;
  // b, c are also weighting factors for points B, C
  w[1] = ((x - v0[0]) * v2y - (y - v0[1]) * v2x) / (v1x * v2y - v1y * v2x);
  if (w[1] < 0. || w[1] > 1.) return false;
  w[2] = ((v0[0] - x) * v1y - (v0[1] - y) * v1x) / (v1x * v2y - v1y * v2x);
  if (w[2] < 0. || w[1] + w[2] > 1.) return false;

  // Weighting factor for point A
  w[0] = 1. - w[1] - w[2];

  return true;
}

void ComponentTcad3d::Reset() {
  Cleanup();
  m_ready = false;
}

}
