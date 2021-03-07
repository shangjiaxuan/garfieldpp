#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Garfield/ComponentTcad2d.hh"
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

bool ComponentTcad2d::LoadData(const std::string& filename) {
  std::ifstream datafile;
  datafile.open(filename.c_str(), std::ios::in);
  if (!datafile) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }

  const size_t nVertices = m_vertices.size();
  std::vector<unsigned int> fillCount(nVertices, 0);

  // Read the file line by line.
  std::string line;
  while (std::getline(datafile, line)) {
    // Strip white space from the beginning of the line.
    ltrim(line);
    // Find data section.
    if (line.substr(0, 8) != "function") continue;
    // Read type of data set.
    const std::string::size_type pEq = line.find('=');
    if (pEq == std::string::npos) {
      // No "=" found.
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading file " << filename << ".\n"
                << "    Line:\n    " << line << "\n";
      datafile.close();
      return false;
    }
    line = line.substr(pEq + 1);
    std::string dataset;
    std::istringstream data;
    data.str(line);
    data >> dataset;
    if (dataset == "ElectrostaticPotential") {
      m_potential.assign(nVertices, 0.);
      if (!ReadDataset(datafile, dataset)) {
        m_potential.clear();
        return false;
      }
    } else if (dataset == "ElectricField") {
      m_efield.assign(nVertices, {0., 0.});
      if (!ReadDataset(datafile, dataset)) {
        m_efield.clear();
        return false;
      }
    } else if (dataset == "eDriftVelocity") {
      m_eVelocity.assign(nVertices, {0., 0.});
      if (!ReadDataset(datafile, dataset)) {
        m_eVelocity.clear();
        return false;
      }
    } else if (dataset == "hDriftVelocity") {
      m_hVelocity.assign(nVertices, {0., 0.});
      if (!ReadDataset(datafile, dataset)) {
        m_hVelocity.clear();
        return false;
      }
    } else if (dataset == "eMobility") {
      m_eMobility.assign(nVertices, 0.);
      if (!ReadDataset(datafile, dataset)) {
        m_eMobility.clear();
        return false;
      }
    } else if (dataset == "hMobility") {
      m_hMobility.assign(nVertices, 0.);
      if (!ReadDataset(datafile, dataset)) {
        m_hMobility.clear();
        return false;
      }
    } else if (dataset == "eLifetime") {
      m_eLifetime.assign(nVertices, 0.);
      if (!ReadDataset(datafile, dataset)) {
        m_eLifetime.clear();
        return false;
      }
    } else if (dataset == "hLifetime") {
      m_hLifetime.assign(nVertices, 0.);
      if (!ReadDataset(datafile, dataset)) {
        m_hLifetime.clear();
        return false;
      }
    } else if (dataset.substr(0, 14) == "TrapOccupation" &&
               dataset.substr(17, 2) == "Do") {
      if (!ReadDataset(datafile, dataset)) return false;
      Defect donor;
      donor.xsece = -1.;
      donor.xsech = -1.;
      donor.conc = -1.;
      m_donors.push_back(donor);
    } else if (dataset.substr(0, 14) == "TrapOccupation" &&
               dataset.substr(17, 2) == "Ac") {
      if (!ReadDataset(datafile, dataset)) return false;
      Defect acceptor;
      acceptor.xsece = -1.;
      acceptor.xsech = -1.;
      acceptor.conc = -1.;
      m_acceptors.push_back(acceptor);
    }
  }
  if (datafile.fail() && !datafile.eof()) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Error reading file " << filename << "\n";
    datafile.close();
    return false;
  }
  datafile.close();
  return true;
}

bool ComponentTcad2d::ReadDataset(std::ifstream& datafile,
                                  const std::string& dataset) {
  enum DataSet {
    ElectrostaticPotential,
    EField,
    eDriftVelocity,
    hDriftVelocity,
    eMobility,
    hMobility,
    eLifetime,
    hLifetime,
    DonorTrapOccupation,
    AcceptorTrapOccupation,
    Unknown
  };
  DataSet ds = Unknown;
  if (dataset == "ElectrostaticPotential") {
    ds = ElectrostaticPotential;
  } else if (dataset == "ElectricField") {
    ds = EField;
  } else if (dataset == "eDriftVelocity") {
    ds = eDriftVelocity;
  } else if (dataset == "hDriftVelocity") {
    ds = hDriftVelocity;
  } else if (dataset == "eMobility") {
    ds = eMobility;
  } else if (dataset == "hMobility") {
    ds = hMobility;
  } else if (dataset == "eLifetime") {
    ds = eLifetime;
  } else if (dataset == "hLifetime") {
    ds = hLifetime;
  } else if (dataset.substr(0, 14) == "TrapOccupation") {
    if (dataset.substr(17, 2) == "Do") {
      ds = DonorTrapOccupation;
    } else if (dataset.substr(17, 2) == "Ac") {
      ds = AcceptorTrapOccupation;
    }
  } else {
    std::cerr << m_className << "::ReadDataset:\n"
              << "    Unexpected dataset " << dataset << ".\n";
    return false;
  }

  bool isVector = false;
  if (ds == EField || ds == eDriftVelocity || ds == hDriftVelocity) {
    isVector = true;
  }

  if (!datafile.is_open()) return false;
  std::string line;
  std::getline(datafile, line);
  std::getline(datafile, line);
  std::getline(datafile, line);
  std::getline(datafile, line);
  // Get the region name (given in brackets).
  if (!ExtractFromSquareBrackets(line)) {
    std::cerr << m_className << "::ReadDataset:\n"
              << "    Cannot extract region name.\n"
              << "    Line:\n    " << line << "\n";
    datafile.close();
    return false;
  }
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
    return false;
  }
  // Get the number of values.
  std::getline(datafile, line);
  if (!ExtractFromBrackets(line)) {
    std::cerr << m_className << "::ReadDataset:\n"
              << "    Cannot extract number of values to be read.\n"
              << "    Line:\n    " << line << "\n";
    datafile.close();
    return false;
  }
  int nValues;
  data.str(line);
  data >> nValues;
  data.clear();
  if (isVector) nValues = nValues / 2;
  // Mark the vertices belonging to this region.
  const size_t nVertices = m_vertices.size();
  std::vector<bool> isInRegion(nVertices, false);
  const size_t nElements = m_elements.size();
  for (size_t j = 0; j < nElements; ++j) {
    if (m_elements[j].region != index) continue;
    for (int k = 0; k <= m_elements[j].type; ++k) {
      isInRegion[m_elements[j].vertex[k]] = true;
    }
  }

  unsigned int ivertex = 0;
  for (int j = 0; j < nValues; ++j) {
    // Read the next value.
    double val1 = 0., val2 = 0.;
    if (isVector) {
      datafile >> val1 >> val2;
    } else {
      datafile >> val1;
    }
    // Find the next vertex belonging to the region.
    while (ivertex < nVertices) {
      if (isInRegion[ivertex]) break;
      ++ivertex;
    }
    // Check if there is a mismatch between the number of vertices
    // and the number of values.
    if (ivertex >= nVertices) {
      std::cerr << m_className << "::ReadDataset:\n"
                << "    Dataset " << dataset
                << " has more values than vertices in region " << name << "\n";
      datafile.close();
      return false;
    }
    switch (ds) {
      case ElectrostaticPotential:
        m_potential[ivertex] = val1;
        break;
      case EField:
        m_efield[ivertex][0] = val1;
        m_efield[ivertex][1] = val2;
        break;
      case eDriftVelocity:
        // Scale from cm/s to cm/ns.
        m_eVelocity[ivertex][0] = val1 * 1.e-9;
        m_eVelocity[ivertex][1] = val2 * 1.e-9;
        break;
      case hDriftVelocity:
        // Scale from cm/s to cm/ns.
        m_hVelocity[ivertex][0] = val1 * 1.e-9;
        m_hVelocity[ivertex][1] = val2 * 1.e-9;
        break;
      case eMobility:
        // Convert from cm2 / (V s) to cm2 / (V ns).
        m_eMobility[ivertex] = val1 * 1.e-9;
        break;
      case hMobility:
        // Convert from cm2 / (V s) to cm2 / (V ns).
        m_hMobility[ivertex] = val1 * 1.e-9;
        break;
      case eLifetime:
        // Convert from s to ns.
        m_eLifetime[ivertex] = val1 * 1.e9;
        break;
      case hLifetime:
        // Convert from s to ns.
        m_hLifetime[ivertex] = val1 * 1.e9;
        break;
      case DonorTrapOccupation:
        m_donorOcc[ivertex].push_back(val1);
        break;
      case AcceptorTrapOccupation:
        m_acceptorOcc[ivertex].push_back(val1);
        break;
      default:
        std::cerr << m_className << "::ReadDataset:\n"
                  << "    Unexpected dataset (" << ds << "). Program bug!\n";
        datafile.close();
        return false;
    }
    ++ivertex;
  }
  return true;
}

bool ComponentTcad2d::LoadGrid(const std::string& filename) {
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
  unsigned int iLine = 0;
  // Get the number of regions.
  size_t nRegions = 0;
  // Read the file line by line.
  std::string line;
  while (std::getline(gridfile, line)) {
    ++iLine;
    ltrim(line);
    // Find entry 'nb_regions'.
    if (line.substr(0, 10) != "nb_regions") continue;
    const std::string::size_type pEq = line.find('=');
    if (pEq == std::string::npos) {
      // No "=" sign found.
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of regions.\n";
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
              << "    " << filename << ".\n";
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    // Error reading from the file.
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
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
    if (!ExtractFromSquareBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read region names.\n";
      gridfile.close();
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
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    // Error reading from the file.
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
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
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of vertices.\n";
      gridfile.close();
      return false;
    }
    std::istringstream data;
    data.str(line);
    data >> nVertices;
    m_vertices.resize(nVertices);
    // Get the coordinates of this vertex.
    for (unsigned int j = 0; j < nVertices; ++j) {
      gridfile >> m_vertices[j][0] >> m_vertices[j][1];
      // Change units from micron to cm.
      m_vertices[j][0] *= 1.e-4;
      m_vertices[j][1] *= 1.e-4;
      ++iLine;
    }
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Vertices' in file\n"
              << "    " << filename << ".\n";
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
    gridfile.close();
    return false;
  }

  // Get the "edges" (lines connecting two vertices).
  size_t nEdges = 0;
  // Temporary arrays for storing edge points.
  std::vector<int> edgeP1;
  std::vector<int> edgeP2;
  while (std::getline(gridfile, line)) {
    ++iLine;
    ltrim(line);
    // Find section 'Edges'.
    if (line.substr(0, 5) != "Edges") continue;
    // Get the number of edges (given in brackets).
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Could not read number of edges.\n";
      gridfile.close();
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
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
    gridfile.close();
    return false;
  }

  for (size_t i = 0; i < nEdges; ++i) {
    // Make sure the indices of the edge endpoints are not out of range.
    if (edgeP1[i] < 0 || edgeP1[i] >= (int)nVertices || edgeP2[i] < 0 ||
        edgeP2[i] >= (int)nVertices) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Vertex index of edge " << i << " out of range.\n";
      gridfile.close();
      return false;
    }
    // Make sure the edge is non-degenerate.
    if (edgeP1[i] == edgeP2[i]) {
      std::cerr << m_className << "::LoadGrid:\n"
                << "    Edge " << i << " is degenerate.\n";
      gridfile.close();
      return false;
    }
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
      gridfile.close();
      return false;
    }
    std::istringstream data;
    data.str(line);
    data >> nElements;
    // Resize array of elements.
    m_elements.resize(nElements);
    // Get type and constituting edges of each element.
    for (size_t j = 0; j < nElements; ++j) {
      for (int k = nMaxVertices; k--;) m_elements[j].vertex[k] = -1;
      ++iLine;
      int type;
      gridfile >> type;
      int p0, p1, p2, p3;
      switch (type) {
        case 0:
          // Point
          gridfile >> p0 ;
          // Make sure the indices are not out of range.
          if (p0 >= (int)nVertices) {
            std::cerr << m_className << "::LoadGrid:\n"
                      << "    Error reading file " << filename << " (line "
                      << iLine << ").\n"
                      << "    Vertex index out of range.\n";
            gridfile.close();
            return false;
          }
          m_elements[j].vertex[0] = p0;
          break;
        case 1:
          // Line
          gridfile >> p0 >> p1;
          if (p0 < 0) p0 = -p0 - 1;
          if (p1 < 0) p1 = -p1 - 1;
          // Make sure the indices are not out of range.
          if (p0 >= (int)nVertices || p1 >= (int)nVertices) {
            std::cerr << m_className << "::LoadGrid:\n"
                      << "    Error reading file " << filename << " (line "
                      << iLine << ").\n"
                      << "    Vertex index out of range.\n";
            gridfile.close();
            return false;
          }
          m_elements[j].vertex[0] = p0;
          m_elements[j].vertex[1] = p1;
          break;
        case 2:
          // Triangle
          gridfile >> p0 >> p1 >> p2;
          // Negative edge index means that the sequence of the two points
          // is supposed to be inverted.
          // The actual index is then given by "-index - 1".
          if (p0 < 0) p0 = -p0 - 1;
          if (p1 < 0) p1 = -p1 - 1;
          if (p2 < 0) p2 = -p2 - 1;
          // Make sure the indices are not out of range.
          if (p0 >= (int)nEdges || p1 >= (int)nEdges || p2 >= (int)nEdges) {
            std::cerr << m_className << "::LoadGrid:\n"
                      << "    Error reading file " << filename << " (line "
                      << iLine << ").\n"
                      << "    Edge index out of range.\n";
            gridfile.close();
            return false;
          }
          m_elements[j].vertex[0] = edgeP1[p0];
          m_elements[j].vertex[1] = edgeP2[p0];
          if (edgeP1[p1] != (int)m_elements[j].vertex[0] &&
              edgeP1[p1] != (int)m_elements[j].vertex[1]) {
            m_elements[j].vertex[2] = edgeP1[p1];
          } else {
            m_elements[j].vertex[2] = edgeP2[p1];
          }
          // Rearrange vertices such that point 0 is on the left.
          while (m_vertices[m_elements[j].vertex[0]][0] >
                 m_vertices[m_elements[j].vertex[1]][0] ||
                 m_vertices[m_elements[j].vertex[0]][0] >
                 m_vertices[m_elements[j].vertex[2]][0]) {
            const int tmp = m_elements[j].vertex[0];
            m_elements[j].vertex[0] = m_elements[j].vertex[1];
            m_elements[j].vertex[1] = m_elements[j].vertex[2];
            m_elements[j].vertex[2] = tmp;
          }
          break;
        case 3:
          // Rectangle
          gridfile >> p0 >> p1 >> p2 >> p3;
          // Make sure the indices are not out of range.
          if (p0 >= (int)nEdges || -p0 - 1 >= (int)nEdges || 
              p1 >= (int)nEdges || -p1 - 1 >= (int)nEdges || 
              p2 >= (int)nEdges || -p2 - 1 >= (int)nEdges ||
              p3 >= (int)nEdges || -p3 - 1 >= (int)nEdges) {
            std::cerr << m_className << "::LoadGrid:\n"
                      << "    Error reading file " << filename << " (line "
                      << iLine << ").\n"
                      << "    Edge index out of range.\n";
            gridfile.close();
            return false;
          }
          if (p0 >= 0)
            m_elements[j].vertex[0] = edgeP1[p0];
          else
            m_elements[j].vertex[0] = edgeP2[-p0 - 1];
          if (p1 >= 0)
            m_elements[j].vertex[1] = edgeP1[p1];
          else
            m_elements[j].vertex[1] = edgeP2[-p1 - 1];
          if (p2 >= 0)
            m_elements[j].vertex[2] = edgeP1[p2];
          else
            m_elements[j].vertex[2] = edgeP2[-p2 - 1];
          if (p3 >= 0)
            m_elements[j].vertex[3] = edgeP1[p3];
          else
            m_elements[j].vertex[3] = edgeP2[-p3 - 1];

          // Rearrange vertices such that point 0 is on the left.
          while (m_vertices[m_elements[j].vertex[0]][0] >
                 m_vertices[m_elements[j].vertex[1]][0] ||
                 m_vertices[m_elements[j].vertex[0]][0] >
                 m_vertices[m_elements[j].vertex[2]][0] ||
                 m_vertices[m_elements[j].vertex[0]][0] >
                 m_vertices[m_elements[j].vertex[3]][0]) {
            const int tmp = m_elements[j].vertex[0];
            m_elements[j].vertex[0] = m_elements[j].vertex[1];
            m_elements[j].vertex[1] = m_elements[j].vertex[2];
            m_elements[j].vertex[2] = m_elements[j].vertex[3];
            m_elements[j].vertex[3] = tmp;
          }
          break;
        default:
          // Other element types are not permitted for 2d grids.
          std::cerr << m_className << "::LoadGrid:\n"
                    << "    Error reading file " << filename << " (line "
                    << iLine << ").\n";
          std::cerr << "    Invalid element type (" << type
                    << ") for 2d mesh.\n";
          gridfile.close();
          return false;
      }
      m_elements[j].type = type;
      m_elements[j].region = m_regions.size();
    }
    break;
  }
  if (gridfile.eof()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Could not find section 'Elements' in file "
              << filename << ".\n";
    gridfile.close();
    return false;
  } else if (gridfile.fail()) {
    std::cerr << m_className << "::LoadGrid:\n"
              << "    Error reading file " << filename << " (line " << iLine
              << ").\n";
    gridfile.close();
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
      Cleanup();
      gridfile.close();
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
      gridfile.close();
      return false;
    }
    int nElementsRegion;
    int iElement;
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
              << "    Error reading file " << filename << ".\n";
    return false;
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

  // The point is not in the previous element nor in the adjacent ones.
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

void ComponentTcad2d::Reset() {
  Cleanup();
  m_hasRangeZ = false;
  m_ready = false;
}

}
