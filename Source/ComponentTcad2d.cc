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

ComponentTcad2d::ComponentTcad2d() : Component("Tcad2d") {
  m_regions.reserve(10);
  m_vertices.reserve(10000);
  m_elements.reserve(10000);
}

bool ComponentTcad2d::SetDonor(const size_t donorNumber,
                               const double eXsec, const double hXsec,
                               const double conc) {
  if (donorNumber >= m_donors.size()) {
    std::cerr << m_className << "::SetDonor: Index out of range.\n";
    return false;
  }
  m_donors[donorNumber].xsece = eXsec;
  m_donors[donorNumber].xsech = hXsec;
  m_donors[donorNumber].conc = conc;

  UpdateAttachment();
  return true;
}

bool ComponentTcad2d::SetAcceptor(const size_t acceptorNumber,
                                  const double eXsec, const double hXsec,
                                  const double conc) {
  if (acceptorNumber >= m_acceptors.size()) {
    std::cerr << m_className << "::SetAcceptor: Index out of range.\n";
    return false;
  }
  m_acceptors[acceptorNumber].xsece = eXsec;
  m_acceptors[acceptorNumber].xsech = hXsec;
  m_acceptors[acceptorNumber].conc = conc;

  UpdateAttachment();
  return true;
}

bool ComponentTcad2d::HasAttachmentMap() const {
  return (m_useAttachmentMap && !(m_acceptors.empty() && m_donors.empty()));
}

bool ComponentTcad2d::ElectronAttachment(const double x, const double y,
                                         const double z, double& eta) {
  Interpolate(x, y, z, m_eAttachment, eta);
  return true;
}

bool ComponentTcad2d::HoleAttachment(const double x, const double y,
                                     const double z, double& eta) {
  Interpolate(x, y, z, m_hAttachment, eta);
  return true;
}

void ComponentTcad2d::WeightingField(const double x, const double y,
                                     const double z, double& wx, double& wy,
                                     double& wz, const std::string& /*label*/) {
  wx = wy = wz = 0.;
  if (m_wf.empty()) {
    std::cerr << m_className << "::WeightingField: Not available.\n";
    return;
  }
  Interpolate(x, y, z, m_wf, wx, wy);
}

double ComponentTcad2d::WeightingPotential(const double x, const double y,
                                           const double z, 
                                           const std::string& /*label*/) {

  if (m_wp.empty()) {
    std::cerr << m_className << "::WeightingPotential: Not available.\n";
    return 0.;
  }
  double v = 0.;
  Interpolate(x, y, z, m_wp, v);
  return v;
}

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
  double x = xin, y = yin;
  bool xmirr = false, ymirr = false;
  MapCoordinates(x, y, xmirr, ymirr);
  // Check if the point is inside the bounding box.
  if (!InBoundingBox(x, y)) {
    status = -6;
    return;
  }

  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x, y, w);
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
  if (xmirr) ex = -ex;
  if (ymirr) ey = -ey;
  m = m_regions[element.region].medium;
  if (!m_regions[element.region].drift || !m) status = -5;
  m_lastElement = i;
}

bool ComponentTcad2d::Interpolate(const double xin, const double yin,
    const double z,
    const std::vector<std::array<double, 2> >& field,
    double& fx, double& fy) {

  fx = fy = 0.;
  if (field.empty()) return false;
  if (m_hasRangeZ && (z < m_bbMin[2] || z > m_bbMax[2])) return false;
  double x = xin, y = yin;
  // In case of periodicity, reduce to the cell volume.
  bool xmirr = false, ymirr = false;
  MapCoordinates(x, y, xmirr, ymirr);
  // Make sure the point is inside the bounding box.
  if (!InBoundingBox(x, y)) return false;

  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x, y, w);
  // Stop if the point is outside the mesh.
  if (i >= m_elements.size()) return false;

  const Element& element = m_elements[i];
  const size_t nVertices = element.type + 1;
  for (size_t j = 0; j < nVertices; ++j) {
    const size_t index = element.vertex[j];
    fx += w[j] * field[index][0];
    fy += w[j] * field[index][1];
  }
  if (xmirr) fx = -fx;
  if (ymirr) fy = -fy;
  m_lastElement = i;
  return true;
}

bool ComponentTcad2d::Interpolate(const double xin, const double yin,
    const double z,
    const std::vector<double>& field, double& f) {

  f = 0.;
  if (field.empty()) return false;
  if (m_hasRangeZ && (z < m_bbMin[2] || z > m_bbMax[2])) return false;
  double x = xin, y = yin;
  // In case of periodicity, reduce to the cell volume.
  bool xmirr = false, ymirr = false;
  MapCoordinates(x, y, xmirr, ymirr);
  if (!InBoundingBox(x, y)) return false;

  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x, y, w);
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

bool ComponentTcad2d::ElectronVelocity(const double x, const double y,
                                       const double z, double& vx, double& vy,
                                       double& vz) {
  vz = 0.;
  return Interpolate(x, y, z, m_eVelocity, vx, vy);
}

bool ComponentTcad2d::HoleVelocity(const double x, const double y,
                                   const double z, double& vx, double& vy,
                                   double& vz) {
  vz = 0.;
  return Interpolate(x, y, z, m_hVelocity, vx, vy);
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
  double x = xin, y = yin;
  // In case of periodicity, reduce to the cell volume.
  bool xmirr = false, ymirr = false;
  MapCoordinates(x, y, xmirr, ymirr);
  // Check if the point is inside the bounding box.
  if (!InBoundingBox(x, y)) return nullptr;

  // Shape functions
  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x, y, w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    return nullptr;
  }
  m_lastElement = i;
  const Element& element = m_elements[i];
  return m_regions[element.region].medium;
}

bool ComponentTcad2d::GetElectronLifetime(const double x, const double y,
                                          const double z, double& tau) {
  return Interpolate(x, y, z, m_eLifetime, tau);
}
bool ComponentTcad2d::GetHoleLifetime(const double x, const double y,
                                      const double z, double& tau) {
  return Interpolate(x, y, z, m_hLifetime, tau);
}

bool ComponentTcad2d::GetElectronMobility(const double x, const double y,
                                          const double z, double& mob) {
  return Interpolate(x, y, z, m_eMobility, mob);
}

bool ComponentTcad2d::GetHoleMobility(const double x, const double y,
                                      const double z, double& mob) {
  return Interpolate(x, y, z, m_hMobility, mob);
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
  m_bbMax[0] = m_bbMin[0] = m_vertices[m_elements[0].vertex[0]][0];
  m_bbMax[1] = m_bbMin[1] = m_vertices[m_elements[0].vertex[0]][1];
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
    element.xmin = xmin - tol;
    element.xmax = xmax + tol;
    element.ymin = ymin - tol;
    element.ymax = ymax + tol;
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
    const double bb[4] = {e.xmin, e.ymin, e.xmax, e.ymax};
    m_tree->InsertMeshElement(bb, i);
  }

  m_ready = true;
  UpdatePeriodicity();
  std::cout << m_className << "::Initialise: Initialisation finished.\n";
  return true;
}

bool ComponentTcad2d::SetWeightingField(const std::string& datfile1,
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
  std::vector<std::array<double, 2> > wf1;
  std::vector<double> wp1;
  if (!LoadWeightingField(datfile1, wf1, wp1)) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Could not import data from " << datfile1 << ".\n";
    return false;
  }
  // Then load the field/potential for the configuration with the potential 
  // at the electrode to be read out increased by small voltage dv. 
  std::vector<std::array<double, 2> > wf2;
  std::vector<double> wy2;
  std::vector<double> wp2;
  if (!LoadWeightingField(datfile2, wf2, wp2)) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Could not import data from " << datfile2 << ".\n";
    return false;
  }
  const unsigned int nVertices = m_vertices.size();
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
    m_wf.assign(nVertices, {0., 0.});
    for (unsigned int i = 0; i < nVertices; ++i) {
      m_wf[i][0] = (wf2[i][0] - wf1[i][0]) * s; 
      m_wf[i][1] = (wf2[i][1] - wf1[i][1]) * s; 
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

bool ComponentTcad2d::GetVoltageRange(double& vmin, double& vmax) {
  if (!m_ready) return false;
  vmin = m_pMin;
  vmax = m_pMax;
  return true;
}

void ComponentTcad2d::PrintRegions() const {

  if (m_regions.empty()) {
    std::cerr << m_className << "::PrintRegions:\n"
              << "    No regions are currently defined.\n";
    return;
  }

  const size_t nRegions = m_regions.size();
  std::cout << m_className << "::PrintRegions:\n"
            << "    Currently " << nRegions << " regions are defined.\n"
            << "      Index  Name       Medium\n";
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

void ComponentTcad2d::GetRegion(const size_t i, std::string& name,
                                bool& active) const {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::GetRegion: Index out of range.\n";
    return;
  }
  name = m_regions[i].name;
  active = m_regions[i].drift;
}

void ComponentTcad2d::SetDriftRegion(const size_t i) {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::SetDriftRegion: Index out of range.\n";
    return;
  }
  m_regions[i].drift = true;
}

void ComponentTcad2d::UnsetDriftRegion(const size_t i) {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::UnsetDriftRegion: Index out of range.\n";
    return;
  }
  m_regions[i].drift = false;
}

void ComponentTcad2d::SetMedium(const size_t i, Medium* medium) {
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

Medium* ComponentTcad2d::GetMedium(const size_t i) const {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::GetMedium: Index out of range.\n";
    return nullptr;
  }

  return m_regions[i].medium;
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

bool ComponentTcad2d::LoadWeightingField(const std::string& datafilename,
    std::vector<std::array<double, 2> >& wf, std::vector<double>& wp) {
  std::ifstream datafile;
  datafile.open(datafilename.c_str(), std::ios::in);
  if (!datafile) {
    std::cerr << m_className << "::LoadWeightingField:\n"
              << "    Could not open file " << datafilename << ".\n";
    return false;
  }

  const unsigned int nVertices = m_vertices.size();
  bool ok = true;
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
    if (dataset != "ElectrostaticPotential" && dataset != "ElectricField") {
      continue;
    }
    bool field = false;
    if (dataset == "ElectricField") {
      wf.assign(nVertices, {0., 0.});
      field = true;
    } else {
      wp.assign(nVertices, 0.);
    }
    std::getline(datafile, line);
    std::getline(datafile, line);
    std::getline(datafile, line);
    std::getline(datafile, line);
    // Get the region name (given in brackets).
    if (!ExtractFromSquareBrackets(line)) {
      std::cerr << m_className << "::LoadWeightingField:\n"
                << "    Cannot extract region name.\n"
                << "    Line:\n    " << line << "\n";
      ok = false;
      break;
    }
    std::string name;
    data.str(line);
    data >> name;
    data.clear();
    // Check if the region name matches one from the mesh file.
    const auto index = FindRegion(name);
    if (index >= m_regions.size()) {
      std::cerr << m_className << "::LoadWeightingField:\n"
                << "    Unknown region " << name << ".\n";
      ok = false;
      break;
    }
    // Get the number of values.
    std::getline(datafile, line);
    if (!ExtractFromBrackets(line)) {
      std::cerr << m_className << "::LoadWeightingField:\n"
                << "    Cannot extract number of values to be read.\n"
                << "    Line:\n    " << line << "\n";
      ok = false;
      break;
    }
    int nValues;
    data.str(line);
    data >> nValues;
    data.clear();
    if (field) nValues /= 2;
    // Mark the vertices belonging to this region.
    std::vector<bool> isInRegion(nVertices, false);
    const unsigned int nElements = m_elements.size();
    for (unsigned int j = 0; j < nElements; ++j) {
      if (m_elements[j].region != index) continue;
      for (int k = 0; k <= m_elements[j].type; ++k) {
        isInRegion[m_elements[j].vertex[k]] = true;
      }
    }
    unsigned int ivertex = 0;
    for (int j = 0; j < nValues; ++j) {
      // Read the next value.
      double val1 = 0., val2 = 0.;
      if (field) {
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
        std::cerr << m_className << "::LoadWeightingField:\n"
                  << "    Dataset " << dataset
                  << " has more values than vertices in region " << name << "\n";
        ok = false;
        break;
      }
      if (field) {
        wf[ivertex] = {val1, val2};
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
          if (edgeP1[p1] != m_elements[j].vertex[0] &&
              edgeP1[p1] != m_elements[j].vertex[1]) {
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

void ComponentTcad2d::Cleanup() {
  // Vertices
  m_vertices.clear();
  // Elements
  m_elements.clear();
  // Regions
  m_regions.clear();
  // Potential and electric field.
  m_potential.clear();
  m_efield.clear();
  // Weighting fields and potentials
  m_wf.clear();
  m_wp.clear();
  // Other data.
  m_eVelocity.clear();
  m_hVelocity.clear();
  m_eMobility.clear();
  m_hMobility.clear();
  m_eLifetime.clear();
  m_hLifetime.clear();
  m_donors.clear();
  m_acceptors.clear();
  m_donorOcc.clear();
  m_acceptorOcc.clear();
  m_eAttachment.clear();
  m_hAttachment.clear();
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
  if (x < element.xmin || x > element.xmax || y < element.ymin ||
      y > element.ymax) {
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

void ComponentTcad2d::UpdatePeriodicity() {
  if (!m_ready) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Field map not available.\n";
    return;
  }

  // Check for conflicts.
  for (size_t i = 0; i < 3; ++i) {
    if (m_periodic[i] && m_mirrorPeriodic[i]) {
      std::cerr << m_className << "::UpdatePeriodicity:\n"
                << "    Both simple and mirror periodicity requested. Reset.\n";
      m_periodic[i] = m_mirrorPeriodic[i] = false;
    }
    if (m_axiallyPeriodic[i]) {
      std::cerr << m_className << "::UpdatePeriodicity:\n"
                 << "    Axial symmetry is not supported. Reset.\n";
      m_axiallyPeriodic.fill(false);
    }
    if (m_rotationSymmetric[i]) {
      std::cerr << m_className << "::UpdatePeriodicity:\n"
                << "    Rotation symmetry is not supported. Reset.\n";
      m_rotationSymmetric.fill(false);
    }
  }
}

size_t ComponentTcad2d::FindRegion(const std::string& name) const {
  const auto nRegions = m_regions.size();
  for (size_t j = 0; j < nRegions; ++j) {
    if (name == m_regions[j].name) return j;
  }
  return m_regions.size();
}

void ComponentTcad2d::MapCoordinates(double& x, double& y, bool& xmirr,
                                     bool& ymirr) const {
  // In case of periodicity, reduce to the cell volume.
  xmirr = false;
  const double cellsx = m_bbMax[0] - m_bbMin[0];
  if (m_periodic[0]) {
    x = m_bbMin[0] + fmod(x - m_bbMin[0], cellsx);
    if (x < m_bbMin[0]) x += cellsx;
  } else if (m_mirrorPeriodic[0]) {
    double xNew = m_bbMin[0] + fmod(x - m_bbMin[0], cellsx);
    if (xNew < m_bbMin[0]) xNew += cellsx;
    const int nx = int(floor(0.5 + (xNew - x) / cellsx));
    if (nx != 2 * (nx / 2)) {
      xNew = m_bbMin[0] + m_bbMax[0] - xNew;
      xmirr = true;
    }
    x = xNew;
  }
  ymirr = false;
  const double cellsy = m_bbMax[1] - m_bbMin[1];
  if (m_periodic[1]) {
    y = m_bbMin[1] + fmod(y - m_bbMin[1], cellsy);
    if (y < m_bbMin[1]) y += cellsy;
  } else if (m_mirrorPeriodic[1]) {
    double yNew = m_bbMin[1] + fmod(y - m_bbMin[1], cellsy);
    if (yNew < m_bbMin[1]) yNew += cellsy;
    const int ny = int(floor(0.5 + (yNew - y) / cellsy));
    if (ny != 2 * (ny / 2)) {
      yNew = m_bbMin[1] + m_bbMax[1] - yNew;
      ymirr = true;
    }
    y = yNew;
  }
}

void ComponentTcad2d::UpdateAttachment() {

  if (m_vertices.empty()) return;
  const size_t nVertices = m_vertices.size();
  m_eAttachment.assign(nVertices, 0.);
  m_hAttachment.assign(nVertices, 0.);

  const size_t nAcceptors = m_acceptors.size();
  for (size_t i = 0; i < nAcceptors; ++i) { 
    const auto& defect = m_acceptors[i];
    if (defect.conc < 0.) continue;
    for (size_t j = 0; j < nVertices; ++j) {
      // Get the occupation probability.
      const double f = m_acceptorOcc[j][i];
      if (defect.xsece > 0.) {
        m_eAttachment[j] += defect.conc * defect.xsece * (1. - f);
      }
      if (defect.xsech > 0.) {
        m_hAttachment[j] += defect.conc * defect.xsech * f;
      }
    }
  }
  const size_t nDonors = m_donors.size();
  for (size_t i = 0; i < nDonors; ++i) {
    const auto& defect = m_donors[i];
    if (defect.conc < 0.) continue;
    for (size_t j = 0; j < nVertices; ++j) { 
      const double f = m_donorOcc[j][i];
      if (defect.xsece > 0.) {
        m_eAttachment[j] += defect.conc * defect.xsece * f;
      }
      if (defect.xsech > 0.) {
        m_hAttachment[j] += defect.conc * defect.xsech * (1. - f);
      }
    }
  }
}

}
