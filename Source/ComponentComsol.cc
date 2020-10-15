#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "Garfield/ComponentComsol.hh"
#include "Garfield/KDTree.hh"

namespace {

bool ends_with(std::string s, std::string t) {
  if (!s.empty() && s.back() == '\r') s.pop_back();
  return s.size() >= t.size() && s.substr(s.size() - t.size(), t.size()) == t;
}

int readInt(std::string s) {
  std::istringstream iss(s);
  int ret;
  iss >> ret;
  return ret;
}

void PrintProgress(const double f) {
  if (f < 0.) return;
  constexpr unsigned int width = 70;
  const unsigned int n = static_cast<unsigned int>(std::floor(width * f));
  std::string bar = "[";
  if (n < 1) {
    bar += std::string(width, ' ');
  } else if (n >= width) {
    bar += std::string(width, '=');
  } else {
    bar += std::string(n, '=') + ">" + std::string(width - n - 1, ' ');
  }
  bar += "]";
  std::cout << bar << "\r" << std::flush;
}

}

namespace Garfield {

ComponentComsol::ComponentComsol() : ComponentFieldMap("Comsol") {}

ComponentComsol::ComponentComsol(const std::string& mesh, 
                                 const std::string& mplist,
                                 const std::string& field, 
                                 const std::string& unit)
    : ComponentComsol() {
  Initialise(mesh, mplist, field, unit);
}

bool ComponentComsol::Initialise(const std::string& mesh, 
                                 const std::string& mplist,
                                 const std::string& field, 
                                 const std::string& unit) {
  m_ready = false;
  m_warning = false;
  m_nWarnings = 0;

  // Get the conversion factor to be applied to the coordinates.
  m_unit = ScalingFactor(unit);
  if (m_unit <= 0.) {
    std::cerr << m_className << "::Initialise:\n    Unknown length unit " 
              << unit << ". Will use default (m).\n";
    m_unit = 100.;
  }
  // Open the materials file.
  m_materials.clear();
  std::ifstream fmplist;
  fmplist.open(mplist.c_str(), std::ios::in);
  if (fmplist.fail()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Could not open materials file " << mplist << ".\n";
    return false;
  }
  unsigned int nMaterials;
  fmplist >> nMaterials;
  for (unsigned int i = 0; i < nMaterials; ++i) {
    Material material;
    material.driftmedium = false;
    material.medium = nullptr;
    material.ohm = -1;
    fmplist >> material.eps;
    m_materials.push_back(std::move(material));
  }
  if (m_materials.empty()) {
    // Add default material
    Material material;
    material.driftmedium = false;
    material.medium = nullptr;
    material.eps = material.ohm = -1;
    m_materials.push_back(std::move(material));
    nMaterials = 1;
  }
  std::map<int, int> domain2material;
  int d2msize;
  fmplist >> d2msize;
  for (int i = 0; i < d2msize; ++i) {
    int domain;
    fmplist >> domain;
    fmplist >> domain2material[domain];
  }
  fmplist.close();

  // Find the lowest epsilon, check for eps = 0, set default drift media.
  double epsmin = -1.;
  unsigned int iepsmin = 0;
  for (unsigned int i = 0; i < nMaterials; ++i) {
    if (m_materials[i].eps < 0) continue;
    if (m_materials[i].eps == 0) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Material " << i 
                << " has been assigned a permittivity equal to zero in\n    "
                << mplist << ".\n";
      m_materials[i].eps = -1.;
    } else if (epsmin < 0. || epsmin > m_materials[i].eps) {
      epsmin = m_materials[i].eps;
      iepsmin = i;
    }
  }

  if (epsmin < 0.) {
    std::cerr << m_className << "::Initialise:\n"
              << "    No material with positive permittivity found \n"
              << "    in material list " << mplist << ".\n";
    return false;
  }
  m_materials[iepsmin].driftmedium = true;

  m_nodes.clear();
  std::ifstream fmesh;
  fmesh.open(mesh.c_str(), std::ios::in);
  if (fmesh.fail()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Could not open nodes file " << mesh << ".\n";
    return false;
  }

  std::string line;
  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Could not read number of nodes from " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (!ends_with(line, "# number of mesh points") && 
           !ends_with(line, "# number of mesh vertices"));
  const int nNodes = readInt(line);

  std::cout << m_className << "::Initialise: " << nNodes << " nodes.\n";
  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("# Mesh point coordinates") == std::string::npos &&
           line.find("# Mesh vertex coordinates") == std::string::npos);
  for (int i = 0; i < nNodes; ++i) {
    Node newNode;
    fmesh >> newNode.x >> newNode.y >> newNode.z;
    newNode.x *= m_unit;
    newNode.y *= m_unit;
    newNode.z *= m_unit;
    m_nodes.push_back(std::move(newNode));
  }

  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("4 tet2 # type name") == std::string::npos);
  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (!ends_with(line, "# number of elements"));
  const int nElements = readInt(line);
  m_elements.clear();
  std::cout << m_className << "::Initialise: " << nElements << " elements.\n";
  std::getline(fmesh, line);
  // Elements 6 & 7 are swapped due to differences in COMSOL and ANSYS
  // representation
  int perm[10] = {0, 1, 2, 3, 4, 5, 7, 6, 8, 9};
  for (int i = 0; i < nElements; ++i) {
    Element newElement;
    newElement.degenerate = false;
    for (int j = 0; j < 10; ++j) {
      fmesh >> newElement.emap[perm[j]];
    }
    m_elements.push_back(std::move(newElement));
  }

  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("# Geometric entity indices") == std::string::npos);
  for (auto& element : m_elements) {
    int domain;
    fmesh >> domain;
    element.matmap = domain2material.count(domain) ? domain2material[domain]
                                                   : nMaterials - 1;
  }
  fmesh.close();

  std::ifstream ffield;
  ffield.open(field.c_str(), std::ios::in);
  if (ffield.fail()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Could not open potentials file " << field << ".\n";
    return false;
  }

  const std::string hdr1 = "% x                       y                        z                        V (V)";

  const std::string hdr2 = "% x             y              z              V (V)";
  do {
    if (!std::getline(ffield, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << field << ".\n";
      ffield.close();
      return false;
    }
  } while (line.find(hdr1) == std::string::npos && 
           line.find(hdr2) == std::string::npos);
  std::istringstream sline(line);
  std::string token;
  sline >> token;  // %
  sline >> token;  // x
  sline >> token;  // y
  sline >> token;  // z
  sline >> token;  // V
  sline >> token;  // (V)
  m_wfields.clear();
  m_wfieldsOk.clear();
  while (sline >> token) {
    std::cout << m_className << "::Initialise:\n"
              << "    Reading data for weighting field " << token << ".\n";
    m_wfields.push_back(token);
    m_wfieldsOk.push_back(true);
    sline >> token;  // (V)
  }
  const size_t nWeightingFields = m_wfields.size();

  const unsigned int nPrint =
      std::pow(10, static_cast<unsigned int>(
                       std::max(std::floor(std::log10(nNodes)) - 1, 1.)));
  std::cout << m_className << "::Initialise: Reading potentials.\n";
  PrintProgress(0.);
  // Build a k-d tree from the node coordinates.
  std::vector<std::vector<double> > points;
  for (const auto& node : m_nodes) {
    std::vector<double> point = {node.x, node.y, node.z};
    points.push_back(std::move(point));
  }
  KDTree kdtree(points);
  std::vector<bool> used(nNodes, false);
  for (int i = 0; i < nNodes; ++i) {
    double x, y, z, v;
    ffield >> x >> y >> z >> v;
    x *= m_unit;
    y *= m_unit;
    z *= m_unit;
    std::vector<double> w;
    for (size_t j = 0; j < nWeightingFields; ++j) {
      double p;
      ffield >> p;
      w.push_back(p);
    }
    const std::vector<double> pt = {x, y, z};
    std::vector<KDTreeResult> res;
    kdtree.n_nearest(pt, 1, res);
    if (res.empty()) {
      std::cerr << std::endl << m_className << "::Initialise:\n"
                << "    Could not find a matching mesh node for point ("
                << x << ", " << y << ", " << z << ")\n.";
      ffield.close();
      return false;
    }
    const size_t k = res[0].idx;
    used[k] = true;
    m_nodes[k].v = v;
    m_nodes[k].w = w;
    if ((i + 1) % nPrint == 0) PrintProgress(double(i + 1) / nNodes);
  }
  PrintProgress(1.);
  ffield.close();
  auto nMissing = std::count(used.begin(), used.end(), false);
  if (nMissing > 0) {
    std::cerr << std::endl << m_className << "::Initialise:\n"
              << "    Missing potentials for " << nMissing << " nodes.\n";
    return false;
  }
  std::cout << std::endl << m_className << "::Initialise: Done.\n";

  m_ready = true;
  // Establish the ranges.
  SetRange();
  UpdatePeriodicity();
  return true;
}

bool ComponentComsol::SetWeightingField(const std::string& field, 
                                        const std::string& label) {

  if (!m_ready) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    No valid field map is present.\n"
              << "    Weighting fields cannot be added.\n";
    return false;
  }

  // Open the voltage list.
  std::ifstream ffield;
  ffield.open(field.c_str(), std::ios::in);
  if (ffield.fail()) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Could not open potentials file " << field << ".\n";
    return false;
  }

  // Check if a weighting field with the same label already exists.
  const size_t iw = GetOrCreateWeightingFieldIndex(label);
  if (iw + 1 != m_wfields.size()) {
    std::cout << m_className << "::SetWeightingField:\n"
              << "    Replacing existing weighting field " << label << ".\n";
  }
  m_wfieldsOk[iw] = false;

  // Build a k-d tree from the node coordinates.
  std::vector<std::vector<double> > points;
  for (const auto& node : m_nodes) {
    std::vector<double> point = {node.x, node.y, node.z};
    points.push_back(std::move(point));
  }
  KDTree kdtree(points);

  const std::string hdr = "% x                       y                        z                        V (V)";
  std::string line;
  do {
    if (!std::getline(ffield, line)) {
      std::cerr << m_className << "::SetWeightingField:\n"
                << "    Error parsing " << field << ".\n";
      ffield.close();
      return false;
    }
  } while (line.find(hdr) == std::string::npos);
  const int nNodes = m_nodes.size();
  for (int i = 0; i < nNodes; ++i) {
    double x, y, z, v;
    ffield >> x >> y >> z >> v;
    x *= m_unit;
    y *= m_unit;
    z *= m_unit;
    // Find the closest mesh node.
    const std::vector<double> pt = {x, y, z};
    std::vector<KDTreeResult> res;
    kdtree.n_nearest(pt, 1, res);
    if (res.empty()) {
      std::cerr << m_className << "::SetWeightingField:\n"
                << "    Could not find a matching mesh node for point ("
                << x << ", " << y << ", " << z << ")\n.";
      ffield.close();
      return false;
    }
    const size_t k = res[0].idx;
    m_nodes[k].w[iw] = v;
  }
  ffield.close();
  return true;
}

void ComponentComsol::ElectricField(const double x, const double y,
                                    const double z, double& ex, double& ey,
                                    double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

void ComponentComsol::ElectricField(const double xin, const double yin,
                                    const double zin, double& ex, double& ey,
                                    double& ez, double& volt, Medium*& m,
                                    int& status) {
  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  // Initial values
  ex = ey = ez = volt = 0.;
  status = 0;
  m = nullptr;

  // Do not proceed if not properly initialised.
  if (!m_ready) {
    status = -10;
    PrintNotReady("ElectricField");
    return;
  }

  if (m_warning) PrintWarning("ElectricField");

  // Find the element that contains this point
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField:\n"
                << "    Point (" << x << ", " << y << ", " << z
                << " not in the mesh.\n";
    }
    status = -6;
    return;
  }

  const Element& element = m_elements[imap];
  if (m_debug) {
    PrintElement("ElectricField", x, y, z, t1, t2, t3, t4, element, 10);
  }
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];
  const Node& n6 = m_nodes[element.emap[6]];
  const Node& n7 = m_nodes[element.emap[7]];
  const Node& n8 = m_nodes[element.emap[8]];
  const Node& n9 = m_nodes[element.emap[9]];
  // Tetrahedral field
  volt = n0.v * t1 * (2 * t1 - 1) + n1.v * t2 * (2 * t2 - 1) +
         n2.v * t3 * (2 * t3 - 1) + n3.v * t4 * (2 * t4 - 1) +
         4 * n4.v * t1 * t2 + 4 * n5.v * t1 * t3 + 4 * n6.v * t1 * t4 +
         4 * n7.v * t2 * t3 + 4 * n8.v * t2 * t4 + 4 * n9.v * t3 * t4;
  ex = -(n0.v * (4 * t1 - 1) * jac[0][1] + n1.v * (4 * t2 - 1) * jac[1][1] +
         n2.v * (4 * t3 - 1) * jac[2][1] + n3.v * (4 * t4 - 1) * jac[3][1] +
         n4.v * (4 * t2 * jac[0][1] + 4 * t1 * jac[1][1]) +
         n5.v * (4 * t3 * jac[0][1] + 4 * t1 * jac[2][1]) +
         n6.v * (4 * t4 * jac[0][1] + 4 * t1 * jac[3][1]) +
         n7.v * (4 * t3 * jac[1][1] + 4 * t2 * jac[2][1]) +
         n8.v * (4 * t4 * jac[1][1] + 4 * t2 * jac[3][1]) +
         n9.v * (4 * t4 * jac[2][1] + 4 * t3 * jac[3][1])) /
       det;
  ey = -(n0.v * (4 * t1 - 1) * jac[0][2] + n1.v * (4 * t2 - 1) * jac[1][2] +
         n2.v * (4 * t3 - 1) * jac[2][2] + n3.v * (4 * t4 - 1) * jac[3][2] +
         n4.v * (4 * t2 * jac[0][2] + 4 * t1 * jac[1][2]) +
         n5.v * (4 * t3 * jac[0][2] + 4 * t1 * jac[2][2]) +
         n6.v * (4 * t4 * jac[0][2] + 4 * t1 * jac[3][2]) +
         n7.v * (4 * t3 * jac[1][2] + 4 * t2 * jac[2][2]) +
         n8.v * (4 * t4 * jac[1][2] + 4 * t2 * jac[3][2]) +
         n9.v * (4 * t4 * jac[2][2] + 4 * t3 * jac[3][2])) /
       det;
  ez = -(n0.v * (4 * t1 - 1) * jac[0][3] + n1.v * (4 * t2 - 1) * jac[1][3] +
         n2.v * (4 * t3 - 1) * jac[2][3] + n3.v * (4 * t4 - 1) * jac[3][3] +
         n4.v * (4 * t2 * jac[0][3] + 4 * t1 * jac[1][3]) +
         n5.v * (4 * t3 * jac[0][3] + 4 * t1 * jac[2][3]) +
         n6.v * (4 * t4 * jac[0][3] + 4 * t1 * jac[3][3]) +
         n7.v * (4 * t3 * jac[1][3] + 4 * t2 * jac[2][3]) +
         n8.v * (4 * t4 * jac[1][3] + 4 * t2 * jac[3][3]) +
         n9.v * (4 * t4 * jac[2][3] + 4 * t3 * jac[3][3])) /
       det;

  // Transform field to global coordinates
  UnmapFields(ex, ey, ez, x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  // Drift medium?
  if (m_debug) {
    std::cout << m_className << "::ElectricField:\n"
              << "    Material " << element.matmap << ", drift flag "
              << m_materials[element.matmap].driftmedium << "\n";
  }
  m = m_materials[element.matmap].medium;
  status = -5;
  if (m_materials[element.matmap].driftmedium) {
    if (m && m->IsDriftable()) status = 0;
  }
}

void ComponentComsol::WeightingField(const double xin, const double yin,
                                     const double zin, double& wx, double& wy,
                                     double& wz, const std::string& label) {
  // Initial values
  wx = wy = wz = 0;

  // Do not proceed if not properly initialised.
  if (!m_ready) return;

  // Look for the label.
  const size_t iw = GetWeightingFieldIndex(label);
  // Do not proceed if the requested weighting field does not exist.
  if (iw == m_wfields.size()) return;
  // Check if the weighting field is properly initialised.
  if (!m_wfieldsOk[iw]) return;

  // Copy the coordinates.
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  if (m_warning) PrintWarning("WeightingField");

  // Find the element that contains this point.
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  // Check if the point is in the mesh.
  if (imap < 0) return;

  const Element& element = m_elements[imap];
  if (m_debug) {
    PrintElement("WeightingField", x, y, z, t1, t2, t3, t4, element, 10, iw);
  }
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];
  const Node& n6 = m_nodes[element.emap[6]];
  const Node& n7 = m_nodes[element.emap[7]];
  const Node& n8 = m_nodes[element.emap[8]];
  const Node& n9 = m_nodes[element.emap[9]];
  // Tetrahedral field
  wx = -(n0.w[iw] * (4 * t1 - 1) * jac[0][1] +
         n1.w[iw] * (4 * t2 - 1) * jac[1][1] +
         n2.w[iw] * (4 * t3 - 1) * jac[2][1] +
         n3.w[iw] * (4 * t4 - 1) * jac[3][1] +
         n4.w[iw] * (4 * t2 * jac[0][1] + 4 * t1 * jac[1][1]) +
         n5.w[iw] * (4 * t3 * jac[0][1] + 4 * t1 * jac[2][1]) +
         n6.w[iw] * (4 * t4 * jac[0][1] + 4 * t1 * jac[3][1]) +
         n7.w[iw] * (4 * t3 * jac[1][1] + 4 * t2 * jac[2][1]) +
         n8.w[iw] * (4 * t4 * jac[1][1] + 4 * t2 * jac[3][1]) +
         n9.w[iw] * (4 * t4 * jac[2][1] + 4 * t3 * jac[3][1])) /
       det;

  wy = -(n0.w[iw] * (4 * t1 - 1) * jac[0][2] +
         n1.w[iw] * (4 * t2 - 1) * jac[1][2] +
         n2.w[iw] * (4 * t3 - 1) * jac[2][2] +
         n3.w[iw] * (4 * t4 - 1) * jac[3][2] +
         n4.w[iw] * (4 * t2 * jac[0][2] + 4 * t1 * jac[1][2]) +
         n5.w[iw] * (4 * t3 * jac[0][2] + 4 * t1 * jac[2][2]) +
         n6.w[iw] * (4 * t4 * jac[0][2] + 4 * t1 * jac[3][2]) +
         n7.w[iw] * (4 * t3 * jac[1][2] + 4 * t2 * jac[2][2]) +
         n8.w[iw] * (4 * t4 * jac[1][2] + 4 * t2 * jac[3][2]) +
         n9.w[iw] * (4 * t4 * jac[2][2] + 4 * t3 * jac[3][2])) /
       det;

  wz = -(n0.w[iw] * (4 * t1 - 1) * jac[0][3] +
         n1.w[iw] * (4 * t2 - 1) * jac[1][3] +
         n2.w[iw] * (4 * t3 - 1) * jac[2][3] +
         n3.w[iw] * (4 * t4 - 1) * jac[3][3] +
         n4.w[iw] * (4 * t2 * jac[0][3] + 4 * t1 * jac[1][3]) +
         n5.w[iw] * (4 * t3 * jac[0][3] + 4 * t1 * jac[2][3]) +
         n6.w[iw] * (4 * t4 * jac[0][3] + 4 * t1 * jac[3][3]) +
         n7.w[iw] * (4 * t3 * jac[1][3] + 4 * t2 * jac[2][3]) +
         n8.w[iw] * (4 * t4 * jac[1][3] + 4 * t2 * jac[3][3]) +
         n9.w[iw] * (4 * t4 * jac[2][3] + 4 * t3 * jac[3][3])) /
       det;

  // Transform field to global coordinates
  UnmapFields(wx, wy, wz, x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);
}

double ComponentComsol::WeightingPotential(const double xin, const double yin,
                                           const double zin,
                                           const std::string& label) {
  // Do not proceed if not properly initialised.
  if (!m_ready) return 0.;

  // Look for the label.
  const size_t iw = GetWeightingFieldIndex(label);
  // Do not proceed if the requested weighting field does not exist.
  if (iw == m_wfields.size()) return 0.;
  // Check if the weighting field is properly initialised.
  if (!m_wfieldsOk[iw]) return 0.;

  // Copy the coordinates.
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates.
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  if (m_warning) PrintWarning("WeightingPotential");

  // Find the element that contains this point.
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) return 0.;

  const Element& element = m_elements[imap];
  if (m_debug) {
    PrintElement("WeightingPotential", x, y, z, t1, t2, t3, t4, element, 10,
                 iw);
  }
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];
  const Node& n6 = m_nodes[element.emap[6]];
  const Node& n7 = m_nodes[element.emap[7]];
  const Node& n8 = m_nodes[element.emap[8]];
  const Node& n9 = m_nodes[element.emap[9]];
  // Tetrahedral field
  return n0.w[iw] * t1 * (2 * t1 - 1) + n1.w[iw] * t2 * (2 * t2 - 1) +
         n2.w[iw] * t3 * (2 * t3 - 1) + n3.w[iw] * t4 * (2 * t4 - 1) +
         4 * n4.w[iw] * t1 * t2 + 4 * n5.w[iw] * t1 * t3 +
         4 * n6.w[iw] * t1 * t4 + 4 * n7.w[iw] * t2 * t3 +
         4 * n8.w[iw] * t2 * t4 + 4 * n9.w[iw] * t3 * t4;
}

Medium* ComponentComsol::GetMedium(const double xin, const double yin,
                                   const double zin) {
  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  // Do not proceed if not properly initialised.
  if (!m_ready) {
    PrintNotReady("GetMedium");
    return nullptr;
  }
  if (m_warning) PrintWarning("GetMedium");

  // Find the element that contains this point
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) {
    if (m_debug) {
      std::cout << m_className << "::GetMedium:\n"
                << "    Point (" << x << ", " << y << ", " << z
                << ") not in the mesh.\n";
    }
    return nullptr;
  }
  const Element& element = m_elements[imap];
  if (element.matmap >= m_materials.size()) {
    if (m_debug) {
      std::cerr << m_className << "::GetMedium:\n"
                << "    Point (" << x << ", " << y << ", " << z
                << ") has out of range material number " << imap << ".\n";
    }
    return nullptr;
  }

  if (m_debug) {
    PrintElement("GetMedium", x, y, z, t1, t2, t3, t4, element, 10);
  }

  return m_materials[element.matmap].medium;
}

double ComponentComsol::GetElementVolume(const unsigned int i) {
  if (i >= m_elements.size()) return 0.;
  const Element& element = m_elements[i];
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];

  // Uses formula V = |a (dot) b x c|/6
  // with a => "3", b => "1", c => "2" and origin "0"
  const double vol =
      fabs((n3.x - n0.x) *
               ((n1.y - n0.y) * (n2.z - n0.z) - (n2.y - n0.y) * (n1.z - n0.z)) +
           (n3.y - n0.y) *
               ((n1.z - n0.z) * (n2.x - n0.x) - (n2.z - n0.z) * (n1.x - n0.x)) +
           (n3.z - n0.z) * ((n1.x - n0.x) * (n2.y - n0.y) -
                            (n3.x - n0.x) * (n1.y - n0.y))) /
      6.;
  return vol;
}

void ComponentComsol::GetAspectRatio(const unsigned int i, double& dmin,
                                     double& dmax) {
  if (i >= m_elements.size()) {
    dmin = dmax = 0.;
    return;
  }

  const Element& element = m_elements[i];
  const int np = 4;
  // Loop over all pairs of vertices.
  for (int j = 0; j < np - 1; ++j) {
    const Node& nj = m_nodes[element.emap[j]];
    for (int k = j + 1; k < np; ++k) {
      const Node& nk = m_nodes[element.emap[k]];
      // Compute distance.
      const double dx = nj.x - nk.x;
      const double dy = nj.y - nk.y;
      const double dz = nj.z - nk.z;
      const double dist = sqrt(dx * dx + dy * dy + dz * dz);
      if (k == 1) {
        dmin = dmax = dist;
      } else {
        if (dist < dmin) dmin = dist;
        if (dist > dmax) dmax = dist;
      }
    }
  }
}

}  // namespace Garfield
