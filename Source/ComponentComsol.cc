// Copied and modified ComponentAnsys123.cc

#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "Garfield/ComponentComsol.hh"

namespace Garfield {

ComponentComsol::ComponentComsol() : ComponentFieldMap() {
  m_className = "ComponentComsol";
}

ComponentComsol::ComponentComsol(std::string mesh, std::string mplist,
                                 std::string field)
    : ComponentFieldMap() {
  m_className = "ComponentComsol";
  Initialise(mesh, mplist, field);
}

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

bool ComponentComsol::Initialise(std::string mesh, std::string mplist,
                                 std::string field) {
  m_ready = false;
  m_warning = false;
  m_nWarnings = 0;

  double unit = 100.0;  // m

  std::string line;

  // Open the materials file.
  materials.clear();
  std::ifstream fmplist;
  fmplist.open(mplist.c_str(), std::ios::in);
  if (fmplist.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open materials file " << mplist
              << " for reading.\n";
    return false;
  }
  fmplist >> m_nMaterials;
  for (unsigned int i = 0; i < m_nMaterials; ++i) {
    Material newMaterial;
    newMaterial.driftmedium = true;
    newMaterial.medium = nullptr;
    newMaterial.ohm = -1;
    fmplist >> newMaterial.eps;
    materials.push_back(newMaterial);
  }
  {
    // add default material
    Material newMaterial;
    newMaterial.driftmedium = false;
    newMaterial.medium = nullptr;
    newMaterial.eps = newMaterial.ohm = -1;
    materials.push_back(newMaterial);
    m_nMaterials++;
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

  nodes.clear();
  std::ifstream fmesh;
  fmesh.open(mesh.c_str(), std::ios::in);
  if (fmesh.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open nodes file " << mesh << " for reading.\n";
    return false;
  }

  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Could not read number of nodes from " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (!ends_with(line, "# number of mesh points"));
  nNodes = readInt(line);

  std::cout << m_className << "::Initialise:\n";
  std::cout << "    Reading " << nNodes << " nodes from " << mesh << ".\n";
  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error reading " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("# Mesh point coordinates") == std::string::npos);
  double minx = +std::numeric_limits<double>::max();
  double maxx = -std::numeric_limits<double>::max();
  double miny = minx, minz = minx;
  double maxy = maxx, maxz = maxx;
  for (int i = 0; i < nNodes; ++i) {
    Node newNode;
    fmesh >> newNode.x >> newNode.y >> newNode.z;
    newNode.x *= unit;
    newNode.y *= unit;
    newNode.z *= unit;
    nodes.push_back(std::move(newNode));
    minx = std::min(minx, newNode.x);
    maxx = std::max(maxx, newNode.x);
    miny = std::min(miny, newNode.y);
    maxy = std::max(maxy, newNode.y);
    minz = std::min(minz, newNode.z);
    maxz = std::max(maxz, newNode.z);
  }
  std::cout << m_className << "::Initialise: Bounding box:\n";
  std::cout << "    " << minx << " < x < " << maxx << "\n";
  std::cout << "    " << miny << " < y < " << maxy << "\n";
  std::cout << "    " << minz << " < z < " << maxz << "\n";

  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error reading " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("4 tet2 # type name") == std::string::npos);
  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error reading " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (!ends_with(line, "# number of elements"));
  nElements = readInt(line);
  elements.clear();
  std::cout << m_className << "::Initialise:\n"
            << "    Reading " << nElements << " elements.\n";
  std::getline(fmesh, line);
  // elements 6 & 7 are swapped due to differences in COMSOL and ANSYS
  // representation
  int perm[10] = {0, 1, 2, 3, 4, 5, 7, 6, 8, 9};
  for (int i = 0; i < nElements; ++i) {
    Element newElement;
    newElement.degenerate = false;
    for (int j = 0; j < 10; ++j) {
      fmesh >> newElement.emap[perm[j]];
    }
    elements.push_back(std::move(newElement));
  }

  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error reading " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("# Geometric entity indices") == std::string::npos);
  for (int i = 0; i < nElements; ++i) {
    int domain;
    fmesh >> domain;
    elements[i].matmap = domain2material.count(domain) ? domain2material[domain]
                                                       : m_nMaterials - 1;
  }
  fmesh.close();

  std::map<Node, std::vector<int>, nodeCmp> nodeIdx;
  for (int i = 0; i < nNodes; ++i) {
    nodeIdx[nodes[i]].push_back(i);
  }
  std::cout << "Map size: " << nodeIdx.size() << std::endl;

  std::ifstream ffield;
  ffield.open(field.c_str(), std::ios::in);
  if (ffield.fail()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Could not open potentials file " << field
              << " for reading.\n";
    return false;
  }
  do {
    if (!std::getline(ffield, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error reading " << field << ".\n";
      ffield.close();
      return false;
    }
  } while (line.substr(0, 81) !=
           "% x                       y                        z               "
           "         V (V)");
  {
    std::istringstream sline(line);
    std::string token;
    sline >> token;  // %
    sline >> token;  // x
    sline >> token;  // y
    sline >> token;  // z
    sline >> token;  // V
    sline >> token;  // (V)
    while (sline >> token) {
      std::cout << m_className << "::Initialise:\n";
      std::cout << "    Reading data for weighting field " << token << ".\n";
      nWeightingFields++;
      wfields.push_back(token);
      wfieldsOk.push_back(true);
      sline >> token;  // (V)
    }
  }
  for (int i = 0; i < nNodes; ++i) {
    Node tmp;
    ffield >> tmp.x >> tmp.y >> tmp.z >> tmp.v;
    tmp.x *= unit;
    tmp.y *= unit;
    tmp.z *= unit;
    for (int j = 0; j < nWeightingFields; ++j) {
      double w;
      ffield >> w;
      tmp.w.push_back(w);
    }
    int closest = -1;
    double closestDist = 1;
    const unsigned int nIdx = nodeIdx[tmp].size();
    // for (int j : nodeIdx[tmp]) {
    for (unsigned int k = 0; k < nIdx; ++k) {
      int j = nodeIdx[tmp][k];
      double dist = (tmp.x - nodes[j].x) * (tmp.x - nodes[j].x) +
                    (tmp.y - nodes[j].y) * (tmp.y - nodes[j].y) +
                    (tmp.z - nodes[j].z) * (tmp.z - nodes[j].z);
      if (dist < closestDist) {
        closestDist = dist;
        closest = j;
      }
    }
    if (closest == -1) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Could not match the node from field potentials file: "
                << tmp.x << " " << tmp.y << " " << tmp.z << "\n.";
      return false;
    }
    nodes[closest].v = tmp.v;
    nodes[closest].w = tmp.w;
  }

  m_ready = true;

  //  for (int i = 0; i < nNodes; ++i) {
  //    double ex, ey, ez, v;
  //    Medium* m;
  //    int status;
  //    ElectricField(nodes[i].x, nodes[i].y, nodes[i].z, ex, ey, ez, v, m,
  // status);
  //    std::cout << "Field at " << nodes[i].x << " " << nodes[i].y << " " <<
  // nodes[i].z << ": " << ex << " " << ey << " " << ez << " " << v << "\n";
  //  }

  // Establish the ranges.
  SetRange();
  UpdatePeriodicity();
  return true;
}

bool ComponentComsol::SetWeightingField(std::string field, std::string label) {
  double unit = 100.0;  // m;

  if (!m_ready) {
    std::cerr << m_className << "::SetWeightingField:\n";
    std::cerr << "    No valid field map is present.\n";
    std::cerr << "    Weighting field cannot be added.\n";
    return false;
  }

  // Open the voltage list.
  std::ifstream ffield;
  ffield.open(field.c_str(), std::ios::in);
  if (ffield.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open field potentials file " << field
              << " for reading.\n";
    return false;
  }

  // Check if a weighting field with the same label alm_ready exists.
  int iw = nWeightingFields;
  for (int i = nWeightingFields; i--;) {
    if (wfields[i] == label) {
      iw = i;
      break;
    }
  }
  if (iw == nWeightingFields) {
    ++nWeightingFields;
    wfields.resize(nWeightingFields);
    wfieldsOk.resize(nWeightingFields);
    for (int j = 0; j < nNodes; ++j) {
      nodes[j].w.resize(nWeightingFields);
    }
  } else {
    std::cout << m_className << "::SetWeightingField:\n";
    std::cout << "    Replacing existing weighting field " << label << ".\n";
  }
  wfields[iw] = label;
  wfieldsOk[iw] = false;
  std::map<Node, std::vector<int>, nodeCmp> nodeIdx;
  for (int i = 0; i < nNodes; ++i) {
    nodeIdx[nodes[i]].push_back(i);
  }
  std::cout << "Map size: " << nodeIdx.size() << std::endl;

  std::string line;
  do {
    std::getline(ffield, line);
  } while (line !=
           "% x                       y                        z               "
           "         V (V)");
  for (int i = 0; i < nNodes; ++i) {
    Node tmp;
    ffield >> tmp.x >> tmp.y >> tmp.z >> tmp.v;
    tmp.x *= unit;
    tmp.y *= unit;
    tmp.z *= unit;
    int closest = -1;
    double closestDist = 1;
    const unsigned int nIdx = nodeIdx[tmp].size();
    // for (int j : nodeIdx[tmp]) {
    for (unsigned int k = 0; k < nIdx; ++k) {
      int j = nodeIdx[tmp][k];
      double dist = (tmp.x - nodes[j].x) * (tmp.x - nodes[j].x) +
                    (tmp.y - nodes[j].y) * (tmp.y - nodes[j].y) +
                    (tmp.z - nodes[j].z) * (tmp.z - nodes[j].z);
      if (dist < closestDist) {
        closestDist = dist;
        closest = j;
      }
    }
    if (closest == -1) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Could not match the node from field potentials file: "
                << tmp.x << " " << tmp.y << " " << tmp.z << "\n.";
      return false;
    }
    nodes[closest].w[iw] = tmp.v;
  }

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
  m = NULL;

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
      std::cout << m_className << "::ElectricField:\n";
      std::cout << "    Point (" << x << ", " << y << ", " << z
                << " not in the mesh.\n";
    }
    status = -6;
    return;
  }

  const Element& element = elements[imap];
  if (m_debug) {
    PrintElement("ElectricField", x, y, z, t1, t2, t3, t4, element, 10);
  }
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  const Node& n4 = nodes[element.emap[4]];
  const Node& n5 = nodes[element.emap[5]];
  const Node& n6 = nodes[element.emap[6]];
  const Node& n7 = nodes[element.emap[7]];
  const Node& n8 = nodes[element.emap[8]];
  const Node& n9 = nodes[element.emap[9]];
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
  //  std::cout << "ef @(" << xin << ", " << yin << ", " << zin << ") = " <<
  // volt << "\n";

  // Drift medium?
  if (m_debug) {
    std::cout << m_className << "::ElectricField:\n";
    std::cout << "    Material " << element.matmap << ", drift flag "
              << materials[element.matmap].driftmedium << "\n";
  }
  m = materials[element.matmap].medium;
  status = -5;
  if (materials[element.matmap].driftmedium) {
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
  int iw = 0;
  bool found = false;
  for (int i = nWeightingFields; i--;) {
    if (wfields[i] == label) {
      iw = i;
      found = true;
      break;
    }
  }

  // Do not proceed if the requested weighting field does not exist.
  if (!found) return;
  // Check if the weighting field is properly initialised.
  if (!wfieldsOk[iw]) return;

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

  const Element& element = elements[imap];
  if (m_debug) {
    PrintElement("WeightingField", x, y, z, t1, t2, t3, t4, element, 10, iw);
  }
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  const Node& n4 = nodes[element.emap[4]];
  const Node& n5 = nodes[element.emap[5]];
  const Node& n6 = nodes[element.emap[6]];
  const Node& n7 = nodes[element.emap[7]];
  const Node& n8 = nodes[element.emap[8]];
  const Node& n9 = nodes[element.emap[9]];
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
  int iw = 0;
  bool found = false;
  for (int i = nWeightingFields; i--;) {
    if (wfields[i] == label) {
      iw = i;
      found = true;
      break;
    }
  }

  // Do not proceed if the requested weighting field does not exist.
  if (!found) return 0.;
  // Check if the weighting field is properly initialised.
  if (!wfieldsOk[iw]) return 0.;

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

  const Element& element = elements[imap];
  if (m_debug) {
    PrintElement("WeightingPotential", x, y, z, t1, t2, t3, t4, element, 10,
                 iw);
  }
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  const Node& n4 = nodes[element.emap[4]];
  const Node& n5 = nodes[element.emap[5]];
  const Node& n6 = nodes[element.emap[6]];
  const Node& n7 = nodes[element.emap[7]];
  const Node& n8 = nodes[element.emap[8]];
  const Node& n9 = nodes[element.emap[9]];
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
      std::cout << m_className << "::GetMedium:\n";
      std::cout << "    Point (" << x << ", " << y << ", " << z
                << ") not in the mesh.\n";
    }
    return nullptr;
  }
  const Element& element = elements[imap];
  if (element.matmap >= m_nMaterials) {
    if (m_debug) {
      std::cerr << m_className << "::GetMedium:\n";
      std::cerr << "    Point (" << x << ", " << y
                << ") has out of range material number " << imap << ".\n";
    }
    return nullptr;
  }

  if (m_debug) {
    PrintElement("GetMedium", x, y, z, t1, t2, t3, t4, element, 10);
  }

  return materials[element.matmap].medium;
}

double ComponentComsol::GetElementVolume(const unsigned int i) {
  if (i >= elements.size()) return 0.;
  const Element& element = elements[i];
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];

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
  if (i >= elements.size()) {
    dmin = dmax = 0.;
    return;
  }

  const Element& element = elements[i];
  const int np = 4;
  // Loop over all pairs of vertices.
  for (int j = 0; j < np - 1; ++j) {
    const Node& nj = nodes[element.emap[j]];
    for (int k = j + 1; k < np; ++k) {
      const Node& nk = nodes[element.emap[k]];
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
