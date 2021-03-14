#include "Garfield/ComponentFieldMap.hh"

#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "Garfield/FundamentalConstants.hh"

namespace Garfield {

ComponentFieldMap::ComponentFieldMap(const std::string& name)
    : Component(name) {}

ComponentFieldMap::~ComponentFieldMap() {}

void ComponentFieldMap::PrintMaterials() {
  // Do not proceed if not properly initialised.
  if (!m_ready) PrintNotReady("PrintMaterials");

  if (m_materials.empty()) {
    std::cerr << m_className << "::PrintMaterials:\n"
              << "    No materials are currently defined.\n";
    return;
  }

  const size_t nMaterials = m_materials.size();
  std::cout << m_className << "::PrintMaterials:\n"
            << "    Currently " << nMaterials << " materials are defined.\n"
            << "      Index Permittivity  Resistivity Notes\n";
  for (size_t i = 0; i < nMaterials; ++i) {
    printf("      %5zu %12g %12g", i, m_materials[i].eps, m_materials[i].ohm);
    if (m_materials[i].medium) {
      std::string name = m_materials[i].medium->GetName();
      std::cout << " " << name;
      if (m_materials[i].medium->IsDriftable()) std::cout << ", drift medium";
      if (m_materials[i].medium->IsIonisable()) std::cout << ", ionisable";
    }
    if (m_materials[i].driftmedium) {
      std::cout << " (drift medium)\n";
    } else {
      std::cout << "\n";
    }
  }
}

void ComponentFieldMap::DriftMedium(const unsigned int imat) {
  // Do not proceed if not properly initialised.
  if (!m_ready) PrintNotReady("DriftMedium");

  // Check value
  if (imat >= m_materials.size()) {
    std::cerr << m_className << "::DriftMedium: Index out of range.\n";
    return;
  }

  // Make drift medium
  m_materials[imat].driftmedium = true;
}

void ComponentFieldMap::NotDriftMedium(const unsigned int imat) {
  // Do not proceed if not properly initialised.
  if (!m_ready) PrintNotReady("NotDriftMedium");

  // Check value
  if (imat >= m_materials.size()) {
    std::cerr << m_className << "::NotDriftMedium: Index out of range.\n";
    return;
  }

  // Make drift medium
  m_materials[imat].driftmedium = false;
}

double ComponentFieldMap::GetPermittivity(const unsigned int imat) const {
  if (imat >= m_materials.size()) {
    std::cerr << m_className << "::GetPermittivity: Index out of range.\n";
    return -1.;
  }

  return m_materials[imat].eps;
}

double ComponentFieldMap::GetConductivity(const unsigned int imat) const {
  if (imat >= m_materials.size()) {
    std::cerr << m_className << "::GetConductivity: Index out of range.\n";
    return -1.;
  }

  return m_materials[imat].ohm;
}

void ComponentFieldMap::SetMedium(const unsigned int imat, Medium* m) {
  if (imat >= m_materials.size()) {
    std::cerr << m_className << "::SetMedium: Index out of range.\n";
    return;
  }

  if (!m) {
    std::cerr << m_className << "::SetMedium: Null pointer.\n";
    return;
  }

  if (m_debug) {
    std::cout << m_className << "::SetMedium:\n    Associated material " << imat
              << " with medium " << m->GetName() << ".\n";
  }

  m_materials[imat].medium = m;
}

Medium* ComponentFieldMap::GetMedium(const unsigned int imat) const {
  if (imat >= m_materials.size()) {
    std::cerr << m_className << "::GetMedium: Index out of range.\n";
    return nullptr;
  }

  return m_materials[imat].medium;
}

bool ComponentFieldMap::GetElement(const size_t i, double& vol,
                                   double& dmin, double& dmax) {
  if (i >= m_elements.size()) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }

  vol = GetElementVolume(i);
  GetAspectRatio(i, dmin, dmax);
  return true;
}

bool ComponentFieldMap::GetElement(const size_t i, size_t& mat, bool& drift,
                                   std::vector<size_t>& nodes) const {
  if (i >= m_elements.size()) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }
  const auto& element = m_elements[i];
  mat = element.matmap;
  drift = m_materials[mat].driftmedium;
  nodes.resize(4);
  for (size_t j = 0; j < 4; ++j) nodes[j] = element.emap[j];
  return true; 
}

bool ComponentFieldMap::GetNode(const size_t i, double& x, double& y,
                                double& z) const {
  if (i >= m_nodes.size()) {
    std::cerr << m_className << "::GetNode: Index out of range.\n";
    return false;
  }
  x = m_nodes[i].x;
  y = m_nodes[i].y;
  z = m_nodes[i].z;
  return true;
}

bool ComponentFieldMap::SetDefaultDriftMedium() {

  // Find lowest epsilon and set drift medium flags.
  const size_t nMaterials = m_materials.size(); 
  double epsmin = -1;
  size_t iepsmin = 0;
  for (size_t i = 0; i < nMaterials; ++i) {
    m_materials[i].driftmedium = false;
    if (m_materials[i].eps < 0) continue;
    // Check for eps == 0.
    if (m_materials[i].eps == 0) {
      std::cerr << m_className << "::SetDefaultDriftMedium:\n"
                << "    Material " << i << " has zero permittivity.\n";
      m_materials[i].eps = -1.;
    } else if (epsmin < 0. || epsmin > m_materials[i].eps) {
      epsmin = m_materials[i].eps;
      iepsmin = i;
    }
  }
  if (epsmin < 0.) {
    std::cerr << m_className << "::SetDefaultDriftMedium:\n"
              << "    Found no material with positive permittivity.\n";
    return false;
  }
  m_materials[iepsmin].driftmedium = true;
  return true;
}

int ComponentFieldMap::FindElement5(const double x, const double y,
                                    double const z, double& t1, double& t2,
                                    double& t3, double& t4, double jac[4][4],
                                    double& det) {
  // Check if bounding boxes of elements have been computed
  if (!m_cacheElemBoundingBoxes) {
    m_cacheElemBoundingBoxes = true;
  }

  // Tetra list in the block that contains the input 3D point.
  std::vector<int> tetList;
  if (m_useTetrahedralTree && m_octree) {
    tetList = m_octree->GetElementsInBlock(Vec3(x, y, z));
  }
  // Backup
  double jacbak[4][4], detbak = 1.;
  double t1bak = 0., t2bak = 0., t3bak = 0., t4bak = 0.;
  int imapbak = -1;

  // Initial values.
  t1 = t2 = t3 = t4 = 0;

  // Check previously used element
  if (m_lastElement > -1 && !m_checkMultipleElement) {
    const Element& element = m_elements[m_lastElement];
    if (element.degenerate) {
      if (Coordinates3(x, y, z, t1, t2, t3, t4, jac, det, element) == 0) {
        if (t1 >= 0 && t1 <= +1 && t2 >= 0 && t2 <= +1 && t3 >= 0 && t3 <= +1) {
          return m_lastElement;
        }
      }
    } else {
      if (Coordinates5(x, y, z, t1, t2, t3, t4, jac, det, element) == 0) {
        if (t1 >= -1 && t1 <= +1 && t2 >= -1 && t2 <= +1) return m_lastElement;
      }
    }
  }

  // Verify the count of volumes that contain the point.
  int nfound = 0;
  int imap = -1;

  // Number of elements to scan.
  // With tetra tree disabled, all elements are scanned.
  const int numElemToSearch =
      m_useTetrahedralTree ? tetList.size() : m_elements.size();
  for (int i = 0; i < numElemToSearch; ++i) {
    const int idxToElemList = m_useTetrahedralTree ? tetList[i] : i;
    const Element& element = m_elements[idxToElemList];
    if (x < element.bbMin[0] || x > element.bbMax[0] || 
        y < element.bbMin[1] || y > element.bbMax[1] || 
        z < element.bbMin[2] || z > element.bbMax[2])
      continue;
    if (element.degenerate) {
      // Degenerate element
      if (Coordinates3(x, y, z, t1, t2, t3, t4, jac, det, element) != 0) {
        continue;
      }
      if (t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1 || t3 < 0 || t3 > 1) continue;
      ++nfound;
      imap = idxToElemList;
      m_lastElement = idxToElemList;
      if (m_debug) {
        std::cout << m_className << "::FindElement5:\n";
        std::cout << "    Found matching degenerate element " << idxToElemList
                  << ".\n";
      }
      if (!m_checkMultipleElement) return idxToElemList;
      for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) jacbak[j][k] = jac[j][k];
      }
      detbak = det;
      t1bak = t1;
      t2bak = t2;
      t3bak = t3;
      t4bak = t4;
      imapbak = imap;
      if (m_debug) {
        PrintElement("FindElement5", x, y, z, t1, t2, t3, t4, element, 6);
      }
    } else {
      // Non-degenerate element
      if (Coordinates5(x, y, z, t1, t2, t3, t4, jac, det, element) != 0) {
        continue;
      }
      if (t1 < -1 || t1 > 1 || t2 < -1 || t2 > 1) continue;
      ++nfound;
      imap = idxToElemList;
      m_lastElement = idxToElemList;
      if (m_debug) {
        std::cout << m_className << "::FindElement5:\n";
        std::cout << "    Found matching non-degenerate element "
                  << idxToElemList << ".\n";
      }
      if (!m_checkMultipleElement) return idxToElemList;
      for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) jacbak[j][k] = jac[j][k];
      }
      detbak = det;
      t1bak = t1;
      t2bak = t2;
      t3bak = t3;
      t4bak = t4;
      imapbak = imap;
      if (m_debug) {
        PrintElement("FindElement5", x, y, z, t1, t2, t3, t4, element, 8);
      }
    }
  }

  // In checking mode, verify the tetrahedron/triangle count.
  if (m_checkMultipleElement) {
    if (nfound < 1) {
      if (m_debug) {
        std::cout << m_className << "::FindElement5:\n";
        std::cout << "    No element matching point (" << x << ", " << y
                  << ") found.\n";
      }
      m_lastElement = -1;
      return -1;
    }
    if (nfound > 1) {
      std::cout << m_className << "::FindElement5:\n";
      std::cout << "    Found " << nfound << " elements matching point (" << x
                << ", " << y << ").\n";
    }
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) jac[j][k] = jacbak[j][k];
    }
    det = detbak;
    t1 = t1bak;
    t2 = t2bak;
    t3 = t3bak;
    t4 = t4bak;
    imap = imapbak;
    m_lastElement = imap;
    return imap;
  }

  if (m_debug) {
    std::cout << m_className << "::FindElement5:\n";
    std::cout << "    No element matching point (" << x << ", " << y
              << ") found.\n";
  }
  return -1;
}

int ComponentFieldMap::FindElement13(const double x, const double y,
                                     const double z, double& t1, double& t2,
                                     double& t3, double& t4, double jac[4][4],
                                     double& det) {
  // Backup
  double jacbak[4][4];
  double detbak = 1.;
  double t1bak = 0., t2bak = 0., t3bak = 0., t4bak = 0.;
  int imapbak = -1;

  // Initial values.
  t1 = t2 = t3 = t4 = 0.;

  // Check previously used element
  if (m_lastElement > -1 && !m_checkMultipleElement) {
    const Element& element = m_elements[m_lastElement];
    if (Coordinates13(x, y, z, t1, t2, t3, t4, jac, det, element) == 0) {
      if (t1 >= 0 && t1 <= +1 && t2 >= 0 && t2 <= +1 && t3 >= 0 && t3 <= +1 &&
          t4 >= 0 && t4 <= +1) {
        return m_lastElement;
      }
    }
  }

  // Tetra list in the block that contains the input 3D point.
  std::vector<int> tetList;
  if (m_useTetrahedralTree && m_octree) {
    tetList = m_octree->GetElementsInBlock(Vec3(x, y, z));
  }
  // Number of elements to scan.
  // With tetra tree disabled, all elements are scanned.
  const int numElemToSearch =
      m_useTetrahedralTree ? tetList.size() : m_elements.size();
  // Verify the count of volumes that contain the point.
  int nfound = 0;
  int imap = -1;

  // Scan all elements
  for (int i = 0; i < numElemToSearch; i++) {
    const int idxToElemList = m_useTetrahedralTree ? tetList[i] : i;
    const Element& element = m_elements[idxToElemList];
    if (x < element.bbMin[0] || x > element.bbMax[0] || 
        y < element.bbMin[1] || y > element.bbMax[1] || 
        z < element.bbMin[2] || z > element.bbMax[2])
      continue;
    if (Coordinates13(x, y, z, t1, t2, t3, t4, jac, det, element) != 0) {
      continue;
    }
    if (t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1 || t3 < 0 || t3 > 1 || t4 < 0 ||
        t4 > 1) {
      continue;
    }
    ++nfound;
    imap = idxToElemList;
    m_lastElement = idxToElemList;
    if (m_debug) {
      std::cout << m_className << "::FindElement13:\n";
      std::cout << "    Found matching element " << i << ".\n";
    }
    if (!m_checkMultipleElement) return idxToElemList;
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) jacbak[j][k] = jac[j][k];
    }
    detbak = det;
    t1bak = t1;
    t2bak = t2;
    t3bak = t3;
    t4bak = t4;
    imapbak = imap;
    if (m_debug) {
      PrintElement("FindElement13", x, y, z, t1, t2, t3, t4, element, 10);
    }
  }

  // In checking mode, verify the tetrahedron/triangle count.
  if (m_checkMultipleElement) {
    if (nfound < 1) {
      if (m_debug) {
        std::cout << m_className << "::FindElement13:\n";
        std::cout << "    No element matching point (" << x << ", " << y << ", "
                  << z << ") found.\n";
      }
      m_lastElement = -1;
      return -1;
    }
    if (nfound > 1) {
      std::cerr << m_className << "::FindElement13:\n";
      std::cerr << "    Found << " << nfound << " elements matching point ("
                << x << ", " << y << ", " << z << ").\n";
    }
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) jac[j][k] = jacbak[j][k];
    }
    det = detbak;
    t1 = t1bak;
    t2 = t2bak;
    t3 = t3bak;
    t4 = t4bak;
    imap = imapbak;
    m_lastElement = imap;
    return imap;
  }

  if (m_debug) {
    std::cout << m_className << "::FindElement13:\n";
    std::cout << "    No element matching point (" << x << ", " << y << ", "
              << z << ") found.\n";
  }
  return -1;
}

int ComponentFieldMap::FindElementCube(const double x, const double y,
                                       const double z, double& t1, double& t2,
                                       double& t3, TMatrixD*& jac,
                                       std::vector<TMatrixD*>& dN) {
  int imap = -1;
  if (m_lastElement >= 0) {
    const Element& element = m_elements[m_lastElement];
    const Node& n3 = m_nodes[element.emap[3]];
    if (x >= n3.x && y >= n3.y && z >= n3.z) {
      const Node& n0 = m_nodes[element.emap[0]];
      const Node& n2 = m_nodes[element.emap[2]];
      const Node& n7 = m_nodes[element.emap[7]];
      if (x < n0.x && y < n2.y && z < n7.z) {
        imap = m_lastElement;
      }
    }
  }

  // Default element loop
  if (imap == -1) {
    const size_t nElements = m_elements.size();
    for (size_t i = 0; i < nElements; ++i) {
      const Element& element = m_elements[i];
      const Node& n3 = m_nodes[element.emap[3]];
      if (x < n3.x || y < n3.y || z < n3.z) continue;
      const Node& n0 = m_nodes[element.emap[0]];
      const Node& n2 = m_nodes[element.emap[2]];
      const Node& n7 = m_nodes[element.emap[7]];
      if (x < n0.x && y < n2.y && z < n7.z) {
        imap = i;
        break;
      }
    }
  }

  if (imap < 0) {
    if (m_debug) {
      std::cout << m_className << "::FindElementCube:\n";
      std::cout << "    Point (" << x << "," << y << "," << z
                << ") not in the mesh, it is background or PEC.\n";
      const Node& first0 = m_nodes[m_elements.front().emap[0]];
      const Node& first2 = m_nodes[m_elements.front().emap[2]];
      const Node& first3 = m_nodes[m_elements.front().emap[3]];
      const Node& first7 = m_nodes[m_elements.front().emap[7]];
      std::cout << "    First node (" << first3.x << "," << first3.y << ","
                << first3.z << ") in the mesh.\n";
      std::cout << "  dx= " << (first0.x - first3.x)
                << ", dy= " << (first2.y - first3.y)
                << ", dz= " << (first7.z - first3.z) << "\n";
      const Node& last0 = m_nodes[m_elements.back().emap[0]];
      const Node& last2 = m_nodes[m_elements.back().emap[2]];
      const Node& last3 = m_nodes[m_elements.back().emap[3]];
      const Node& last5 = m_nodes[m_elements.back().emap[5]];
      const Node& last7 = m_nodes[m_elements.back().emap[7]];
      std::cout << "    Last node (" << last5.x << "," << last5.y << ","
                << last5.z << ") in the mesh.\n";
      std::cout << "  dx= " << (last0.x - last3.x)
                << ", dy= " << (last2.y - last3.y)
                << ", dz= " << (last7.z - last3.z) << "\n";
    }
    return -1;
  }
  CoordinatesCube(x, y, z, t1, t2, t3, jac, dN, m_elements[imap]);
  if (m_debug) {
    PrintElement("FindElementCube", x, y, z, t1, t2, t3, 0., m_elements[imap],
                 8);
  }
  return imap;
}

void ComponentFieldMap::Jacobian3(const Element& element, const double u,
                                  const double v, const double w, double& det,
                                  double jac[4][4]) const {
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];

  // Shorthands.
  const double fouru = 4 * u;
  const double fourv = 4 * v;
  const double fourw = 4 * w;

  const double ax = (-1 + fourv) * n1.x + fouru * n3.x + fourw * n5.x;
  const double ay = (-1 + fourv) * n1.y + fouru * n3.y + fourw * n5.y;
  const double bx = (-1 + fourw) * n2.x + fouru * n4.x + fourv * n5.x;
  const double by = (-1 + fourw) * n2.y + fouru * n4.y + fourv * n5.y;
  const double cx = (-1 + fouru) * n0.x + fourv * n3.x + fourw * n4.x;
  const double cy = (-1 + fouru) * n0.y + fourv * n3.y + fourw * n4.y;
  // Determinant of the quadratic triangular Jacobian
  det = -(ax - bx) * cy - (cx - ax) * by + (cx - bx) * ay;

  // Terms of the quadratic triangular Jacobian
  jac[0][0] = ax * by - bx * ay;
  jac[0][1] = ay - by;
  jac[0][2] = bx - ax;
  jac[1][0] = bx * cy - cx * by;
  jac[1][1] = by - cy;
  jac[1][2] = cx - bx;
  jac[2][0] = -ax * cy + cx * ay;
  jac[2][1] = cy - ay;
  jac[2][2] = ax - cx;
}

void ComponentFieldMap::Jacobian5(const Element& element, const double u,
                                  const double v, double& det,
                                  double jac[4][4]) const {
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];
  const Node& n6 = m_nodes[element.emap[6]];
  const Node& n7 = m_nodes[element.emap[7]];
  const double u2 = u * u;
  const double v2 = v * v;
  const double twou = 2 * u;
  const double twov = 2 * v;
  const double two0x = 2 * n0.x;
  const double two0y = 2 * n0.y;
  const double two1x = 2 * n1.x;
  const double two1y = 2 * n1.y;
  const double two2x = 2 * n2.x;
  const double two2y = 2 * n2.y;
  const double two3x = 2 * n3.x;
  const double two3y = 2 * n3.y;
  const double two4x = 2 * n4.x;
  const double two4y = 2 * n4.y;
  const double two5x = 2 * n5.x;
  const double two5y = 2 * n5.y;
  const double two6x = 2 * n6.x;
  const double two6y = 2 * n6.y;
  const double two7x = 2 * n7.x;
  const double two7y = 2 * n7.y;
  // Determinant of the quadrilateral serendipity Jacobian
  det =
      (-twou * u2 *
           ((n2.x + n3.x - two6x) * (n0.y + n1.y - two4y) -
            (n0.x + n1.x - two4x) * (n2.y + n3.y - two6y)) +
       twov * v2 *
           (-((n0.x + n3.x - two7x) * (n1.y + n2.y - two5y)) +
            (n1.x + n2.x - two5x) * (n0.y + n3.y - two7y)) +
       2 * (-((n5.x - n7.x) * (n4.y - n6.y)) + (n4.x - n6.x) * (n5.y - n7.y)) +
       v * (-(n6.x * n0.y) - two7x * n0.y + n6.x * n1.y - two7x * n1.y -
            n6.x * n2.y - two7x * n2.y + n4.x * (n0.y - n1.y + n2.y - n3.y) +
            n6.x * n3.y - two7x * n3.y - n0.x * n4.y + n1.x * n4.y -
            n2.x * n4.y + n3.x * n4.y - two0x * n5.y - two1x * n5.y -
            two2x * n5.y - two3x * n5.y + 8 * n7.x * n5.y + n0.x * n6.y -
            n1.x * n6.y + n2.x * n6.y - n3.x * n6.y +
            two5x * (n0.y + n1.y + n2.y + n3.y - 4 * n7.y) +
            2 * (n0.x + n1.x + n2.x + n3.x) * n7.y) +
       v2 * (-(n4.x * n0.y) + two5x * n0.y + n6.x * n0.y + two7x * n0.y +
             n4.x * n1.y - two5x * n1.y - n6.x * n1.y - two7x * n1.y +
             n4.x * n2.y + two5x * n2.y - n6.x * n2.y + two7x * n2.y -
             n4.x * n3.y - two5x * n3.y + n6.x * n3.y - two7x * n3.y +
             two2x * (n1.y + n3.y) - n2.x * n4.y + two5x * n4.y - two7x * n4.y -
             two2x * n5.y - two4x * n5.y + two6x * n5.y + n2.x * n6.y -
             two5x * n6.y + two7x * n6.y +
             n0.x * (two1y + two3y + n4.y - two5y - n6.y - two7y) -
             2 * (n2.x - n4.x + n6.x) * n7.y +
             n3.x * (-two0y - two2y + n4.y + two5y - n6.y + two7y) +
             n1.x * (-two0y - two2y - n4.y + two5y + n6.y + two7y)) +
       u * (n5.x * n0.y - two6x * n0.y - n7.x * n0.y - n5.x * n1.y -
            two6x * n1.y + n7.x * n1.y + n5.x * n2.y - two6x * n2.y -
            n7.x * n2.y - n5.x * n3.y - two6x * n3.y + n7.x * n3.y -
            two1x * n4.y - two2x * n4.y - two3x * n4.y + 8 * n6.x * n4.y +
            n1.x * n5.y - n2.x * n5.y + n3.x * n5.y +
            two4x * (n0.y + n1.y + n2.y + n3.y - 4 * n6.y) + two1x * n6.y +
            two2x * n6.y + two3x * n6.y - (n1.x - n2.x + n3.x) * n7.y +
            n0.x * (-two4y - n5.y + two6y + n7.y) +
            v2 * (4 * n4.x * n0.y - 3 * n5.x * n0.y - 4 * n6.x * n0.y -
                  5 * n7.x * n0.y + 4 * n4.x * n1.y - 5 * n5.x * n1.y -
                  4 * n6.x * n1.y - 3 * n7.x * n1.y + 4 * n4.x * n2.y +
                  5 * n5.x * n2.y - 4 * n6.x * n2.y + 3 * n7.x * n2.y +
                  4 * n4.x * n3.y + 3 * n5.x * n3.y - 4 * n6.x * n3.y +
                  5 * n7.x * n3.y + 8 * n5.x * n4.y + 8 * n7.x * n4.y -
                  8 * n4.x * n5.y + 8 * n6.x * n5.y - 8 * n5.x * n6.y -
                  8 * n7.x * n6.y +
                  n3.x * (5 * n0.y + 3 * n1.y - 4 * n4.y - 3 * n5.y + 4 * n6.y -
                          5 * n7.y) +
                  n2.x * (3 * n0.y + 5 * n1.y - 4 * n4.y - 5 * n5.y + 4 * n6.y -
                          3 * n7.y) -
                  8 * n4.x * n7.y + 8 * n6.x * n7.y +
                  n1.x * (-5 * n2.y - 3 * n3.y - 4 * n4.y + 5 * n5.y +
                          4 * n6.y + 3 * n7.y) +
                  n0.x * (-3 * n2.y - 5 * n3.y - 4 * n4.y + 3 * n5.y +
                          4 * n6.y + 5 * n7.y)) -
            twov * (n6.x * n0.y - 3 * n7.x * n0.y + n6.x * n1.y - n7.x * n1.y +
                    3 * n6.x * n2.y - n7.x * n2.y + 3 * n6.x * n3.y -
                    3 * n7.x * n3.y - 3 * n0.x * n4.y - 3 * n1.x * n4.y -
                    n2.x * n4.y - n3.x * n4.y + 4 * n7.x * n4.y + n0.x * n5.y +
                    3 * n1.x * n5.y + 3 * n2.x * n5.y + n3.x * n5.y -
                    4 * n6.x * n5.y - n0.x * n6.y - n1.x * n6.y -
                    3 * n2.x * n6.y - 3 * n3.x * n6.y + 4 * n7.x * n6.y -
                    n5.x * (n0.y + 3 * n1.y + 3 * n2.y + n3.y -
                            4 * (n4.y + n6.y)) +
                    (3 * n0.x + n1.x + n2.x + 3 * n3.x - 4 * n6.x) * n7.y +
                    n4.x * (3 * n0.y + 3 * n1.y + n2.y + n3.y -
                            4 * (n5.y + n7.y)))) +
       u2 *
           (two3x * n0.y - two4x * n0.y - n5.x * n0.y - two6x * n0.y +
            n7.x * n0.y - two0x * n1.y + two4x * n1.y - n5.x * n1.y +
            two6x * n1.y + n7.x * n1.y + two3x * n2.y - two4x * n2.y +
            n5.x * n2.y - two6x * n2.y - n7.x * n2.y + two4x * n3.y +
            n5.x * n3.y + two6x * n3.y - n7.x * n3.y - two3x * n4.y +
            two5x * n4.y - two7x * n4.y - n3.x * n5.y - two4x * n5.y +
            two6x * n5.y - two3x * n6.y - two5x * n6.y + two7x * n6.y +
            n0.x * (-two3y + two4y + n5.y + two6y - n7.y) +
            (n3.x + two4x - two6x) * n7.y +
            n2.x * (-two1y - two3y + two4y - n5.y + two6y + n7.y) -
            3 * v2 *
                (n5.x * n0.y - n6.x * n0.y - n7.x * n0.y + n5.x * n1.y +
                 n6.x * n1.y - n7.x * n1.y - n5.x * n2.y + n6.x * n2.y +
                 n7.x * n2.y - n5.x * n3.y - n6.x * n3.y + n7.x * n3.y -
                 two5x * n4.y + two7x * n4.y - two6x * n5.y + two5x * n6.y -
                 two7x * n6.y +
                 n4.x * (n0.y - n1.y - n2.y + n3.y + two5y - two7y) +
                 n3.x * (n0.y - n2.y - n4.y + n5.y + n6.y - n7.y) +
                 two6x * n7.y +
                 (n0.x - n2.x) * (n1.y - n3.y - n4.y - n5.y + n6.y + n7.y)) +
            v * (4 * n5.x * n0.y + 3 * n6.x * n0.y - 4 * n7.x * n0.y +
                 4 * n5.x * n1.y - 3 * n6.x * n1.y - 4 * n7.x * n1.y +
                 4 * n5.x * n2.y - 5 * n6.x * n2.y - 4 * n7.x * n2.y +
                 4 * n5.x * n3.y + 5 * n6.x * n3.y - 4 * n7.x * n3.y -
                 8 * n5.x * n4.y + 8 * n7.x * n4.y + 8 * n6.x * n5.y -
                 8 * n5.x * n6.y + 8 * n7.x * n6.y +
                 n4.x * (5 * n0.y - 5 * n1.y - 3 * n2.y + 3 * n3.y + 8 * n5.y -
                         8 * n7.y) -
                 8 * n6.x * n7.y +
                 n3.x * (3 * n1.y + 5 * n2.y - 3 * n4.y - 4 * n5.y - 5 * n6.y +
                         4 * n7.y) +
                 n0.x * (5 * n1.y + 3 * n2.y - 5 * n4.y - 4 * n5.y - 3 * n6.y +
                         4 * n7.y) +
                 n2.x * (-3 * n0.y - 5 * n3.y + 3 * n4.y - 4 * n5.y + 5 * n6.y +
                         4 * n7.y)) +
            n1.x * ((-1 + v) * (-2 + 3 * v) * n0.y + two2y - two4y + n5.y -
                    two6y - n7.y +
                    v * (-3 * n3.y + 5 * n4.y - 4 * n5.y + 3 * n6.y + 4 * n7.y -
                         3 * v * (n2.y + n4.y - n5.y - n6.y + n7.y))))) /
      8;
  // Jacobian terms
  jac[0][0] =
      (u2 * (-n0.y - n1.y + n2.y + n3.y + two4y - two6y) +
       2 * (-n4.y + n6.y + v * (n0.y + n1.y + n2.y + n3.y - two5y - two7y)) +
       u * (n0.y - twov * n0.y - n1.y + twov * n1.y + n2.y + twov * n2.y -
            n3.y - twov * n3.y - twov * two5y + twov * two7y)) /
      4;
  jac[0][1] =
      (u2 * (n0.x + n1.x - n2.x - n3.x - two4x + two6x) -
       2 * (-n4.x + n6.x + v * (n0.x + n1.x + n2.x + n3.x - two5x - two7x)) +
       u * ((-1 + twov) * n0.x + n1.x - twov * n1.x - n2.x - twov * n2.x +
            n3.x + twov * n3.x + twov * two5x - twov * two7x)) /
      4;
  jac[1][0] =
      (v * (-n0.y + n1.y - n2.y + n3.y) - two5y +
       twou * ((-1 + v) * n0.y + (-1 + v) * n1.y - n2.y - v * n2.y - n3.y -
               v * n3.y + two4y - twov * n4.y + two6y + twov * n6.y) +
       v2 * (n0.y - n1.y - n2.y + n3.y + two5y - two7y) + two7y) /
      4;
  jac[1][1] =
      (v * (n0.x - n1.x + n2.x - n3.x) +
       twou * (n0.x - v * n0.x + n1.x - v * n1.x + n2.x + v * n2.x + n3.x +
               v * n3.x - two4x + twov * n4.x - two6x - twov * n6.x) +
       two5x - two7x + v2 * (-n0.x + n1.x + n2.x - n3.x - two5x + two7x)) /
      4;
}

void ComponentFieldMap::Jacobian13(const Element& element, const double t,
                                   const double u, const double v,
                                   const double w, double& det,
                                   double jac[4][4]) const {
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

  // Shorthands.
  const double fourt = 4 * t;
  const double fouru = 4 * u;
  const double fourv = 4 * v;
  const double fourw = 4 * w;

  const double x0 = (-1 + fourt) * n0.x;
  const double y0 = (-1 + fourt) * n0.y;
  const double z0 = (-1 + fourt) * n0.z;
  const double xp = fouru * n4.x + fourv * n5.x + fourw * n6.x;
  const double yp = fouru * n4.y + fourv * n5.y + fourw * n6.y;
  const double zp = fouru * n4.z + fourv * n5.z + fourw * n6.z;

  const double tx = x0 + xp;
  const double ty = y0 + yp;
  const double tz = z0 + zp;

  const double x1 = (-1 + fouru) * n1.x;
  const double y1 = (-1 + fouru) * n1.y;
  const double z1 = (-1 + fouru) * n1.z;
  const double xq = fourt * n4.x + fourv * n7.x + fourw * n8.x;
  const double yq = fourt * n4.y + fourv * n7.y + fourw * n8.y;
  const double zq = fourt * n4.z + fourv * n7.z + fourw * n8.z;

  const double ux = x1 + xq;
  const double uy = y1 + yq;
  const double uz = z1 + zq;

  const double x2 = (-1 + fourv) * n2.x;
  const double y2 = (-1 + fourv) * n2.y;
  const double z2 = (-1 + fourv) * n2.z;
  const double xr = fourt * n5.x + fouru * n7.x + fourw * n9.x;
  const double yr = fourt * n5.y + fouru * n7.y + fourw * n9.y;
  const double zr = fourt * n5.z + fouru * n7.z + fourw * n9.z;

  const double vx = x2 + xr;
  const double vy = y2 + yr;
  const double vz = z2 + zr;

  const double x3 = (-1 + fourw) * n3.x;
  const double y3 = (-1 + fourw) * n3.y;
  const double z3 = (-1 + fourw) * n3.z;
  const double xs = fourt * n6.x + fouru * n8.x + fourv * n9.x;
  const double ys = fourt * n6.y + fouru * n8.y + fourv * n9.y;
  const double zs = fourt * n6.z + fouru * n8.z + fourv * n9.z;

  const double wx = x3 + xs;
  const double wy = y3 + ys;
  const double wz = z3 + zs;

  const double ax = x1 - x3 + xq - xs;
  const double ay = y1 - y3 + yq - ys;

  const double bx = x1 - x2 + xq - xr;
  const double by = y1 - y2 + yq - yr;

  const double cx = x2 - x3 + xr - xs;
  const double cy = y2 - y3 + yr - ys;

  const double dx = x0 - x3 + xp - xs;
  const double dy = y0 - y3 + yp - ys;

  const double ex = x0 - x2 + xp - xr;
  const double ey = y0 - y2 + yp - yr;

  const double fx = x0 - x1 + xp - xq;
  const double fy = y0 - y1 + yp - yq;

  // Determinant of the quadrilateral serendipity Jacobian
  det = (-ax * vy + bx * wy + cx * uy) * tz -
        (-ax * ty - fx * wy + dx * uy) * vz +
        (-bx * ty - fx * vy + ex * uy) * wz +
        (-cx * ty + dx * vy - ex * wy) * uz;

  const double tu = tx * uy - ux * ty;
  const double tv = tx * vy - vx * ty;
  const double tw = tx * wy - wx * ty;
  const double uv = ux * vy - vx * uy;
  const double uw = ux * wy - wx * uy;
  const double vw = vx * wy - wx * vy;

  jac[0][0] = -uw * vz + uv * wz + vw * uz;
  jac[1][0] = -vw * tz + tw * vz - tv * wz;
  jac[2][0] =  uw * tz + tu * wz - tw * uz;
  jac[3][0] = -uv * tz - tu * vz + tv * uz;

  jac[0][1] = -ay * vz + by * wz + cy * uz;
  jac[0][2] =  ax * vz - bx * wz - cx * uz;
  jac[0][3] = -ax * vy + bx * wy + cx * uy;

  jac[1][1] = -cy * tz + dy * vz - ey * wz;
  jac[1][2] =  cx * tz - dx * vz + ex * wz;
  jac[1][3] = -cx * ty + dx * vy - ex * wy;

  jac[2][1] =  ay * tz + fy * wz - dy * uz;
  jac[2][2] = -ax * tz - fx * wz + dx * uz;
  jac[2][3] =  ax * ty + fx * wy - dx * uy;

  jac[3][1] = -by * tz - fy * vz + ey * uz;
  jac[3][2] =  bx * tz + fx * vz - ex * uz;
  jac[3][3] = -bx * ty - fx * vy + ex * uy;
}

void ComponentFieldMap::JacobianCube(const Element& element, const double t1,
                                     const double t2, const double t3,
                                     TMatrixD*& jac,
                                     std::vector<TMatrixD*>& dN) const {
  if (!jac) {
    std::cerr << m_className << "::JacobianCube:\n";
    std::cerr << "    Pointer to Jacobian matrix is empty!\n";
    return;
  }
  dN.clear();

  // Here the partial derivatives of the 8 shaping functions are calculated
  double N1[3] = {-1 * (1 - t2) * (1 - t3), (1 - t1) * -1 * (1 - t3),
                  (1 - t1) * (1 - t2) * -1};
  double N2[3] = {+1 * (1 - t2) * (1 - t3), (1 + t1) * -1 * (1 - t3),
                  (1 + t1) * (1 - t2) * -1};
  double N3[3] = {+1 * (1 + t2) * (1 - t3), (1 + t1) * +1 * (1 - t3),
                  (1 + t1) * (1 + t2) * -1};
  double N4[3] = {-1 * (1 + t2) * (1 - t3), (1 - t1) * +1 * (1 - t3),
                  (1 - t1) * (1 + t2) * -1};
  double N5[3] = {-1 * (1 - t2) * (1 + t3), (1 - t1) * -1 * (1 + t3),
                  (1 - t1) * (1 - t2) * +1};
  double N6[3] = {+1 * (1 - t2) * (1 + t3), (1 + t1) * -1 * (1 + t3),
                  (1 + t1) * (1 - t2) * +1};
  double N7[3] = {+1 * (1 + t2) * (1 + t3), (1 + t1) * +1 * (1 + t3),
                  (1 + t1) * (1 + t2) * +1};
  double N8[3] = {-1 * (1 + t2) * (1 + t3), (1 - t1) * +1 * (1 + t3),
                  (1 - t1) * (1 + t2) * +1};
  // Partial derivatives are stored in dN
  TMatrixD* m_N1 = new TMatrixD(3, 1, N1);
  *m_N1 = (1. / 8. * (*m_N1));
  dN.push_back(m_N1);
  TMatrixD* m_N2 = new TMatrixD(3, 1, N2);
  *m_N2 = (1. / 8. * (*m_N2));
  dN.push_back(m_N2);
  TMatrixD* m_N3 = new TMatrixD(3, 1, N3);
  *m_N3 = (1. / 8. * (*m_N3));
  dN.push_back(m_N3);
  TMatrixD* m_N4 = new TMatrixD(3, 1, N4);
  *m_N4 = (1. / 8. * (*m_N4));
  dN.push_back(m_N4);
  TMatrixD* m_N5 = new TMatrixD(3, 1, N5);
  *m_N5 = (1. / 8. * (*m_N5));
  dN.push_back(m_N5);
  TMatrixD* m_N6 = new TMatrixD(3, 1, N6);
  *m_N6 = (1. / 8. * (*m_N6));
  dN.push_back(m_N6);
  TMatrixD* m_N7 = new TMatrixD(3, 1, N7);
  *m_N7 = (1. / 8. * (*m_N7));
  dN.push_back(m_N7);
  TMatrixD* m_N8 = new TMatrixD(3, 1, N8);
  *m_N8 = (1. / 8. * (*m_N8));
  dN.push_back(m_N8);
  // Calculation of the jacobian using dN
  for (int j = 0; j < 8; ++j) {
    const Node& node = m_nodes[element.emap[j]];
    (*jac)(0, 0) += node.x * ((*dN.at(j))(0, 0));
    (*jac)(0, 1) += node.y * ((*dN.at(j))(0, 0));
    (*jac)(0, 2) += node.z * ((*dN.at(j))(0, 0));
    (*jac)(1, 0) += node.x * ((*dN.at(j))(1, 0));
    (*jac)(1, 1) += node.y * ((*dN.at(j))(1, 0));
    (*jac)(1, 2) += node.z * ((*dN.at(j))(1, 0));
    (*jac)(2, 0) += node.x * ((*dN.at(j))(2, 0));
    (*jac)(2, 1) += node.y * ((*dN.at(j))(2, 0));
    (*jac)(2, 2) += node.z * ((*dN.at(j))(2, 0));
  }

  // compute determinant
  if (m_debug) {
    std::cout << m_className << "::JacobianCube:" << std::endl;
    std::cout << "   Det.: " << jac->Determinant() << std::endl;
    std::cout << "   Jacobian matrix.: " << std::endl;
    jac->Print("%11.10g");
    std::cout << "   Hexahedral coordinates (t, u, v) = (" << t1 << "," << t2
              << "," << t3 << ")" << std::endl;
    std::cout << "   Node xyzV" << std::endl;
    for (int j = 0; j < 8; ++j) {
      const Node& node = m_nodes[element.emap[j]];
      std::cout << "         " << element.emap[j] << "          " << node.x
                << "         " << node.y << "         " << node.z << "         "
                << node.v << std::endl;
    }
  }
}

int ComponentFieldMap::Coordinates3(const double x, const double y,
                                    const double z, double& t1, double& t2,
                                    double& t3, double& t4, double jac[4][4],
                                    double& det, const Element& element) const {
  if (m_debug) {
    std::cout << m_className << "::Coordinates3:\n";
    std::cout << "   Point (" << x << ", " << y << ", " << z << ")\n";
  }

  // Provisional values
  t1 = t2 = t3 = t4 = 0;

  // Make a first order approximation, using the linear triangle.
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const double d1 =
      (n0.x - n1.x) * (n2.y - n1.y) - (n2.x - n1.x) * (n0.y - n1.y);
  const double d2 =
      (n1.x - n2.x) * (n0.y - n2.y) - (n0.x - n2.x) * (n1.y - n2.y);
  const double d3 =
      (n2.x - n0.x) * (n1.y - n0.y) - (n1.x - n0.x) * (n2.y - n0.y);
  if (d1 == 0 || d2 == 0 || d3 == 0) {
    std::cerr << m_className << "::Coordinates3:\n";
    std::cerr << "    Calculation of linear coordinates failed; abandoned.\n";
    return 1;
  }
  t1 = ((x - n1.x) * (n2.y - n1.y) - (y - n1.y) * (n2.x - n1.x)) / d1;
  t2 = ((x - n2.x) * (n0.y - n2.y) - (y - n2.y) * (n0.x - n2.x)) / d2;
  t3 = ((x - n0.x) * (n1.y - n0.y) - (y - n0.y) * (n1.x - n0.x)) / d3;

  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];

  // Start iterative refinement.
  double td1 = t1, td2 = t2, td3 = t3;
  bool converged = false;
  for (int iter = 0; iter < 10; iter++) {
    if (m_debug) {
      std::cout << m_className << "::Coordinates3:\n";
      std::cout << "    Iteration " << iter << ":     (u, v, w) = (" << td1
                << ", " << td2 << ", " << td3 << "), sum = " << td1 + td2 + td3
                << "\n";
    }
    // Evaluate the shape functions.
    const double f0 = td1 * (2 * td1 - 1);
    const double f1 = td2 * (2 * td2 - 1);
    const double f2 = td3 * (2 * td3 - 1);
    const double f3 = 4 * td1 * td2;
    const double f4 = 4 * td1 * td3;
    const double f5 = 4 * td2 * td3;
    // Re-compute the (x,y,z) position for this coordinate.
    const double xr =
        n0.x * f0 + n1.x * f1 + n2.x * f2 + n3.x * f3 + n4.x * f4 + n5.x * f5;
    const double yr =
        n0.y * f0 + n1.y * f1 + n2.y * f2 + n3.y * f3 + n4.y * f4 + n5.y * f5;
    const double sr = td1 + td2 + td3;
    // Compute the Jacobian.
    Jacobian3(element, td1, td2, td3, det, jac);
    // Compute the difference vector.
    const double diff[3] = {1 - sr, x - xr, y - yr};
    // Update the estimate.
    const double invdet = 1. / det;
    double corr[3] = {0., 0., 0.};
    for (size_t l = 0; l < 3; l++) {
      for (size_t k = 0; k < 3; k++) {
        corr[l] += jac[l][k] * diff[k];
      }
      corr[l] *= invdet;
    }
    // Debugging
    if (m_debug) {
      std::cout << m_className << "::Coordinates3:\n";
      std::cout << "    Difference vector:  (1, x, y)  = (" << diff[0] << ", "
                << diff[1] << ", " << diff[2] << ").\n";
      std::cout << "    Correction vector:  (u, v, w) = (" << corr[0] << ", "
                << corr[1] << ", " << corr[2] << ").\n";
    }
    // Update the vector.
    td1 += corr[0];
    td2 += corr[1];
    td3 += corr[2];
    // Check for convergence.
    constexpr double tol = 1.e-5;
    if (fabs(corr[0]) < tol && fabs(corr[1]) < tol && fabs(corr[2]) < tol) {
      if (m_debug) {
        std::cout << m_className << "::Coordinates3: Convergence reached.";
      }
      converged = true;
      break;
    }
  }
  // No convergence reached
  if (!converged) {
    const double xmin = std::min({n0.x, n1.x, n2.x});
    const double xmax = std::max({n0.x, n1.x, n2.x});
    const double ymin = std::min({n0.y, n1.y, n2.y});
    const double ymax = std::max({n0.y, n1.y, n2.y});
    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
      if (m_printConvergenceWarnings) {
        std::cout << m_className << "::Coordinates3:\n"
                  << "    No convergence achieved "
                  << "when refining internal isoparametric coordinates\n"
                  << "    at position (" << x << ", " << y << ").\n";
      }
      t1 = t2 = t3 = t4 = 0;
      return 1;
    }
  }

  // Convergence reached.
  t1 = td1;
  t2 = td2;
  t3 = td3;
  t4 = 0;
  if (m_debug) {
    std::cout << m_className << "::Coordinates3:\n";
    std::cout << "    Convergence reached at (t1, t2, t3) = (" << t1 << ", "
              << t2 << ", " << t3 << ").\n";
    // For debugging purposes, show position
    const double f0 = td1 * (2 * td1 - 1);
    const double f1 = td2 * (2 * td2 - 1);
    const double f2 = td3 * (2 * td3 - 1);
    const double f3 = 4 * td1 * td2;
    const double f4 = 4 * td1 * td3;
    const double f5 = 4 * td2 * td3;
    const double xr =
        n0.x * f0 + n1.x * f1 + n2.x * f2 + n3.x * f3 + n4.x * f4 + n5.x * f5;
    const double yr =
        n0.y * f0 + n1.y * f1 + n2.y * f2 + n3.y * f3 + n4.y * f4 + n5.y * f5;
    const double sr = td1 + td2 + td3;
    std::cout << m_className << "::Coordinates3:\n";
    std::cout << "    Position requested:     (" << x << ", " << y << ")\n";
    std::cout << "    Reconstructed:          (" << xr << ", " << yr << ")\n";
    std::cout << "    Difference:             (" << x - xr << ", " << y - yr
              << ")\n";
    std::cout << "    Checksum - 1:           " << sr - 1 << "\n";
  }

  // Success
  return 0;
}

int ComponentFieldMap::Coordinates4(const double x, const double y,
                                    const double z, double& t1, double& t2,
                                    double& t3, double& t4, double& det,
                                    const Element& element) const {
  // Debugging
  if (m_debug) {
    std::cout << m_className << "::Coordinates4:\n";
    std::cout << "   Point (" << x << ", " << y << ", " << z << ")\n";
  }

  // Failure flag
  int ifail = 1;

  // Provisional values
  t1 = t2 = t3 = t4 = 0.;

  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  // Compute determinant.
  const double dd = -(n0.x * n1.y) + n3.x * n2.y - n2.x * n3.y +
                    x * (-n0.y + n1.y - n2.y + n3.y) + n1.x * (n0.y - y) +
                    (n0.x + n2.x - n3.x) * y;
  det = -(-((n0.x - n3.x) * (n1.y - n2.y)) + (n1.x - n2.x) * (n0.y - n3.y)) *
            (2 * x * (-n0.y + n1.y + n2.y - n3.y) -
             (n0.x + n3.x) * (n1.y + n2.y - 2 * y) +
             n1.x * (n0.y + n3.y - 2 * y) + n2.x * (n0.y + n3.y - 2 * y)) +
        dd * dd;

  // Check that the determinant is non-negative
  // (this can happen if the point is out of range).
  if (det < 0) {
    if (m_debug) {
      std::cerr << m_className << "::Coordinates4:\n"
                << "    No solution found for isoparametric coordinates\n"
                << "    because the determinant " << det << " is < 0.\n";
    }
    return ifail;
  }

  // Vector products for evaluation of T1.
  double prod = ((n2.x - n3.x) * (n0.y - n1.y) - (n0.x - n1.x) * (n2.y - n3.y));
  if (prod * prod >
      1.0e-12 *
          ((n0.x - n1.x) * (n0.x - n1.x) + (n0.y - n1.y) * (n0.y - n1.y)) *
          ((n2.x - n3.x) * (n2.x - n3.x) + (n2.y - n3.y) * (n2.y - n3.y))) {
    t1 = (-(n3.x * n0.y) + x * n0.y + n2.x * n1.y - x * n1.y - n1.x * n2.y +
          x * n2.y + n0.x * n3.y - x * n3.y - n0.x * y + n1.x * y - n2.x * y +
          n3.x * y + sqrt(det)) /
         prod;
  } else {
    double xp = n0.y - n1.y;
    double yp = n1.x - n0.x;
    double dn = sqrt(xp * xp + yp * yp);
    if (dn <= 0) {
      std::cerr << m_className << "::Coordinates4:\n"
                << "    Element appears to be degenerate in the 1 - 2 axis.\n";
      return ifail;
    }
    xp = xp / dn;
    yp = yp / dn;
    double dpoint = xp * (x - n0.x) + yp * (y - n0.y);
    double dbox = xp * (n3.x - n0.x) + yp * (n3.y - n0.y);
    if (dbox == 0) {
      std::cerr << m_className << "::Coordinates4:\n"
                << "    Element appears to be degenerate in the 1 - 3 axis.\n";
      return ifail;
    }
    double t = -1 + 2 * dpoint / dbox;
    double xt1 = n0.x + 0.5 * (t + 1) * (n3.x - n0.x);
    double yt1 = n0.y + 0.5 * (t + 1) * (n3.y - n0.y);
    double xt2 = n1.x + 0.5 * (t + 1) * (n2.x - n1.x);
    double yt2 = n1.y + 0.5 * (t + 1) * (n2.y - n1.y);
    dn = (xt1 - xt2) * (xt1 - xt2) + (yt1 - yt2) * (yt1 - yt2);
    if (dn <= 0) {
      std::cout << m_className << "::Coordinates4:\n";
      std::cout
          << "    Coordinate requested at convergence point of element.\n";
      return ifail;
    }
    t1 = -1 + 2 * ((x - xt1) * (xt2 - xt1) + (y - yt1) * (yt2 - yt1)) / dn;
  }

  // Vector products for evaluation of T2.
  prod = ((n0.x - n3.x) * (n1.y - n2.y) - (n1.x - n2.x) * (n0.y - n3.y));
  if (prod * prod >
      1.0e-12 *
          ((n0.x - n3.x) * (n0.x - n3.x) + (n0.y - n3.y) * (n0.y - n3.y)) *
          ((n1.x - n2.x) * (n1.x - n2.x) + (n1.y - n2.y) * (n1.y - n2.y))) {
    t2 = (-(n1.x * n0.y) + x * n0.y + n0.x * n1.y - x * n1.y - n3.x * n2.y +
          x * n2.y + n2.x * n3.y - x * n3.y - n0.x * y + n1.x * y - n2.x * y +
          n3.x * y - sqrt(det)) /
         prod;
  } else {
    double xp = n0.y - n3.y;
    double yp = n3.x - n0.x;
    double dn = sqrt(xp * xp + yp * yp);
    if (dn <= 0) {
      std::cerr << m_className << "Coordinates4:\n"
                << "    Element appears to be degenerate in the 1 - 4 axis.\n";
      return ifail;
    }
    xp = xp / dn;
    yp = yp / dn;
    double dpoint = xp * (x - n0.x) + yp * (y - n0.y);
    double dbox = xp * (n1.x - n0.x) + yp * (n1.y - n0.y);
    if (dbox == 0) {
      std::cerr << m_className << "::Coordinates4:\n"
                << "    Element appears to be degenerate in the 1 - 2 axis.\n";
      return ifail;
    }
    double t = -1 + 2 * dpoint / dbox;
    double xt1 = n0.x + 0.5 * (t + 1) * (n1.x - n0.x);
    double yt1 = n0.y + 0.5 * (t + 1) * (n1.y - n0.y);
    double xt2 = n3.x + 0.5 * (t + 1) * (n2.x - n3.x);
    double yt2 = n3.y + 0.5 * (t + 1) * (n2.y - n3.y);
    dn = (xt1 - xt2) * (xt1 - xt2) + (yt1 - yt2) * (yt1 - yt2);
    if (dn <= 0) {
      std::cout
          << m_className << "::Coordinates4:\n"
          << "    Coordinate requested at convergence point of element.\n";
      return ifail;
    }
    t2 = -1 + 2 * ((x - xt1) * (xt2 - xt1) + (y - yt1) * (yt2 - yt1)) / dn;
  }
  if (m_debug) {
    std::cout << m_className << "::Coordinates4:\n";
    std::cout << "    Isoparametric (u, v):   (" << t1 << ", " << t2 << ").\n";
    // Re-compute the (x,y,z) position for this coordinate.
    const double f0 = (1 - t1) * (1 - t2) * 0.25;
    const double f1 = (1 + t1) * (1 - t2) * 0.25;
    const double f2 = (1 + t1) * (1 + t2) * 0.25;
    const double f3 = (1 - t1) * (1 + t2) * 0.25;
    const double xr = n0.x * f0 + n1.x * f1 + n2.x * f2 + n3.x * f3;
    const double yr = n0.y * f0 + n1.y * f1 + n2.y * f2 + n3.y * f3;
    std::cout << m_className << "::Coordinates4: \n";
    std::cout << "    Position requested:     (" << x << ", " << y << ")\n";
    std::cout << "    Reconstructed:          (" << xr << ", " << yr << ")\n";
    std::cout << "    Difference:             (" << x - xr << ", " << y - yr
              << ")\n";
  }

  // This should have worked if we get this far.
  ifail = 0;
  return ifail;
}

int ComponentFieldMap::Coordinates5(const double x, const double y,
                                    const double z, double& t1, double& t2,
                                    double& t3, double& t4, double jac[4][4],
                                    double& det, const Element& element) const {
  // Debugging
  if (m_debug) {
    std::cout << m_className << "::Coordinates5:\n";
    std::cout << "   Point (" << x << ", " << y << ", " << z << ")\n";
  }

  // Failure flag
  int ifail = 1;

  // Provisional values
  t1 = t2 = t3 = t4 = 0;

  // Degenerate elements should have been treated as triangles.
  if (element.degenerate) {
    std::cerr << m_className << "::Coordinates5: Degenerate element.\n";
    return ifail;
  }

  // Make a first order approximation.
  if (Coordinates4(x, y, z, t1, t2, t3, t4, det, element) > 0) {
    if (m_debug) {
      std::cout << m_className << "::Coordinates5:\n";
      std::cout << "    Failure to obtain linear estimate of isoparametric "
                   "coordinates\n.";
    }
    return ifail;
  }

  // Set tolerance parameter.
  constexpr double f = 0.5;
  // Check whether the point is far outside.
  if (t1 < -(1 + f) || t1 > (1 + f) || t2 < -(1 + f) || t2 > (1 + f)) {
    if (m_debug) {
      std::cout << m_className << "::Coordinates5:\n";
      std::cout << "    Point far outside, (t1,t2) = (" << t1 << ", " << t2
                << ").\n";
    }
    return ifail;
  }

  // Start iteration
  double td1 = t1, td2 = t2;
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];
  const Node& n6 = m_nodes[element.emap[6]];
  const Node& n7 = m_nodes[element.emap[7]];
  bool converged = false;
  for (int iter = 0; iter < 10; iter++) {
    if (m_debug) {
      std::cout << m_className << "::Coordinates5:\n";
      std::cout << "    Iteration " << iter << ":     (t1, t2) = (" << td1
                << ", " << td2 << ").\n";
    }
    // Re-compute the (x,y,z) position for this coordinate.
    const double r0 = (-(1 - td1) * (1 - td2) * (1 + td1 + td2)) * 0.25;
    const double r1 = (-(1 + td1) * (1 - td2) * (1 - td1 + td2)) * 0.25;
    const double r2 = (-(1 + td1) * (1 + td2) * (1 - td1 - td2)) * 0.25;
    const double r3 = (-(1 - td1) * (1 + td2) * (1 + td1 - td2)) * 0.25;
    const double r4 = (1 - td1) * (1 + td1) * (1 - td2) * 0.5;
    const double r5 = (1 + td1) * (1 + td2) * (1 - td2) * 0.5;
    const double r6 = (1 - td1) * (1 + td1) * (1 + td2) * 0.5;
    const double r7 = (1 - td1) * (1 + td2) * (1 - td2) * 0.5;
    double xr = n0.x * r0 + n1.x * r1 + n2.x * r2 + n3.x * r3 + n4.x * r4 +
                n5.x * r5 + n6.x * r6 + n7.x * r7;
    double yr = n0.y * r0 + n1.y * r1 + n2.y * r2 + n3.y * r3 + n4.y * r4 +
                n5.y * r5 + n6.y * r6 + n7.y * r7;
    // Compute the Jacobian.
    Jacobian5(element, td1, td2, det, jac);
    // Compute the difference vector.
    double diff[2] = {x - xr, y - yr};
    // Update the estimate.
    double corr[2] = {0., 0.};
    const double invdet = 1. / det;
    for (size_t l = 0; l < 2; ++l) {
      for (size_t k = 0; k < 2; ++k) {
        corr[l] += jac[l][k] * diff[k];
      }
      corr[l] *= invdet;
    }
    // Debugging
    if (m_debug) {
      std::cout << m_className << "::Coordinates5:\n";
      std::cout << "    Difference vector: (x, y)   = (" << diff[0] << ", "
                << diff[1] << ").\n";
      std::cout << "    Correction vector: (t1, t2) = (" << corr[0] << ", "
                << corr[1] << ").\n";
    }
    // Update the vector.
    td1 += corr[0];
    td2 += corr[1];
    // Check for convergence.
    constexpr double tol = 1.e-5;
    if (fabs(corr[0]) < tol && fabs(corr[1]) < tol) {
      if (m_debug) {
        std::cout << m_className << "::Coordinates5: Convergence reached.\n";
      }
      converged = true;
      break;
    }
  }
  // No convergence reached.
  if (!converged) {
    double xmin = std::min({n0.x, n1.x, n2.x, n3.x, n4.x, n5.x, n6.x, n7.x});
    double xmax = std::max({n0.x, n1.x, n2.x, n3.x, n4.x, n5.x, n6.x, n7.x});
    double ymin = std::min({n0.y, n1.y, n2.y, n3.y, n4.y, n5.y, n6.y, n7.y});
    double ymax = std::max({n0.y, n1.y, n2.y, n3.y, n4.y, n5.y, n6.y, n7.y});
    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
      if (m_printConvergenceWarnings) {
        std::cout << m_className << "::Coordinates5:\n"
                  << "    No convergence achieved "
                  << "when refining internal isoparametric coordinates\n"
                  << "    at position (" << x << ", " << y << ").\n";
      }
      t1 = t2 = 0;
      return ifail;
    }
  }

  // Convergence reached.
  t1 = td1;
  t2 = td2;
  t3 = 0;
  t4 = 0;
  if (m_debug) {
    std::cout << m_className << "::Coordinates5:\n";
    std::cout << "    Convergence reached at (t1, t2) = (" << t1 << ", " << t2
              << ").\n";
    // For debugging purposes, show position.
    const double r0 = (-(1 - td1) * (1 - td2) * (1 + td1 + td2)) * 0.25;
    const double r1 = (-(1 + td1) * (1 - td2) * (1 - td1 + td2)) * 0.25;
    const double r2 = (-(1 + td1) * (1 + td2) * (1 - td1 - td2)) * 0.25;
    const double r3 = (-(1 - td1) * (1 + td2) * (1 + td1 - td2)) * 0.25;
    const double r4 = (1 - td1) * (1 + td1) * (1 - td2) * 0.5;
    const double r5 = (1 + td1) * (1 + td2) * (1 - td2) * 0.5;
    const double r6 = (1 - td1) * (1 + td1) * (1 + td2) * 0.5;
    const double r7 = (1 - td1) * (1 + td2) * (1 - td2) * 0.5;
    double xr = n0.x * r0 + n1.x * r1 + n2.x * r2 + n3.x * r3 + n4.x * r4 +
                n5.x * r5 + n6.x * r6 + n7.x * r7;
    double yr = n0.y * r0 + n1.y * r1 + n2.y * r2 + n3.y * r3 + n4.y * r4 +
                n5.y * r5 + n6.y * r6 + n7.y * r7;
    std::cout << "    Position requested:     (" << x << ", " << y << ")\n";
    std::cout << "    Reconstructed:          (" << xr << ", " << yr << ")\n";
    std::cout << "    Difference:             (" << x - xr << ", " << y - yr
              << ")\n";
  }

  // Success
  ifail = 0;
  return ifail;
}

void ComponentFieldMap::Coordinates12(const double x, const double y,
                                      const double z, double& t1, double& t2,
                                      double& t3, double& t4,
                                      const Element& element) const {
  if (m_debug) {
    std::cout << m_className << "::Coordinates12:\n"
              << "   Point (" << x << ", " << y << ", " << z << ").\n";
  }

  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  // Compute tetrahedral coordinates.
  const double f1x =
      (n2.y - n1.y) * (n3.z - n1.z) - (n3.y - n1.y) * (n2.z - n1.z);
  const double f1y =
      (n2.z - n1.z) * (n3.x - n1.x) - (n3.z - n1.z) * (n2.x - n1.x);
  const double f1z =
      (n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y);
  t1 = (x - n1.x) * f1x + (y - n1.y) * f1y + (z - n1.z) * f1z;
  t1 = t1 / ((n0.x - n1.x) * f1x + (n0.y - n1.y) * f1y + (n0.z - n1.z) * f1z);
  const double f2x =
      (n0.y - n2.y) * (n3.z - n2.z) - (n3.y - n2.y) * (n0.z - n2.z);
  const double f2y =
      (n0.z - n2.z) * (n3.x - n2.x) - (n3.z - n2.z) * (n0.x - n2.x);
  const double f2z =
      (n0.x - n2.x) * (n3.y - n2.y) - (n3.x - n2.x) * (n0.y - n2.y);
  t2 = (x - n2.x) * f2x + (y - n2.y) * f2y + (z - n2.z) * f2z;
  t2 = t2 / ((n1.x - n2.x) * f2x + (n1.y - n2.y) * f2y + (n1.z - n2.z) * f2z);
  const double f3x =
      (n0.y - n3.y) * (n1.z - n3.z) - (n1.y - n3.y) * (n0.z - n3.z);
  const double f3y =
      (n0.z - n3.z) * (n1.x - n3.x) - (n1.z - n3.z) * (n0.x - n3.x);
  const double f3z =
      (n0.x - n3.x) * (n1.y - n3.y) - (n1.x - n3.x) * (n0.y - n3.y);
  t3 = (x - n3.x) * f3x + (y - n3.y) * f3y + (z - n3.z) * f3z;
  t3 = t3 / ((n2.x - n3.x) * f3x + (n2.y - n3.y) * f3y + (n2.z - n3.z) * f3z);
  const double f4x =
      (n2.y - n0.y) * (n1.z - n0.z) - (n1.y - n0.y) * (n2.z - n0.z);
  const double f4y =
      (n2.z - n0.z) * (n1.x - n0.x) - (n1.z - n0.z) * (n2.x - n0.x);
  const double f4z =
      (n2.x - n0.x) * (n1.y - n0.y) - (n1.x - n0.x) * (n2.y - n0.y);
  t4 = (x - n0.x) * f4x + (y - n0.y) * f4y + (z - n0.z) * f4z;
  t4 = t4 / ((n3.x - n0.x) * f4x + (n3.y - n0.y) * f4y + (n3.z - n0.z) * f4z);

  // Result
  if (m_debug) {
    std::cout << m_className << "::Coordinates12:\n";
    std::cout << "    Tetrahedral coordinates (t, u, v, w) = (" << t1 << ", "
              << t2 << ", " << t3 << ", " << t4
              << ") sum = " << t1 + t2 + t3 + t4 << ".\n";
    // Re-compute the (x,y,z) position for this coordinate.
    const double xr = n0.x * t1 + n1.x * t2 + n2.x * t3 + n3.x * t4;
    const double yr = n0.y * t1 + n1.y * t2 + n2.y * t3 + n3.y * t4;
    const double zr = n0.z * t1 + n1.z * t2 + n2.z * t3 + n3.z * t4;
    const double sr = t1 + t2 + t3 + t4;
    std::cout << "    Position requested:     (" << x << ", " << y << ", " << z
              << ")\n";
    std::cout << "    Reconstructed:          (" << xr << ", " << yr << ", "
              << zr << ")\n";
    std::cout << "    Difference:             (" << x - xr << ", " << y - yr
              << ", " << z - zr << ")\n";
    std::cout << "    Checksum - 1:           " << sr - 1 << "\n";
  }
}

int ComponentFieldMap::Coordinates13(const double x, const double y,
                                     const double z, double& t1, double& t2,
                                     double& t3, double& t4, double jac[4][4],
                                     double& det,
                                     const Element& element) const {
  if (m_debug) {
    std::cout << m_className << "::Coordinates13:\n";
    std::cout << "   Point (" << x << ", " << y << ", " << z << ")\n";
  }

  // Provisional values
  t1 = t2 = t3 = t4 = 0.;

  // Make a first order approximation.
  Coordinates12(x, y, z, t1, t2, t3, t4, element);

  // Set tolerance parameter.
  constexpr double f = 0.5;
  if (t1 < -f || t2 < -f || t3 < -f || t4 < -f || t1 > 1 + f || t2 > 1 + f ||
      t3 > 1 + f || t4 > 1 + f) {
    if (m_debug) {
      std::cout << m_className << "::Coordinates13:\n";
      std::cout << "    Linear isoparametric coordinates more than\n";
      std::cout << "    f (" << f << ") out of range.\n";
    }
    return 0;
  }

  // Start iteration.
  double td1 = t1, td2 = t2, td3 = t3, td4 = t4;
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

  // Loop
  bool converged = false;
  double diff[4], corr[4];
  for (int iter = 0; iter < 10; iter++) {
    if (m_debug) {
      std::cout << m_className << "::Coordinates13:\n";
      std::printf("    Iteration %4u: t = (%15.8f, %15.8f %15.8f %15.8f)\n",
                  iter, td1, td2, td3, td4);
    }
    const double f0 = td1 * (2 * td1 - 1);
    const double f1 = td2 * (2 * td2 - 1);
    const double f2 = td3 * (2 * td3 - 1);
    const double f3 = td4 * (2 * td4 - 1);
    const double f4 = 4 * td1 * td2;
    const double f5 = 4 * td1 * td3;
    const double f6 = 4 * td1 * td4;
    const double f7 = 4 * td2 * td3;
    const double f8 = 4 * td2 * td4;
    const double f9 = 4 * td3 * td4;
    // Re-compute the (x,y,z) position for this coordinate.
    const double xr = n0.x * f0 + n1.x * f1 + n2.x * f2 + n3.x * f3 +
                      n4.x * f4 + n5.x * f5 + n6.x * f6 + n7.x * f7 +
                      n8.x * f8 + n9.x * f9;
    const double yr = n0.y * f0 + n1.y * f1 + n2.y * f2 + n3.y * f3 +
                      n4.y * f4 + n5.y * f5 + n6.y * f6 + n7.y * f7 +
                      n8.y * f8 + n9.y * f9;
    const double zr = n0.z * f0 + n1.z * f1 + n2.z * f2 + n3.z * f3 +
                      n4.z * f4 + n5.z * f5 + n6.z * f6 + n7.z * f7 +
                      n8.z * f8 + n9.z * f9;
    const double sr = td1 + td2 + td3 + td4;

    // Compute the Jacobian.
    Jacobian13(element, td1, td2, td3, td4, det, jac);
    // Compute the difference vector.
    diff[0] = 1 - sr;
    diff[1] = x - xr;
    diff[2] = y - yr;
    diff[3] = z - zr;
    // Update the estimate.
    const double invdet = 1. / det;
    for (int l = 0; l < 4; ++l) {
      corr[l] = 0;
      for (int k = 0; k < 4; ++k) {
        corr[l] += jac[l][k] * diff[k];
      }
      corr[l] *= invdet;
    }

    // Debugging
    if (m_debug) {
      std::cout << m_className << "::Coordinates13:\n";
      std::cout << "    Difference vector:  (1, x, y, z)  = (" << diff[0]
                << ", " << diff[1] << ", " << diff[2] << ", " << diff[3]
                << ").\n";
      std::cout << "    Correction vector:  (t1,t2,t3,t4) = (" << corr[0]
                << ", " << corr[1] << ", " << corr[2] << ", " << corr[3]
                << ").\n";
    }

    // Update the vector.
    td1 += corr[0];
    td2 += corr[1];
    td3 += corr[2];
    td4 += corr[3];

    // Check for convergence.
    constexpr double tol = 1.e-5;
    if (fabs(corr[0]) < tol && fabs(corr[1]) < tol && fabs(corr[2]) < tol &&
        fabs(corr[3]) < tol) {
      if (m_debug) {
        std::cout << m_className << "::Coordinates13: Convergence reached.\n";
      }
      converged = true;
      break;
    }
  }

  // No convergence reached.
  if (!converged) {
    const double xmin = std::min({n0.x, n1.x, n2.x, n3.x});
    const double xmax = std::max({n0.x, n1.x, n2.x, n3.x});
    const double ymin = std::min({n0.y, n1.y, n2.y, n3.y});
    const double ymax = std::max({n0.y, n1.y, n2.y, n3.y});
    const double zmin = std::min({n0.z, n1.z, n2.z, n3.z});
    const double zmax = std::max({n0.z, n1.z, n2.z, n3.z});
    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax && z >= zmin &&
        z <= zmax) {
      if (m_printConvergenceWarnings) {
        std::cout << m_className << "::Coordinates13:\n"
                  << "    No convergence achieved "
                  << "when refining internal isoparametric coordinates\n"
                  << "    at position (" << x << ", " << y << ", " << z
                  << ").\n";
      }
      t1 = t2 = t3 = t4 = -1;
      return 1;
    }
  }

  // Convergence reached.
  t1 = td1;
  t2 = td2;
  t3 = td3;
  t4 = td4;
  if (m_debug) {
    std::cout << m_className << "::Coordinates13:\n";
    std::cout << "    Convergence reached at (t1, t2, t3, t4) = (" << t1 << ", "
              << t2 << ", " << t3 << ", " << t4 << ").\n";
    // Re-compute the (x,y,z) position for this coordinate.
    const double f0 = td1 * (2 * td1 - 1);
    const double f1 = td2 * (2 * td2 - 1);
    const double f2 = td3 * (2 * td3 - 1);
    const double f3 = td4 * (2 * td4 - 1);
    const double f4 = 4 * td1 * td2;
    const double f5 = 4 * td1 * td3;
    const double f6 = 4 * td1 * td4;
    const double f7 = 4 * td2 * td3;
    const double f8 = 4 * td2 * td4;
    const double f9 = 4 * td3 * td4;
    double xr = n0.x * f0 + n1.x * f1 + n2.x * f2 + n3.x * f3 + n4.x * f4 +
                n5.x * f5 + n6.x * f6 + n7.x * f7 + n8.x * f8 + n9.x * f9;
    double yr = n0.y * f0 + n1.y * f1 + n2.y * f2 + n3.y * f3 + n4.y * f4 +
                n5.y * f5 + n6.y * f6 + n7.y * f7 + n8.y * f8 + n9.y * f9;
    double zr = n0.z * f0 + n1.z * f1 + n2.z * f2 + n3.z * f3 + n4.z * f4 +
                n5.z * f5 + n6.z * f6 + n7.z * f7 + n8.z * f8 + n9.z * f9;
    double sr = td1 + td2 + td3 + td4;
    std::cout << "    Position requested:     (" << x << ", " << y << ", " << z
              << ")\n";
    std::cout << "    Reconstructed:          (" << xr << ", " << yr << ", "
              << zr << ")\n";
    std::cout << "    Difference:             (" << x - xr << ", " << y - yr
              << ", " << z - zr << ")\n";
    std::cout << "    Checksum - 1:           " << sr - 1 << "\n";
  }

  // Success
  return 0;
}

int ComponentFieldMap::CoordinatesCube(const double x, const double y,
                                       const double z, double& t1, double& t2,
                                       double& t3, TMatrixD*& jac,
                                       std::vector<TMatrixD*>& dN,
                                       const Element& element) const {
  /*
  global coordinates   7__ _ _ 6     t3    t2
                      /       /|     ^   /|
    ^ z              /       / |     |   /
    |               4_______5  |     |  /
    |              |        |  |     | /
    |              |  3     |  2     |/     t1
     ------->      |        | /       ------->
    /      y       |        |/       local coordinates
   /               0--------1
  /
 v x
 */

  // Failure flag
  int ifail = 1;

  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n7 = m_nodes[element.emap[7]];

  // Compute hexahedral coordinates (t1->[-1,1],t2->[-1,1],t3->[-1,1]) and
  // t1 (zeta) is in y-direction
  // t2 (eta)  is in opposite x-direction
  // t3 (mu)   is in z-direction
  // Nodes are set in that way, that node [0] has always lowest x,y,z!
  t2 = (2. * (x - n3.x) / (n0.x - n3.x) - 1) * -1.;
  t1 = 2. * (y - n3.y) / (n2.y - n3.y) - 1;
  t3 = 2. * (z - n3.z) / (n7.z - n3.z) - 1;
  // Re-compute the (x,y,z) position for this coordinate.
  if (m_debug) {
    double n[8];
    n[0] = 1. / 8 * (1 - t1) * (1 - t2) * (1 - t3);
    n[1] = 1. / 8 * (1 + t1) * (1 - t2) * (1 - t3);
    n[2] = 1. / 8 * (1 + t1) * (1 + t2) * (1 - t3);
    n[3] = 1. / 8 * (1 - t1) * (1 + t2) * (1 - t3);
    n[4] = 1. / 8 * (1 - t1) * (1 - t2) * (1 + t3);
    n[5] = 1. / 8 * (1 + t1) * (1 - t2) * (1 + t3);
    n[6] = 1. / 8 * (1 + t1) * (1 + t2) * (1 + t3);
    n[7] = 1. / 8 * (1 - t1) * (1 + t2) * (1 + t3);

    double xr = 0;
    double yr = 0;
    double zr = 0;

    for (int i = 0; i < 8; i++) {
      const Node& node = m_nodes[element.emap[i]];
      xr += node.x * n[i];
      yr += node.y * n[i];
      zr += node.z * n[i];
    }
    double sr = n[0] + n[1] + n[2] + n[3] + n[4] + n[5] + n[6] + n[7];
    std::cout << m_className << "::CoordinatesCube:\n";
    std::cout << "    Position requested:     (" << x << "," << y << "," << z
              << ")\n";
    std::cout << "    Position reconstructed: (" << xr << "," << yr << "," << zr
              << ")\n";
    std::cout << "    Difference:             (" << (x - xr) << "," << (y - yr)
              << "," << (z - zr) << ")\n";
    std::cout << "    Hexahedral coordinates (t, u, v) = (" << t1 << "," << t2
              << "," << t3 << ")\n";
    std::cout << "    Checksum - 1:           " << (sr - 1) << "\n";
  }
  if (jac != 0) JacobianCube(element, t1, t2, t3, jac, dN);
  // This should always work.
  ifail = 0;
  return ifail;
}

void ComponentFieldMap::Reset() {

  m_ready = false;

  m_elements.clear();
  m_nodes.clear();
  m_materials.clear();
  m_wfields.clear();
  m_wfieldsOk.clear();
  m_hasBoundingBox = false;
  m_warning = false;
  m_nWarnings = 0;

  m_octree.reset(nullptr);
  m_cacheElemBoundingBoxes = false;
  m_lastElement = -1;
}

void ComponentFieldMap::Prepare() {

  // Establish the ranges.
  SetRange();
  UpdatePeriodicity();
  std::cout << m_className << "::Prepare:\n"
            << "    Caching the bounding boxes of all elements...";
  CalculateElementBoundingBoxes();
  std::cout << " done.\n";
  // Initialize the tetrahedral tree.
  InitializeTetrahedralTree();
}

void ComponentFieldMap::UpdatePeriodicityCommon() {
  // Check the required data is available.
  if (!m_ready) {
    std::cerr << m_className << "::UpdatePeriodicityCommon:\n";
    std::cerr << "    No valid field map available.\n";
    return;
  }

  for (unsigned int i = 0; i < 3; ++i) {
    // No regular and mirror periodicity at the same time.
    if (m_periodic[i] && m_mirrorPeriodic[i]) {
      std::cerr << m_className << "::UpdatePeriodicityCommon:\n"
                << "    Both simple and mirror periodicity requested. Reset.\n";
      m_periodic[i] = false;
      m_mirrorPeriodic[i] = false;
      m_warning = true;
    }
    // In case of axial periodicity,
    // the range must be an integral part of two pi.
    if (m_axiallyPeriodic[i]) {
      if (m_mapamin[i] >= m_mapamax[i]) {
        m_mapna[i] = 0;
      } else {
        m_mapna[i] = TwoPi / (m_mapamax[i] - m_mapamin[i]);
      }
      if (fabs(m_mapna[i] - int(0.5 + m_mapna[i])) > 0.001 ||
          m_mapna[i] < 1.5) {
        std::cerr
            << m_className << "::UpdatePeriodicityCommon:\n"
            << "    Axial symmetry has been requested but the map\n"
            << "    does not cover an integral fraction of 2 pi. Reset.\n";
        m_axiallyPeriodic[i] = false;
        m_warning = true;
      }
    }
  }

  // Not more than 1 rotational symmetry
  if ((m_rotationSymmetric[0] && m_rotationSymmetric[1]) ||
      (m_rotationSymmetric[0] && m_rotationSymmetric[2]) ||
      (m_rotationSymmetric[1] && m_rotationSymmetric[2])) {
    std::cerr << m_className << "::UpdatePeriodicityCommon:\n";
    std::cerr << "    Only 1 rotational symmetry allowed; reset.\n";
    m_rotationSymmetric.fill(false);
    m_warning = true;
  }

  // No rotational symmetry as well as axial periodicity
  if ((m_rotationSymmetric[0] || m_rotationSymmetric[1] ||
       m_rotationSymmetric[2]) &&
      (m_axiallyPeriodic[0] || m_axiallyPeriodic[1] || m_axiallyPeriodic[2])) {
    std::cerr << m_className << "::UpdatePeriodicityCommon:\n";
    std::cerr << "    Not allowed to combine rotational symmetry\n";
    std::cerr << "    and axial periodicity; reset.\n";
    m_axiallyPeriodic.fill(false);
    m_rotationSymmetric.fill(false);
    m_warning = true;
  }

  // In case of rotational symmetry, the x-range should not straddle 0.
  if (m_rotationSymmetric[0] || m_rotationSymmetric[1] ||
      m_rotationSymmetric[2]) {
    if (m_mapmin[0] * m_mapmax[0] < 0) {
      std::cerr << m_className << "::UpdatePeriodicityCommon:\n";
      std::cerr << "    Rotational symmetry requested, \n";
      std::cerr << "    but x-range straddles 0; reset.\n";
      m_rotationSymmetric.fill(false);
      m_warning = true;
    }
  }

  // Recompute the cell ranges.
  for (unsigned int i = 0; i < 3; ++i) {
    m_minBoundingBox[i] = m_mapmin[i];
    m_maxBoundingBox[i] = m_mapmax[i];
    m_cells[i] = fabs(m_mapmax[i] - m_mapmin[i]);
  }
  if (m_rotationSymmetric[0]) {
    m_minBoundingBox[0] = m_mapmin[1];
    m_maxBoundingBox[0] = m_mapmax[1];
    m_minBoundingBox[1] = -std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_maxBoundingBox[1] = +std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_minBoundingBox[2] = -std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_maxBoundingBox[2] = +std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
  } else if (m_rotationSymmetric[1]) {
    m_minBoundingBox[0] = -std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_maxBoundingBox[0] = +std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_minBoundingBox[1] = m_mapmin[1];
    m_maxBoundingBox[1] = m_mapmax[1];
    m_minBoundingBox[2] = -std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_maxBoundingBox[2] = +std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
  } else if (m_rotationSymmetric[2]) {
    m_minBoundingBox[0] = -std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_maxBoundingBox[0] = +std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_minBoundingBox[1] = -std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_maxBoundingBox[1] = +std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0]));
    m_minBoundingBox[2] = m_mapmin[1];
    m_maxBoundingBox[2] = m_mapmax[1];
  }

  if (m_axiallyPeriodic[0]) {
    m_minBoundingBox[1] =
        -std::max(std::max(fabs(m_mapmin[1]), fabs(m_mapmax[1])),
                  std::max(fabs(m_mapmin[2]), fabs(m_mapmax[2])));
    m_maxBoundingBox[1] =
        +std::max(std::max(fabs(m_mapmin[1]), fabs(m_mapmax[1])),
                  std::max(fabs(m_mapmin[2]), fabs(m_mapmax[2])));
    m_minBoundingBox[2] =
        -std::max(std::max(fabs(m_mapmin[1]), fabs(m_mapmax[1])),
                  std::max(fabs(m_mapmin[2]), fabs(m_mapmax[2])));
    m_maxBoundingBox[2] =
        +std::max(std::max(fabs(m_mapmin[1]), fabs(m_mapmax[1])),
                  std::max(fabs(m_mapmin[2]), fabs(m_mapmax[2])));
  } else if (m_axiallyPeriodic[1]) {
    m_minBoundingBox[0] =
        -std::max(std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0])),
                  std::max(fabs(m_mapmin[2]), fabs(m_mapmax[2])));
    m_maxBoundingBox[0] =
        +std::max(std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0])),
                  std::max(fabs(m_mapmin[2]), fabs(m_mapmax[2])));
    m_minBoundingBox[2] =
        -std::max(std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0])),
                  std::max(fabs(m_mapmin[2]), fabs(m_mapmax[2])));
    m_maxBoundingBox[2] =
        +std::max(std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0])),
                  std::max(fabs(m_mapmin[2]), fabs(m_mapmax[2])));
  } else if (m_axiallyPeriodic[2]) {
    m_minBoundingBox[0] =
        -std::max(std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0])),
                  std::max(fabs(m_mapmin[1]), fabs(m_mapmax[1])));
    m_maxBoundingBox[0] =
        +std::max(std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0])),
                  std::max(fabs(m_mapmin[1]), fabs(m_mapmax[1])));
    m_minBoundingBox[1] =
        -std::max(std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0])),
                  std::max(fabs(m_mapmin[1]), fabs(m_mapmax[1])));
    m_maxBoundingBox[1] =
        +std::max(std::max(fabs(m_mapmin[0]), fabs(m_mapmax[0])),
                  std::max(fabs(m_mapmin[1]), fabs(m_mapmax[1])));
  }

  for (unsigned int i = 0; i < 3; ++i) {
    if (m_periodic[i] || m_mirrorPeriodic[i]) {
      m_minBoundingBox[i] = -INFINITY;
      m_maxBoundingBox[i] = +INFINITY;
    }
  }

  // Display the range if requested.
  if (m_debug) PrintRange();
}

void ComponentFieldMap::UpdatePeriodicity2d() {
  // Check the required data is available.
  if (!m_ready) {
    std::cerr << m_className << "::UpdatePeriodicity2d:\n";
    std::cerr << "    No valid field map available.\n";
    return;
  }

  // No z-periodicity in 2d
  if (m_periodic[2] || m_mirrorPeriodic[2]) {
    std::cerr << m_className << "::UpdatePeriodicity2d:\n";
    std::cerr << "    Simple or mirror periodicity along z\n";
    std::cerr << "    requested for a 2d map; reset.\n";
    m_periodic[2] = false;
    m_mirrorPeriodic[2] = false;
    m_warning = true;
  }

  // Only z-axial periodicity in 2d maps
  if (m_axiallyPeriodic[0] || m_axiallyPeriodic[1]) {
    std::cerr << m_className << "::UpdatePeriodicity2d:\n";
    std::cerr << "    Axial symmetry has been requested \n";
    std::cerr << "    around x or y for a 2D map; reset.\n";
    m_axiallyPeriodic[0] = false;
    m_axiallyPeriodic[1] = false;
    m_warning = true;
  }
}

void ComponentFieldMap::SetRange() {
  // Initial values
  m_mapmin.fill(0.);
  m_mapmax.fill(0.);
  m_mapamin.fill(0.);
  m_mapamax.fill(0.);
  m_mapvmin = m_mapvmax = 0.;
  m_setang.fill(false);

  // Make sure the required data is available.
  if (!m_ready || m_nodes.empty()) {
    std::cerr << m_className << "::SetRange: Field map not yet set.\n";
    return;
  }

  // Loop over the nodes.
  m_mapmin[0] = m_mapmax[0] = m_nodes[0].x;
  m_mapmin[1] = m_mapmax[1] = m_nodes[0].y;
  m_mapmin[2] = m_mapmax[2] = m_nodes[0].z;
  m_mapvmin = m_mapvmax = m_nodes[0].v;

  for (const auto& node : m_nodes) {
    const std::array<double, 3> pos = {{node.x, node.y, node.z}};
    for (unsigned int i = 0; i < 3; ++i) {
      m_mapmin[i] = std::min(m_mapmin[i], pos[i]);
      m_mapmax[i] = std::max(m_mapmax[i], pos[i]);
    }
    m_mapvmin = std::min(m_mapvmin, node.v);
    m_mapvmax = std::max(m_mapvmax, node.v);

    if (node.y != 0 || node.z != 0) {
      const double ang = atan2(node.z, node.y);
      if (m_setang[0]) {
        m_mapamin[0] = std::min(m_mapamin[0], ang);
        m_mapamax[0] = std::max(m_mapamax[0], ang);
      } else {
        m_mapamin[0] = m_mapamax[0] = ang;
        m_setang[0] = true;
      }
    }

    if (node.z != 0 || node.x != 0) {
      const double ang = atan2(node.x, node.z);
      if (m_setang[1]) {
        m_mapamin[1] = std::min(m_mapamin[1], ang);
        m_mapamax[1] = std::max(m_mapamax[1], ang);
      } else {
        m_mapamin[1] = m_mapamax[1] = ang;
        m_setang[1] = true;
      }
    }

    if (node.x != 0 || node.y != 0) {
      const double ang = atan2(node.y, node.x);
      if (m_setang[2]) {
        m_mapamin[2] = std::min(m_mapamin[2], ang);
        m_mapamax[2] = std::max(m_mapamax[2], ang);
      } else {
        m_mapamin[2] = m_mapamax[2] = ang;
        m_setang[2] = true;
      }
    }
  }

  // Fix the angular ranges.
  for (unsigned int i = 0; i < 3; ++i) {
    if (m_mapamax[i] - m_mapamin[i] > Pi) {
      const double aux = m_mapamin[i];
      m_mapamin[i] = m_mapamax[i];
      m_mapamax[i] = aux + TwoPi;
    }
  }

  // Set provisional cell dimensions.
  m_minBoundingBox[0] = m_mapmin[0];
  m_maxBoundingBox[0] = m_mapmax[0];
  m_minBoundingBox[1] = m_mapmin[1];
  m_maxBoundingBox[1] = m_mapmax[1];
  if (m_is3d) {
    m_minBoundingBox[2] = m_mapmin[2];
    m_maxBoundingBox[2] = m_mapmax[2];
  } else {
    m_mapmin[2] = m_minBoundingBox[2];
    m_mapmax[2] = m_maxBoundingBox[2];
  }
  m_hasBoundingBox = true;

  // Display the range if requested.
  if (m_debug) PrintRange();
}

void ComponentFieldMap::PrintRange() {
  std::cout << m_className << "::PrintRange:\n";
  std::cout << "        Dimensions of the elementary block\n";
  printf("            %15g < x < %-15g cm,\n", m_mapmin[0], m_mapmax[0]);
  printf("            %15g < y < %-15g cm,\n", m_mapmin[1], m_mapmax[1]);
  printf("            %15g < z < %-15g cm,\n", m_mapmin[2], m_mapmax[2]);
  printf("            %15g < V < %-15g V.\n", m_mapvmin, m_mapvmax);

  std::cout << "        Periodicities\n";
  const std::array<std::string, 3> axes = {{"x", "y", "z"}};
  for (unsigned int i = 0; i < 3; ++i) {
    std::cout << "            " << axes[i] << ":";
    if (m_periodic[i]) {
      std::cout << " simple with length " << m_cells[i] << " cm";
    }
    if (m_mirrorPeriodic[i]) {
      std::cout << " mirror with length " << m_cells[i] << " cm";
    }
    if (m_axiallyPeriodic[i]) {
      std::cout << " axial " << int(0.5 + m_mapna[i]) << "-fold repetition";
    }
    if (m_rotationSymmetric[i]) std::cout << " rotational symmetry";
    if (!(m_periodic[i] || m_mirrorPeriodic[i] || m_axiallyPeriodic[i] ||
          m_rotationSymmetric[i]))
      std::cout << " none";
    std::cout << "\n";
  }
}

bool ComponentFieldMap::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                       double& xmax, double& ymax,
                                       double& zmax) {
  if (!m_ready) return false;

  xmin = m_minBoundingBox[0];
  xmax = m_maxBoundingBox[0];
  ymin = m_minBoundingBox[1];
  ymax = m_maxBoundingBox[1];
  zmin = m_minBoundingBox[2];
  zmax = m_maxBoundingBox[2];
  return true;
}

void ComponentFieldMap::MapCoordinates(double& xpos, double& ypos, double& zpos,
                                       bool& xmirrored, bool& ymirrored,
                                       bool& zmirrored, double& rcoordinate,
                                       double& rotation) const {
  // Initial values
  rotation = 0;

  // If chamber is periodic, reduce to the cell volume.
  xmirrored = false;
  if (m_periodic[0]) {
    const double xrange = m_mapmax[0] - m_mapmin[0];
    xpos = m_mapmin[0] + fmod(xpos - m_mapmin[0], xrange);
    if (xpos < m_mapmin[0]) xpos += xrange;
  } else if (m_mirrorPeriodic[0]) {
    const double xrange = m_mapmax[0] - m_mapmin[0];
    double xnew = m_mapmin[0] + fmod(xpos - m_mapmin[0], xrange);
    if (xnew < m_mapmin[0]) xnew += xrange;
    int nx = int(floor(0.5 + (xnew - xpos) / xrange));
    if (nx != 2 * (nx / 2)) {
      xnew = m_mapmin[0] + m_mapmax[0] - xnew;
      xmirrored = true;
    }
    xpos = xnew;
  }
  if (m_axiallyPeriodic[0] && (zpos != 0 || ypos != 0)) {
    const double auxr = sqrt(zpos * zpos + ypos * ypos);
    double auxphi = atan2(zpos, ypos);
    const double phirange = m_mapamax[0] - m_mapamin[0];
    const double phim = 0.5 * (m_mapamin[0] + m_mapamax[0]);
    rotation = phirange * floor(0.5 + (auxphi - phim) / phirange);
    if (auxphi - rotation < m_mapamin[0]) rotation -= phirange;
    if (auxphi - rotation > m_mapamax[0]) rotation += phirange;
    auxphi = auxphi - rotation;
    ypos = auxr * cos(auxphi);
    zpos = auxr * sin(auxphi);
  }

  ymirrored = false;
  if (m_periodic[1]) {
    const double yrange = m_mapmax[1] - m_mapmin[1];
    ypos = m_mapmin[1] + fmod(ypos - m_mapmin[1], yrange);
    if (ypos < m_mapmin[1]) ypos += yrange;
  } else if (m_mirrorPeriodic[1]) {
    const double yrange = m_mapmax[1] - m_mapmin[1];
    double ynew = m_mapmin[1] + fmod(ypos - m_mapmin[1], yrange);
    if (ynew < m_mapmin[1]) ynew += yrange;
    int ny = int(floor(0.5 + (ynew - ypos) / yrange));
    if (ny != 2 * (ny / 2)) {
      ynew = m_mapmin[1] + m_mapmax[1] - ynew;
      ymirrored = true;
    }
    ypos = ynew;
  }
  if (m_axiallyPeriodic[1] && (xpos != 0 || zpos != 0)) {
    const double auxr = sqrt(xpos * xpos + zpos * zpos);
    double auxphi = atan2(xpos, zpos);
    const double phirange = (m_mapamax[1] - m_mapamin[1]);
    const double phim = 0.5 * (m_mapamin[1] + m_mapamax[1]);
    rotation = phirange * floor(0.5 + (auxphi - phim) / phirange);
    if (auxphi - rotation < m_mapamin[1]) rotation -= phirange;
    if (auxphi - rotation > m_mapamax[1]) rotation += phirange;
    auxphi = auxphi - rotation;
    zpos = auxr * cos(auxphi);
    xpos = auxr * sin(auxphi);
  }

  zmirrored = false;
  if (m_periodic[2]) {
    const double zrange = m_mapmax[2] - m_mapmin[2];
    zpos = m_mapmin[2] + fmod(zpos - m_mapmin[2], zrange);
    if (zpos < m_mapmin[2]) zpos += zrange;
  } else if (m_mirrorPeriodic[2]) {
    const double zrange = m_mapmax[2] - m_mapmin[2];
    double znew = m_mapmin[2] + fmod(zpos - m_mapmin[2], zrange);
    if (znew < m_mapmin[2]) znew += zrange;
    int nz = int(floor(0.5 + (znew - zpos) / zrange));
    if (nz != 2 * (nz / 2)) {
      znew = m_mapmin[2] + m_mapmax[2] - znew;
      zmirrored = true;
    }
    zpos = znew;
  }
  if (m_axiallyPeriodic[2] && (ypos != 0 || xpos != 0)) {
    const double auxr = sqrt(ypos * ypos + xpos * xpos);
    double auxphi = atan2(ypos, xpos);
    const double phirange = m_mapamax[2] - m_mapamin[2];
    const double phim = 0.5 * (m_mapamin[2] + m_mapamax[2]);
    rotation = phirange * floor(0.5 + (auxphi - phim) / phirange);
    if (auxphi - rotation < m_mapamin[2]) rotation -= phirange;
    if (auxphi - rotation > m_mapamax[2]) rotation += phirange;
    auxphi = auxphi - rotation;
    xpos = auxr * cos(auxphi);
    ypos = auxr * sin(auxphi);
  }

  // If we have a rotationally symmetric field map, store coordinates.
  rcoordinate = 0;
  double zcoordinate = 0;
  if (m_rotationSymmetric[0]) {
    rcoordinate = sqrt(ypos * ypos + zpos * zpos);
    zcoordinate = xpos;
  } else if (m_rotationSymmetric[1]) {
    rcoordinate = sqrt(xpos * xpos + zpos * zpos);
    zcoordinate = ypos;
  } else if (m_rotationSymmetric[2]) {
    rcoordinate = sqrt(xpos * xpos + ypos * ypos);
    zcoordinate = zpos;
  }

  if (m_rotationSymmetric[0] || m_rotationSymmetric[1] ||
      m_rotationSymmetric[2]) {
    xpos = rcoordinate;
    ypos = zcoordinate;
    zpos = 0;
  }
}

void ComponentFieldMap::UnmapFields(double& ex, double& ey, double& ez,
                                    double& xpos, double& ypos, double& zpos,
                                    bool& xmirrored, bool& ymirrored,
                                    bool& zmirrored, double& rcoordinate,
                                    double& rotation) const {
  // Apply mirror imaging.
  if (xmirrored) ex = -ex;
  if (ymirrored) ey = -ey;
  if (zmirrored) ez = -ez;

  // Rotate the field.
  double er, theta;
  if (m_axiallyPeriodic[0]) {
    er = sqrt(ey * ey + ez * ez);
    theta = atan2(ez, ey);
    theta += rotation;
    ey = er * cos(theta);
    ez = er * sin(theta);
  }
  if (m_axiallyPeriodic[1]) {
    er = sqrt(ez * ez + ex * ex);
    theta = atan2(ex, ez);
    theta += rotation;
    ez = er * cos(theta);
    ex = er * sin(theta);
  }
  if (m_axiallyPeriodic[2]) {
    er = sqrt(ex * ex + ey * ey);
    theta = atan2(ey, ex);
    theta += rotation;
    ex = er * cos(theta);
    ey = er * sin(theta);
  }

  // Take care of symmetry.
  double eaxis;
  er = ex;
  eaxis = ey;

  // Rotational symmetry
  if (m_rotationSymmetric[0]) {
    if (rcoordinate <= 0) {
      ex = eaxis;
      ey = 0;
      ez = 0;
    } else {
      ex = eaxis;
      ey = er * ypos / rcoordinate;
      ez = er * zpos / rcoordinate;
    }
  }
  if (m_rotationSymmetric[1]) {
    if (rcoordinate <= 0) {
      ex = 0;
      ey = eaxis;
      ez = 0;
    } else {
      ex = er * xpos / rcoordinate;
      ey = eaxis;
      ez = er * zpos / rcoordinate;
    }
  }
  if (m_rotationSymmetric[2]) {
    if (rcoordinate <= 0) {
      ex = 0;
      ey = 0;
      ez = eaxis;
    } else {
      ex = er * xpos / rcoordinate;
      ey = er * ypos / rcoordinate;
      ez = eaxis;
    }
  }
}

double ComponentFieldMap::ScalingFactor(std::string unit) {
  std::transform(unit.begin(), unit.end(), unit.begin(), toupper);
  if (unit == "MUM" || unit == "MICRON" || unit == "MICROMETER") {
    return 0.0001;
  } else if (unit == "MM" || unit == "MILLIMETER") {
    return 0.1;
  } else if (unit == "CM" || unit == "CENTIMETER") {
    return 1.0;
  } else if (unit == "M" || unit == "METER") {
    return 100.0;
  }
  return -1.;
}

int ComponentFieldMap::ReadInteger(char* token, int def, bool& error) {
  if (!token) {
    error = true;
    return def;
  }

  return atoi(token);
}

double ComponentFieldMap::ReadDouble(char* token, double def, bool& error) {
  if (!token) {
    error = true;
    return def;
  }
  return atof(token);
}

void ComponentFieldMap::CalculateElementBoundingBoxes() {
  // Do not proceed if not properly initialised.
  if (!m_ready) {
    PrintNotReady("CalculateElementBoundingBoxes");
    return;
  }

  // Calculate the bounding boxes of all elements
  for (auto& element : m_elements) {
    const Node& n0 = m_nodes[element.emap[0]];
    const Node& n1 = m_nodes[element.emap[1]];
    const Node& n2 = m_nodes[element.emap[2]];
    const Node& n3 = m_nodes[element.emap[3]];
    element.bbMin[0] = std::min({n0.x, n1.x, n2.x, n3.x});
    element.bbMax[0] = std::max({n0.x, n1.x, n2.x, n3.x});
    element.bbMin[1] = std::min({n0.y, n1.y, n2.y, n3.y});
    element.bbMax[1] = std::max({n0.y, n1.y, n2.y, n3.y});
    element.bbMin[2] = std::min({n0.z, n1.z, n2.z, n3.z});
    element.bbMax[2] = std::max({n0.z, n1.z, n2.z, n3.z});
    // Add tolerances.
    constexpr float f = 0.2;
    const float tolx = f * (element.bbMax[0] - element.bbMin[0]);
    element.bbMin[0] -= tolx;
    element.bbMax[0] += tolx;
    const float toly = f * (element.bbMax[1] - element.bbMin[1]);
    element.bbMin[1] -= toly;
    element.bbMax[1] += toly;
    const float tolz = f * (element.bbMax[2] - element.bbMin[2]);
    element.bbMin[2] -= tolz;
    element.bbMax[2] += tolz;
  }
}

bool ComponentFieldMap::InitializeTetrahedralTree() {
  // Do not proceed if not properly initialised.
  if (!m_ready) {
    PrintNotReady("InitializeTetrahedralTree");
    return false;
  }

  if (m_debug) {
    std::cout << m_className << "::InitializeTetrahedralTree:\n"
              << "    About to initialize the tetrahedral tree.\n";
  }

  // Cache the bounding boxes if it has not been done yet.
  if (!m_cacheElemBoundingBoxes) CalculateElementBoundingBoxes();

  if (m_nodes.empty()) {
    std::cerr << m_className << "::InitializeTetrahedralTree: Empty mesh.\n";
    return false;
  }

  // Determine the bounding box
  double xmin = m_nodes.front().x;
  double ymin = m_nodes.front().y;
  double zmin = m_nodes.front().z;
  double xmax = xmin;
  double ymax = ymin;
  double zmax = zmin;
  for (const auto& node : m_nodes) {
    xmin = std::min(xmin, node.x);
    xmax = std::max(xmax, node.x);
    ymin = std::min(ymin, node.y);
    ymax = std::max(ymax, node.y);
    zmin = std::min(zmin, node.z);
    zmax = std::max(zmax, node.z);
  }

  if (m_debug) {
    std::cout << "    Bounding box:\n"
              << std::scientific << "\tx: " << xmin << " -> " << xmax << "\n"
              << std::scientific << "\ty: " << ymin << " -> " << ymax << "\n"
              << std::scientific << "\tz: " << zmin << " -> " << zmax << "\n";
  }

  const double hx = 0.5 * (xmax - xmin);
  const double hy = 0.5 * (ymax - ymin);
  const double hz = 0.5 * (zmax - zmin);
  m_octree.reset(new TetrahedralTree(Vec3(xmin + hx, ymin + hy, zmin + hz),
                                     Vec3(hx, hy, hz)));

  if (m_debug) std::cout << "    Tree instantiated.\n";

  // Insert all mesh nodes in the tree
  for (unsigned int i = 0; i < m_nodes.size(); i++) {
    const Node& n = m_nodes[i];
    m_octree->InsertMeshNode(Vec3(n.x, n.y, n.z), i);
  }

  if (m_debug) std::cout << "    Tree nodes initialized successfully.\n";

  // Insert all mesh elements (tetrahedrons) in the tree
  for (unsigned int i = 0; i < m_elements.size(); i++) {
    const Element& e = m_elements[i];
    const double bb[6] = {e.bbMin[0], e.bbMin[1], e.bbMin[2], 
                          e.bbMax[0], e.bbMax[1], e.bbMax[2]};
    m_octree->InsertMeshElement(bb, i);
  }

  std::cout << m_className << "::InitializeTetrahedralTree: Success.\n";
  return true;
}

size_t ComponentFieldMap::GetWeightingFieldIndex(
    const std::string& label) const {
  const size_t nWeightingFields = m_wfields.size();
  for (size_t i = 0; i < nWeightingFields; ++i) {
    if (m_wfields[i] == label) return i;
  }
  return nWeightingFields;
}

size_t ComponentFieldMap::GetOrCreateWeightingFieldIndex(
    const std::string& label) {
  // Check if a weighting field with the same label already exists.
  size_t nWeightingFields = m_wfields.size();
  for (size_t i = 0; i < nWeightingFields; ++i) {
    if (m_wfields[i] == label) return i;
  }
  ++nWeightingFields;
  m_wfields.resize(nWeightingFields);
  m_wfieldsOk.resize(nWeightingFields);
  for (auto& node : m_nodes) {
    node.w.resize(nWeightingFields);
    node.dw.resize(nWeightingFields);
  }
  m_wfields.back() = label;
  return nWeightingFields - 1;
}

void ComponentFieldMap::PrintWarning(const std::string& header) {

  if (!m_warning || m_nWarnings > 10) return;
  std::cerr << m_className << "::" << header << ":\n"
            << "    Warnings have been issued for this field map.\n";
  ++m_nWarnings;
}

void ComponentFieldMap::PrintNotReady(const std::string& header) const {
  std::cerr << m_className << "::" << header << ":\n"
            << "    Field map not yet initialised.\n";
}

void ComponentFieldMap::PrintCouldNotOpen(const std::string& header, 
                                          const std::string& filename) const {
  std::cerr << m_className << "::" << header << ":\n"
            << "    Could not open file " << filename << " for reading.\n"
            << "    The file perhaps does not exist.\n";
}

void ComponentFieldMap::PrintElement(const std::string& header, const double x,
                                     const double y, const double z,
                                     const double t1, const double t2,
                                     const double t3, const double t4,
                                     const Element& element,
                                     const unsigned int n, const int iw) const {
  std::cout << m_className << "::" << header << ":\n"
            << "    Global = (" << x << ", " << y << ", " << z << ")\n"
            << "    Local = (" << t1 << ", " << t2 << ", " << t3 << ", " << t4
            << ")\n";
  if (element.degenerate) std::cout << "    Element is degenerate.\n";
  std::cout << " Node             x            y            z            V\n";
  for (unsigned int ii = 0; ii < n; ++ii) {
    const Node& node = m_nodes[element.emap[ii]];
    const double v = iw < 0 ? node.v : node.w[iw];
    printf("      %-5d %12g %12g %12g %12g\n", element.emap[ii], node.x, node.y,
           node.z, v);
  }
}
}  // namespace Garfield
