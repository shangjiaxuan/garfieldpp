#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

#include "Garfield/ComponentAnsys121.hh"

namespace Garfield {

ComponentAnsys121::ComponentAnsys121() : ComponentFieldMap("Ansys121") {
  m_is3d = false;
  m_elementType = ElementType::Serendipity;
  // Default bounding box
  m_minBoundingBox[2] = -50;
  m_maxBoundingBox[2] = 50;
}

bool ComponentAnsys121::Initialise(const std::string& elist, 
                                   const std::string& nlist,
                                   const std::string& mplist, 
                                   const std::string& prnsol,
                                   const std::string& unit) {
  Reset();
  // Keep track of the success.
  bool ok = true;

  // Buffer for reading
  constexpr int size = 100;
  char line[size];

  // Open the material list.
  std::ifstream fmplist(mplist);
  if (!fmplist) {
    PrintCouldNotOpen("Initialise", mplist);
    return false;
  }

  // Read the material list.
  int il = 0;
  unsigned int icurrmat = 0;
  bool readerror = false;
  while (fmplist.getline(line, size, '\n')) {
    il++;
    // Skip page feed.
    if (strcmp(line, "1") == 0) {
      for (size_t k = 0; k < 5; ++k) fmplist.getline(line, size, '\n');
      il += 5;
      continue;
    }
    // Split the line in tokens.
    char* token = strtok(line, " ");
    // Skip blank lines and headers.
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        strcmp(token, "TEMPERATURE") == 0 || strcmp(token, "PROPERTY=") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13) {
      continue;
    }
    // Read number of materials and initialise the list.
    if (strcmp(token, "LIST") == 0) {
      token = strtok(nullptr, " ");
      token = strtok(nullptr, " ");
      token = strtok(nullptr, " ");
      token = strtok(nullptr, " ");
      const int nMaterials = ReadInteger(token, -1, readerror);
      if (readerror || nMaterials < 0) {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        return false;
      }
      m_materials.resize(nMaterials);
      for (auto& material : m_materials) {
        material.ohm = -1;
        material.eps = -1;
        material.medium = nullptr;
      }
      if (m_debug) {
        std::cout << m_className << "::Initialise: " << nMaterials
                  << " materials.\n";
      }
    } else if (strcmp(token, "MATERIAL") == 0) {
      // Version 12 format: read material number
      token = strtok(nullptr, " ");
      token = strtok(nullptr, " ");
      const int imat = ReadInteger(token, -1, readerror);
      if (readerror || imat < 0) {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        return false;
      }
      icurrmat = imat;
    } else if (strcmp(token, "TEMP") == 0) {
      // Version 12 format: read property tag and value
      token = strtok(nullptr, " ");
      int itype = 0;
      if (strncmp(token, "PERX", 4) == 0) {
        itype = 1;
      } else if (strncmp(token, "RSVX", 4) == 0) {
        itype = 2;
      } else {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Unknown material property flag " << token << "\n"
                  << "    in material properties file " << mplist << " (line "
                  << il << ").\n";
        ok = false;
        break;
      }
      fmplist.getline(line, size, '\n');
      il++;
      token = nullptr;
      token = strtok(line, " ");
      if (icurrmat < 1 || icurrmat > m_materials.size()) {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Found out-of-range current material index "
                  << icurrmat << "\n"
                  << "    in material properties file " << mplist << ".\n";
        ok = false;
        break;
      }
      if (itype == 1) {
        m_materials[icurrmat - 1].eps = ReadDouble(token, -1, readerror);
      } else if (itype == 2) {
        m_materials[icurrmat - 1].ohm = ReadDouble(token, -1, readerror);
      }
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        return false;
      }
    } else if (strcmp(token, "PROPERTY") == 0) {
      // Version 11 format
      token = strtok(nullptr, " ");
      token = strtok(nullptr, " ");
      int itype = 0;
      if (strcmp(token, "PERX") == 0) {
        itype = 1;
      } else if (strcmp(token, "RSVX") == 0) {
        itype = 2;
      } else {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Unknown material property flag " << token << "\n"
                  << "    in material properties file " << mplist << " (line "
                  << il << ").\n";
        ok = false;
      }
      token = strtok(nullptr, " ");
      token = strtok(nullptr, " ");
      int imat = ReadInteger(token, -1, readerror);
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        return false;
      } else if (imat < 1 || imat > (int)m_materials.size()) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Found out-of-range current material index " << imat
                  << "\n";
        std::cerr << "    in material properties file " << mplist << ".\n";
        ok = false;
        break;
      } else {
        fmplist.getline(line, size, '\n');
        il++;
        fmplist.getline(line, size, '\n');
        il++;
        token = strtok(line, " ");
        token = strtok(nullptr, " ");
        if (itype == 1) {
          m_materials[imat - 1].eps = ReadDouble(token, -1, readerror);
        } else if (itype == 2) {
          m_materials[imat - 1].ohm = ReadDouble(token, -1, readerror);
        }
        if (readerror) {
          std::cerr << m_className << "::Initialise:\n"
                    << "    Error reading file " << mplist << " (line " << il
                    << ").\n";
          fmplist.close();
          return false;
        }
      }
    }
  }
  // Close the file
  fmplist.close();
  if (!ok) return false;

  // Find lowest epsilon, check for eps = 0, set default drift medium.
  if (!SetDefaultDriftMedium()) return false;

  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read properties of " << m_materials.size()
            << " materials from file " << mplist << ".\n";
  if (m_debug) PrintMaterials();

  // Open the element list
  std::ifstream felist(elist);
  if (!felist) {
    PrintCouldNotOpen("Initialise", elist);
    return false;
  }
  // Read the element list
  int ndegenerate = 0;
  int nbackground = 0;
  il = 0;
  int highestnode = 0;
  while (felist.getline(line, size, '\n')) {
    il++;
    // Skip page feed in batch page title.
    if (strstr(line, "VERSION") != nullptr) {
      for (size_t k = 0; k < 2; ++k) felist.getline(line, size, '\n');
      il += 2;
      continue;
    }
    // Split the line in tokens.
    char* token = strtok(line, " ");
    // Skip blank lines and headers.
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "LIST") == 0 || strcmp(token, "ELEM") == 0 ||
        strcmp(token, "ANSYS") == 0 || strcmp(token, "***") == 0 ||
        strcmp(token, "VERSION") == 0) {
      continue;
    }
    // Read the element.
    int ielem = ReadInteger(token, -1, readerror);
    token = strtok(nullptr, " ");
    int imat = ReadInteger(token, -1, readerror);
    token = strtok(nullptr, " ");
    token = strtok(nullptr, " ");
    token = strtok(nullptr, " ");
    token = strtok(nullptr, " ");
    std::vector<int> inode;
    for (size_t k = 0; k < 8; ++k) {
      token = strtok(nullptr, " ");
      const int in = ReadInteger(token, -1, readerror);
      if (!readerror) inode.push_back(in);
    }

    if (inode.size() != 8) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error reading file " << elist << " (line " << il
                << ").\n"
                << "    Read " << inode.size() << " node indices for element"
                << ielem << " (expected 8).\n";
      felist.close();
      return false;
    }
    // Check synchronisation.
    if (ielem - 1 != (int)m_elements.size() + nbackground) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Synchronisation lost on file " << elist << " (line "
                << il << ").\n"
                << "    Element: " << ielem << " (expected " 
                << m_elements.size() << ").\n";
      ok = false;
      break;
    }
    // Check the material number and ensure that epsilon is non-negative
    if (imat < 1 || imat > (int)m_materials.size()) {
      std::cerr << m_className << "::Initialise:\n"
                << "   Out-of-range material number on file " << elist
                << " (line " << il << ").\n"
                << "    Element: " << ielem << ", material: " << imat << ".\n";
      ok = false;
      break;
    }
    if (m_materials[imat - 1].eps < 0) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Element " << ielem << " in " << elist << "\n"
                << "    uses material " << imat << " which does not have\n" 
                << "    a positive permittivity in " << mplist << ".\n";
      ok = false;
      break;
    }
    // Check the node numbers.
    for (size_t k = 0; k < 8; ++k) {
      if (inode[k] < 1) {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Found a node number < 1 in " << elist << " (line "
                  << il << ").\n"
                  << "    Element: " << ielem << ", material: " << imat << "\n";
        ok = false;
        break;
      }
      if (inode[k] > highestnode) highestnode = inode[k];
    }
    // Skip quadrilaterals which are background.
    if (m_deleteBackground && m_materials[imat - 1].ohm == 0) {
      nbackground++;
      continue;
    }
    // Store the element, degeneracy
    Element element;
    if (inode[2] == inode[3] && inode[3] == inode[6]) {
      ndegenerate++;
      element.degenerate = true;
    } else {
      element.degenerate = false;
    }
    // Store the material reference
    element.matmap = imat - 1;
    // Node references
    if (element.degenerate) {
      element.emap[0] = inode[0] - 1;
      element.emap[1] = inode[1] - 1;
      element.emap[2] = inode[2] - 1;
      element.emap[3] = inode[4] - 1;
      element.emap[4] = inode[7] - 1;
      element.emap[5] = inode[5] - 1;
      element.emap[6] = inode[3] - 1;
      element.emap[7] = inode[6] - 1;
    } else {
      for (size_t k = 0; k < 8; ++k) element.emap[k] = inode[k] - 1;
    }
    m_elements.push_back(std::move(element));
  }
  // Close the file
  felist.close();
  if (!ok) return false;

  if (m_elements.empty()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Found no valid elements in file " << elist << ".\n";
    return false;
  } 

  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read " << m_elements.size() << " elements from file "
            << elist << ",\n";
  std::cout << "    highest node number: " << highestnode << ",\n";
  std::cout << "    degenerate elements: " << ndegenerate << ",\n";
  std::cout << "    background elements skipped: " << nbackground << ".\n";
  // Check the value of the unit.
  double funit = ScalingFactor(unit);
  if (funit <= 0.) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Unknown length unit " << unit << ". Will use cm.\n";
    funit = 1.0;
  }
  if (m_debug) {
    std::cout << m_className << "::Initialise: Unit scaling factor = "
              << funit << ".\n";
  }

  // Open the node list.
  std::ifstream fnlist(nlist);
  if (!fnlist) {
    PrintCouldNotOpen("Initialise", nlist);
    return false;
  }
  // Read the node list.
  il = 0;
  while (fnlist.getline(line, size, '\n')) {
    il++;
    // Skip page feed in batch page title.
    if (strstr(line, "VERSION") != nullptr) {
      for (size_t k = 0; k < 2; ++k) fnlist.getline(line, size, '\n');
      il += 2;
      continue;
    }
     // Split the line in tokens.
    char* token = strtok(line, " ");
    // Skip blank lines and headers.
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "LIST") == 0 || strcmp(token, "NODE") == 0 ||
        strcmp(token, "ANSYS") == 0 || strcmp(token, "***") == 0 ||
        strcmp(token, "FILE") == 0 || strcmp(token, "Electric") == 0 ||
        strcmp(token, "VERSION") == 0) {
      continue;
    }
    // Read the element.
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(nullptr, " ");
    double xnode = ReadDouble(token, -1, readerror);
    token = strtok(nullptr, " ");
    double ynode = ReadDouble(token, -1, readerror);
    token = strtok(nullptr, " ");
    double znode = ReadDouble(token, -1, readerror);
    // Check syntax.
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error reading file " << nlist << " (line " << il
                << ").\n";
      fnlist.close();
      return false;
    }
    // Check synchronisation.
    if (inode - 1 != (int)m_nodes.size()) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Synchronisation lost on file " << nlist << " (line "
                << il << ").\n"
                << "    Node: " << inode << " (expected " << m_nodes.size()
                << "), (x,y,z) = (" << xnode << ", " << ynode << ", " << znode
                << ")\n";
      ok = false;
      break;
    }
    Node node;
    node.w.clear();
    // Store the point coordinates
    node.x = xnode * funit;
    node.y = ynode * funit;
    node.z = znode * funit;
    m_nodes.push_back(std::move(node));
  }
  // Close the file
  fnlist.close();
  if (!ok) return false;

  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read " << m_nodes.size() << " nodes.\n";
  // Check number of nodes
  if ((int)m_nodes.size() != highestnode) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Number of nodes read (" << m_nodes.size() 
              << ") on " << nlist << "\n"
              << "    does not match element list (" << highestnode << ").\n";
    return false;
  }

  // Open the voltage list
  std::ifstream fprnsol(prnsol);
  if (!fprnsol) {
    PrintCouldNotOpen("Initialise", prnsol);
    return false;
  }
  // Read the voltage list
  il = 0;
  unsigned int nread = 0;
  readerror = false;
  while (fprnsol.getline(line, size, '\n')) {
    il++;
    // Skip page feed in batch page title
    if (strstr(line, "VERSION") != nullptr) {
      for (size_t k = 0; k < 2; ++k) fprnsol.getline(line, size, '\n');
      il += 2;
      continue;
    }
    // Split the line in tokens.
    char* token = strtok(line, " ");
    // Skip blank lines and headers.
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "PRINT") == 0 || strcmp(token, "ANSYS") == 0 || 
        strcmp(token, "VERSION") == 0 || strcmp(token, "NODAL") == 0 || 
        strcmp(token, "FILE") == 0 || strcmp(token, "*****") == 0 ||
        strcmp(token, "***") == 0 || strcmp(token, "LOAD") == 0 ||
        strcmp(token, "TIME=") == 0 || strcmp(token, "MAXIMUM") == 0 ||
        strcmp(token, "VALUE") == 0 || strcmp(token, "NODE") == 0) {
      continue;
    }
    // Read the element.
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(nullptr, " ");
    double volt = ReadDouble(token, -1, readerror);
    // Check syntax.
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error reading file " << prnsol << " (line " << il
                << ").\n";
      fprnsol.close();
      return false;
    }
    // Check node number and store if OK.
    if (inode < 1 || inode > (int)m_nodes.size()) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Node number " << inode << " out of range\n"
                << "    on potential file " << prnsol << " (line " << il
                << ").\n";
      ok = false;
      break;
    } else {
      m_nodes[inode - 1].v = volt;
      nread++;
    }
  }
  // Close the file
  fprnsol.close();
  if (!ok) return false;

  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read " << nread << " potentials from file " << prnsol
            << ".\n";
  // Check number of nodes.
  if (nread != m_nodes.size()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Number of nodes read (" << nread << ") on potential file "
              << prnsol << " does not\n"
              << "    match the node list (" << m_nodes.size() << ").\n";
    return false;
  }
  // Set the ready flag.
  m_ready = true;
  Prepare();
  return true;
}

bool ComponentAnsys121::SetWeightingField(const std::string& prnsol,
                                          const std::string& label) {
  if (!m_ready) {
    PrintNotReady("SetWeightingField");
    std::cerr << "    Weighting field cannot be added.\n";
    return false;
  }

  // Open the voltage list.
  std::ifstream fprnsol(prnsol);
  if (!fprnsol) {
    PrintCouldNotOpen("SetWeightingField", prnsol);
    return false;
  }

  // Check if a weighting field with the same label already exists.
  const size_t iw = GetOrCreateWeightingFieldIndex(label);
  if (iw + 1 != m_wfields.size()) {
    std::cout << m_className << "::SetWeightingField:\n"
              << "    Replacing existing weighting field " << label << ".\n";
  }

  // Buffer for reading
  constexpr int size = 100;
  char line[size];

  bool ok = true;
  // Read the voltage list.
  int il = 0;
  unsigned int nread = 0;
  bool readerror = false;
  while (fprnsol.getline(line, size, '\n')) {
    il++;
    // Skip page feed in batch page title
    if (strstr(line, "VERSION") != nullptr) {
      for (size_t k = 0; k < 2; ++k) fprnsol.getline(line, size, '\n');
      il += 2;
      continue;
    }
    // Split the line in tokens.
    char* token = strtok(line, " ");
    // Skip blank lines and headers.
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "PRINT") == 0 || strcmp(token, "*****") == 0 ||
        strcmp(token, "LOAD") == 0 || strcmp(token, "TIME=") == 0 ||
        strcmp(token, "MAXIMUM") == 0 || strcmp(token, "VALUE") == 0 ||
        strcmp(token, "NODE") == 0)
      continue;
    // Read the element.
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(nullptr, " ");
    double volt = ReadDouble(token, -1, readerror);
    // Check the syntax.
    if (readerror) {
      std::cerr << m_className << "::SetWeightingField:\n"
                << "    Error reading file " << prnsol << " (line " << il
                << ").\n";
      fprnsol.close();
      return false;
    }
    // Check node number and store if OK.
    if (inode < 1 || inode > (int)m_nodes.size()) {
      std::cerr << m_className << "::SetWeightingField:\n"
                << "    Node number " << inode << " out of range\n"
                << "    on potential file " << prnsol << " (line " << il
                << ").\n";
      ok = false;
      break;
    } else {
      m_nodes[inode - 1].w[iw] = volt;
      nread++;
    }
  }
  // Close the file.
  fprnsol.close();
  if (!ok) return false;

  // Tell how many lines read.
  std::cout << m_className << "::SetWeightingField:\n"
            << "    Read " << nread << " potentials from file " << prnsol
            << ".\n";
  // Check the number of nodes.
  if (nread != m_nodes.size()) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Number of nodes read from potential file " << prnsol
              << " (" << nread << ")\n    does not match the node list (" 
              << m_nodes.size() << ").\n";
    return false;
  }

  // Set the ready flag.
  m_wfieldsOk[iw] = true;
  return true;
}

void ComponentAnsys121::SetRangeZ(const double zmin, const double zmax) {
  if (fabs(zmax - zmin) <= 0.) {
    std::cerr << m_className << "::SetRangeZ: Zero range is not permitted.\n";
    return;
  }
  m_minBoundingBox[2] = m_mapmin[2] = std::min(zmin, zmax);
  m_maxBoundingBox[2] = m_mapmax[2] = std::max(zmin, zmax);
}

}
