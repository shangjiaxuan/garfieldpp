#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

#include "Garfield/ComponentAnsys121.hh"

namespace Garfield {

ComponentAnsys121::ComponentAnsys121() : ComponentFieldMap("Ansys121") {
  m_is3d = false;
  // Default bounding box
  m_minBoundingBox[2] = -50;
  m_maxBoundingBox[2] = 50;
}

bool ComponentAnsys121::Initialise(std::string elist, std::string nlist,
                                   std::string mplist, std::string prnsol,
                                   std::string unit) {
  m_ready = false;
  m_warning = false;
  m_nWarnings = 0;
  // Keep track of the success.
  bool ok = true;

  // Buffer for reading
  constexpr int size = 100;
  char line[size];

  // Open the material list.
  std::ifstream fmplist;
  fmplist.open(mplist.c_str(), std::ios::in);
  if (fmplist.fail()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Could not open material file " << mplist
              << " for reading. The file perhaps does not exist.\n";
    return false;
  }

  // Read the material list.
  int il = 0;
  unsigned int icurrmat = 0;
  bool readerror = false;
  while (fmplist.getline(line, size, '\n')) {
    il++;
    // Skip page feed
    if (strcmp(line, "1") == 0) {
      fmplist.getline(line, size, '\n');
      il++;
      fmplist.getline(line, size, '\n');
      il++;
      fmplist.getline(line, size, '\n');
      il++;
      fmplist.getline(line, size, '\n');
      il++;
      fmplist.getline(line, size, '\n');
      il++;
      continue;
    }
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        strcmp(token, "TEMPERATURE") == 0 || strcmp(token, "PROPERTY=") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13)
      continue;
    // Read number of materials and initialise the list.
    if (strcmp(token, "LIST") == 0) {
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      const unsigned int nMaterials = ReadInteger(token, -1, readerror);
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        ok = false;
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
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      const int imat = ReadInteger(token, -1, readerror);
      if (readerror || imat < 0) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        ok = false;
        return false;
      }
      icurrmat = imat;
    } else if (strcmp(token, "TEMP") == 0) {
      // Version 12 format: read property tag and value
      token = strtok(NULL, " ");
      int itype = 0;
      if (strncmp(token, "PERX", 4) == 0) {
        itype = 1;
      } else if (strncmp(token, "RSVX", 4) == 0) {
        itype = 2;
      } else {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Found unknown material property flag " << token
                  << "\n";
        std::cerr << "    on material properties file " << mplist << " (line "
                  << il << ").\n";
        ok = false;
      }
      fmplist.getline(line, size, '\n');
      il++;
      token = NULL;
      token = strtok(line, " ");
      if (icurrmat < 1 || icurrmat > m_materials.size()) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Found out-of-range current material index "
                  << icurrmat << "\n";
        std::cerr << "    in material properties file " << mplist << ".\n";
        ok = false;
        readerror = false;
      } else if (itype == 1) {
        m_materials[icurrmat - 1].eps = ReadDouble(token, -1, readerror);
      } else if (itype == 2) {
        m_materials[icurrmat - 1].ohm = ReadDouble(token, -1, readerror);
      }
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        ok = false;
        return false;
      }
    } else if (strcmp(token, "PROPERTY") == 0) {
      // Version 11 format
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      int itype = 0;
      if (strcmp(token, "PERX") == 0) {
        itype = 1;
      } else if (strcmp(token, "RSVX") == 0) {
        itype = 2;
      } else {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Found unknown material property flag " << token
                  << "\n";
        std::cerr << "    on material properties file " << mplist << " (line "
                  << il << ").\n";
        ok = false;
      }
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      int imat = ReadInteger(token, -1, readerror);
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        ok = false;
        return false;
      } else if (imat < 1 || imat > (int)m_materials.size()) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Found out-of-range current material index " << imat
                  << "\n";
        std::cerr << "    in material properties file " << mplist << ".\n";
        ok = false;
      } else {
        fmplist.getline(line, size, '\n');
        il++;
        fmplist.getline(line, size, '\n');
        il++;
        token = NULL;
        token = strtok(line, " ");
        token = strtok(NULL, " ");
        if (itype == 1) {
          m_materials[imat - 1].eps = ReadDouble(token, -1, readerror);
        } else if (itype == 2) {
          m_materials[imat - 1].ohm = ReadDouble(token, -1, readerror);
        }
        if (readerror) {
          std::cerr << m_className << "::Initialise:\n";
          std::cerr << "    Error reading file " << mplist << " (line " << il
                    << ").\n";
          fmplist.close();
          ok = false;
          return false;
        }
      }
    }
  }

  // Close the file
  fmplist.close();

  // Find the lowest epsilon, check for eps = 0, set default drift media
  double epsmin = -1;
  unsigned int iepsmin = 0;
  for (unsigned int imat = 0; imat < m_materials.size(); ++imat) {
    if (m_materials[imat].eps < 0) continue;
    if (m_materials[imat].eps == 0) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Material " << imat
                << " has been assigned a permittivity\n";
      std::cerr << "    equal to zero in " << mplist << ".\n";
      ok = false;
    } else if (epsmin < 0. || epsmin > m_materials[imat].eps) {
      epsmin = m_materials[imat].eps;
      iepsmin = imat;
    }
  }

  if (epsmin < 0.) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    No material with positive permittivity found \n";
    std::cerr << "    in material list " << mplist << ".\n";
    ok = false;
  } else {
    for (unsigned int imat = 0; imat < m_materials.size(); ++imat) {
      if (imat == iepsmin) {
        m_materials[imat].driftmedium = true;
      } else {
        m_materials[imat].driftmedium = false;
      }
    }
  }

  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read properties of " << m_materials.size()
            << " materials from file " << mplist << ".\n";
  if (m_debug) PrintMaterials();

  // Open the element list
  std::ifstream felist;
  felist.open(elist.c_str(), std::ios::in);
  if (felist.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open element file " << elist
              << " for reading.\n";
    std::cerr << "    The file perhaps does not exist.\n";
    return false;
  }
  // Read the element list
  m_elements.clear();
  int ndegenerate = 0;
  int nbackground = 0;
  il = 0;
  int highestnode = 0;
  while (felist.getline(line, size, '\n')) {
    il++;
    // Skip page feed in batch page title
    if (strstr(line,"VERSION") != NULL) {
      felist.getline(line, size, '\n');
      il++;
      felist.getline(line, size, '\n');
      il++;
      continue;
    }
    // Split the line in tokens
    char* token = NULL;
    // Split into tokens
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "LIST") == 0 || strcmp(token, "ELEM") == 0 ||
        strcmp(token, "ANSYS") == 0 || strcmp(token, "***") == 0 ||
        strcmp(token, "VERSION") == 0) {
      continue;
    }
    // Read the element
    int ielem = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int imat = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    if (!token) std::cerr << "No token available\n";
    int in0 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in1 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in2 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in3 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in4 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in5 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in6 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in7 = ReadInteger(token, -1, readerror);

    // Check synchronisation
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << elist << " (line " << il
                << ").\n";
      felist.close();
      ok = false;
      return false;
    } else if (ielem - 1 != (int)m_elements.size() + nbackground) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Synchronisation lost on file " << elist << " (line "
                << il << ").\n";
      std::cerr << "    Element: " << ielem << " (expected " 
                << m_elements.size() << "), material: " << imat 
                << ", nodes: (" << in0 << ", " << in1
                << ", " << in2 << ", " << in3 << ", " << in4 << ", " << in5
                << ", " << in6 << ", " << in7 << ")\n";
      ok = false;
    }
    // Check the material number and ensure that epsilon is non-negative
    if (imat < 1 || imat > (int)m_materials.size()) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "   Out-of-range material number on file " << elist
                << " (line " << il << ").\n";
      std::cerr << "    Element: " << ielem << ", material: " << imat
                << ", nodes: (" << in0 << ", " << in1 << ", " << in2 << ", "
                << in3 << ", " << in4 << ", " << in5 << ", " << in6 << ", "
                << in7 << ")\n";
      ok = false;
    }
    if (m_materials[imat - 1].eps < 0) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Element " << ielem << " in element list " << elist
                << " uses material " << imat << " which\n"
                << "    has not been assigned a positive permittivity\n";
      std::cerr << "    in material list " << mplist << ".\n";
      ok = false;
    }
    // Check the node numbers
    if (in0 < 1 || in1 < 1 || in2 < 1 || in3 < 1 || in4 < 1 || in5 < 1 ||
        in6 < 1 || in7 < 1) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Found a node number < 1 on file " << elist << " (line "
                << il << ").\n";
      std::cerr << "    Element: " << ielem << ", material: " << imat << "\n";
      ok = false;
    }
    if (in0 > highestnode) highestnode = in0;
    if (in1 > highestnode) highestnode = in1;
    if (in2 > highestnode) highestnode = in2;
    if (in3 > highestnode) highestnode = in3;
    if (in4 > highestnode) highestnode = in4;
    if (in5 > highestnode) highestnode = in5;
    if (in6 > highestnode) highestnode = in6;
    if (in7 > highestnode) highestnode = in7;
    // Skip quadrilaterals which are background.
    if (m_deleteBackground && m_materials[imat - 1].ohm == 0) {
      nbackground++;
      continue;
    }
    // Store the element, degeneracy
    Element newElement;
    if (in2 == in3 && in3 == in6) {
      ndegenerate++;
      newElement.degenerate = true;
    } else {
      newElement.degenerate = false;
    }
    // Store the material reference
    newElement.matmap = imat - 1;
    // Node references
    if (newElement.degenerate) {
      newElement.emap[0] = in0 - 1;
      newElement.emap[1] = in1 - 1;
      newElement.emap[2] = in2 - 1;
      newElement.emap[3] = in4 - 1;
      newElement.emap[4] = in7 - 1;
      newElement.emap[5] = in5 - 1;
      newElement.emap[6] = in3 - 1;
      newElement.emap[7] = in6 - 1;
    } else {
      newElement.emap[0] = in0 - 1;
      newElement.emap[1] = in1 - 1;
      newElement.emap[2] = in2 - 1;
      newElement.emap[3] = in3 - 1;
      newElement.emap[4] = in4 - 1;
      newElement.emap[5] = in5 - 1;
      newElement.emap[6] = in6 - 1;
      newElement.emap[7] = in7 - 1;
    }
    m_elements.push_back(std::move(newElement));
  }
  // Close the file
  felist.close();
  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read " << m_elements.size() << " elements from file "
            << elist << ",\n";
  std::cout << "    highest node number: " << highestnode << ",\n";
  std::cout << "    degenerate elements: " << ndegenerate << ",\n";
  std::cout << "    background elements skipped: " << nbackground << ".\n";
  // Check the value of the unit
  double funit = ScalingFactor(unit);
  if (funit <= 0.) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Unknown length unit " << unit << ".\n";
    ok = false;
    funit = 1.0;
  }
  if (m_debug) {
    std::cout << m_className << "::Initialise: Unit scaling factor = "
              << funit << ".\n";
  }

  // Open the node list
  std::ifstream fnlist;
  fnlist.open(nlist.c_str(), std::ios::in);
  if (fnlist.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open nodes file " << nlist << " for reading.\n";
    std::cerr << "    The file perhaps does not exist.\n";
    return false;
  }
  // Read the node list
  m_nodes.clear();
  il = 0;
  while (fnlist.getline(line, size, '\n')) {
    il++;
    // Skip page feed in batch page title
    if (strstr(line,"VERSION") != NULL) {
      fnlist.getline(line, size, '\n');
      il++;
      fnlist.getline(line, size, '\n');
      il++;
      continue;
    }
     // Split the line in tokens
    char* token = NULL;
    // Split into tokens
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "LIST") == 0 || strcmp(token, "NODE") == 0 ||
        strcmp(token, "ANSYS") == 0 || strcmp(token, "***") == 0 ||
        strcmp(token, "FILE") == 0 || strcmp(token, "Electric") == 0 ||
        strcmp(token, "VERSION") == 0) {
      continue;
    }
    // Read the element
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    double xnode = ReadDouble(token, -1, readerror);
    token = strtok(NULL, " ");
    double ynode = ReadDouble(token, -1, readerror);
    token = strtok(NULL, " ");
    double znode = ReadDouble(token, -1, readerror);
    // Check syntax
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << nlist << " (line " << il
                << ").\n";
      fnlist.close();
      ok = false;
      return false;
    }
    // Check synchronisation
    if (inode - 1 != (int)m_nodes.size()) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Synchronisation lost on file " << nlist << " (line "
                << il << ").\n";
      std::cerr << "    Node: " << inode << " (expected " << m_nodes.size()
                << "), (x,y,z) = (" << xnode << ", " << ynode << ", " << znode
                << ")\n";
      ok = false;
    }
    Node newNode;
    newNode.w.clear();
    // Store the point coordinates
    newNode.x = xnode * funit;
    newNode.y = ynode * funit;
    newNode.z = znode * funit;
    m_nodes.push_back(std::move(newNode));
  }
  // Close the file
  fnlist.close();
  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read " << m_nodes.size() << " nodes.\n";
  // Check number of nodes
  if ((int)m_nodes.size() != highestnode) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Number of nodes read (" << m_nodes.size() 
              << ") on " << nlist << "\n"
              << "    does not match element list (" << highestnode << ").\n";
    ok = false;
  }

  // Open the voltage list
  std::ifstream fprnsol;
  fprnsol.open(prnsol.c_str(), std::ios::in);
  if (fprnsol.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open potential file " << prnsol
              << " for reading.\n";
    std::cerr << "    The file perhaps does not exist.\n";
    return false;
  }
  // Read the voltage list
  il = 0;
  unsigned int nread = 0;
  readerror = false;
  while (fprnsol.getline(line, size, '\n')) {
    il++;
    // Skip page feed in batch page title
    if (strstr(line,"VERSION") != NULL) {
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      continue;
    }
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
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
    // Read the element
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    double volt = ReadDouble(token, -1, readerror);
    // Check syntax
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << prnsol << " (line " << il
                << ").\n";
      fprnsol.close();
      ok = false;
      return false;
    }
    // Check node number and store if OK
    if (inode < 1 || inode > (int)m_nodes.size()) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Node number " << inode << " out of range\n";
      std::cerr << "    on potential file " << prnsol << " (line " << il
                << ").\n";
      ok = false;
    } else {
      m_nodes[inode - 1].v = volt;
      nread++;
    }
  }
  // Close the file
  fprnsol.close();
  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n";
  std::cout << "    Read " << nread << " potentials from file " << prnsol
            << ".\n";
  // Check number of nodes
  if (nread != m_nodes.size()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Number of nodes read (" << nread << ") on potential file "
              << prnsol << " does not\n"
              << "    match the node list (" << m_nodes.size() << ").\n";
    ok = false;
  }
  // Set the ready flag
  if (ok) {
    m_ready = true;
  } else {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr
        << "    Field map could not be read and cannot be interpolated.\n";
    return false;
  }

  // Remove weighting fields (if any).
  m_wfields.clear();
  m_wfieldsOk.clear();

  // Establish the ranges
  SetRange();
  UpdatePeriodicity();
  return true;
}

bool ComponentAnsys121::SetWeightingField(std::string prnsol,
                                          std::string label) {
  if (!m_ready) {
    PrintNotReady("SetWeightingField");
    std::cerr << "    Weighting field cannot be added.\n";
    return false;
  }

  // Open the voltage list.
  std::ifstream fprnsol;
  fprnsol.open(prnsol.c_str(), std::ios::in);
  if (fprnsol.fail()) {
    std::cerr << m_className << "::SetWeightingField:\n";
    std::cerr << "    Could not open potential file " << prnsol
              << " for reading.\n";
    std::cerr << "    The file perhaps does not exist.\n";
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
    if (strstr(line,"VERSION") != NULL) {
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      continue;
    }
    // Split the line in tokens.
    char* token = NULL;
    token = strtok(line, " ");
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
    token = strtok(NULL, " ");
    double volt = ReadDouble(token, -1, readerror);
    // Check the syntax.
    if (readerror) {
      std::cerr << m_className << "::SetWeightingField:\n";
      std::cerr << "    Error reading file " << prnsol << " (line " << il
                << ").\n";
      fprnsol.close();
      return false;
    }
    // Check node number and store if OK.
    if (inode < 1 || inode > (int)m_nodes.size()) {
      std::cerr << m_className << "::SetWeightingField:\n";
      std::cerr << "    Node number " << inode << " out of range\n";
      std::cerr << "    on potential file " << prnsol << " (line " << il
                << ").\n";
      ok = false;
    } else {
      m_nodes[inode - 1].w[iw] = volt;
      nread++;
    }
  }
  // Close the file.
  fprnsol.close();
  // Tell how many lines read.
  std::cout << m_className << "::SetWeightingField:\n";
  std::cout << "    Read " << nread << " potentials from file " << prnsol
            << ".\n";
  // Check the number of nodes.
  if (nread != m_nodes.size()) {
    std::cerr << m_className << "::SetWeightingField:\n";
    std::cerr << "    Number of nodes read (" << nread << ") "
              << "    on potential file " << prnsol << "\n";
    std::cerr << "    does not match the node list (" << m_nodes.size() << ").\n";
    ok = false;
  }

  // Set the ready flag.
  m_wfieldsOk[iw] = ok;
  if (!ok) {
    std::cerr << m_className << "::SetWeightingField:\n";
    std::cerr << "    Field map could not be read "
              << "and cannot be interpolated.\n";
    return false;
  }
  return true;
}

void ComponentAnsys121::ElectricField(const double x, const double y,
                                      const double z, double& ex, double& ey,
                                      double& ez, Medium*& m, int& status) {
  double v;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

void ComponentAnsys121::ElectricField(const double xin, const double yin,
                                      const double zin, double& ex, double& ey,
                                      double& ez, double& volt, Medium*& m,
                                      int& status) {
  // Copy the coordinates
  double x = xin, y = yin, z = 0.;

  // Map the coordinates onto field map coordinates
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  // Initial values
  ex = ey = ez = volt = 0;
  m = nullptr;
  status = 0;

  // Do not proceed if not properly initialised.
  if (!m_ready) {
    status = -10;
    PrintNotReady("ElectricField");
    return;
  }

  if (m_warning) PrintWarning("ElectricField");

  if (zin < m_minBoundingBox[2] || zin > m_maxBoundingBox[2]) {
    status = -5;
    return;
  }

  // Find the element that contains this point
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement5(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField:\n";
      std::cout << "    Point (" << x << ", " << y << ") not in the mesh.\n";
    }
    status = -6;
    return;
  }

  const Element& element = m_elements[imap];
  if (m_debug) {
    PrintElement("ElectricField", x, y, z, t1, t2, t3, t4, element, 8);
  }

  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];
  // Calculate quadrilateral field, which can degenerate to a triangular field
  const double invdet = 1. / det;
  if (element.degenerate) {
    volt = n0.v * t1 * (2 * t1 - 1) + n1.v * t2 * (2 * t2 - 1) +
           n2.v * t3 * (2 * t3 - 1) + 4 * n3.v * t1 * t2 + 4 * n4.v * t1 * t3 +
           4 * n5.v * t2 * t3;
    ex = -(n0.v * (4 * t1 - 1) * jac[0][1] + n1.v * (4 * t2 - 1) * jac[1][1] +
           n2.v * (4 * t3 - 1) * jac[2][1] +
           n3.v * (4 * t2 * jac[0][1] + 4 * t1 * jac[1][1]) +
           n4.v * (4 * t3 * jac[0][1] + 4 * t1 * jac[2][1]) +
           n5.v * (4 * t3 * jac[1][1] + 4 * t2 * jac[2][1])) *
         invdet;
    ey = -(n0.v * (4 * t1 - 1) * jac[0][2] + n1.v * (4 * t2 - 1) * jac[1][2] +
           n2.v * (4 * t3 - 1) * jac[2][2] +
           n3.v * (4 * t2 * jac[0][2] + 4 * t1 * jac[1][2]) +
           n4.v * (4 * t3 * jac[0][2] + 4 * t1 * jac[2][2]) +
           n5.v * (4 * t3 * jac[1][2] + 4 * t2 * jac[2][2])) *
         invdet;
  } else {
    const Node& n6 = m_nodes[element.emap[6]];
    const Node& n7 = m_nodes[element.emap[7]];
    volt = -n0.v * (1 - t1) * (1 - t2) * (1 + t1 + t2) * 0.25 -
           n1.v * (1 + t1) * (1 - t2) * (1 - t1 + t2) * 0.25 -
           n2.v * (1 + t1) * (1 + t2) * (1 - t1 - t2) * 0.25 -
           n3.v * (1 - t1) * (1 + t2) * (1 + t1 - t2) * 0.25 +
           n4.v * (1 - t1) * (1 + t1) * (1 - t2) * 0.5 +
           n5.v * (1 + t1) * (1 + t2) * (1 - t2) * 0.5 +
           n6.v * (1 - t1) * (1 + t1) * (1 + t2) * 0.5 +
           n7.v * (1 - t1) * (1 + t2) * (1 - t2) * 0.5;
    ex = -(n0.v * ((1 - t2) * (2 * t1 + t2) * jac[0][0] +
                   (1 - t1) * (t1 + 2 * t2) * jac[1][0]) *
               0.25 +
           n1.v * ((1 - t2) * (2 * t1 - t2) * jac[0][0] -
                   (1 + t1) * (t1 - 2 * t2) * jac[1][0]) *
               0.25 +
           n2.v * ((1 + t2) * (2 * t1 + t2) * jac[0][0] +
                   (1 + t1) * (t1 + 2 * t2) * jac[1][0]) *
               0.25 +
           n3.v * ((1 + t2) * (2 * t1 - t2) * jac[0][0] -
                   (1 - t1) * (t1 - 2 * t2) * jac[1][0]) *
               0.25 +
           n4.v * (t1 * (t2 - 1) * jac[0][0] +
                   (t1 - 1) * (t1 + 1) * jac[1][0] * 0.5) +
           n5.v * ((1 - t2) * (1 + t2) * jac[0][0] * 0.5 -
                   (1 + t1) * t2 * jac[1][0]) +
           n6.v * (-t1 * (1 + t2) * jac[0][0] +
                   (1 - t1) * (1 + t1) * jac[1][0] * 0.5) +
           n7.v * ((t2 - 1) * (t2 + 1) * jac[0][0] * 0.5 +
                   (t1 - 1) * t2 * jac[1][0])) *
         invdet;
    ey = -(n0.v * ((1 - t2) * (2 * t1 + t2) * jac[0][1] +
                   (1 - t1) * (t1 + 2 * t2) * jac[1][1]) *
               0.25 +
           n1.v * ((1 - t2) * (2 * t1 - t2) * jac[0][1] -
                   (1 + t1) * (t1 - 2 * t2) * jac[1][1]) *
               0.25 +
           n2.v * ((1 + t2) * (2 * t1 + t2) * jac[0][1] +
                   (1 + t1) * (t1 + 2 * t2) * jac[1][1]) *
               0.25 +
           n3.v * ((1 + t2) * (2 * t1 - t2) * jac[0][1] -
                   (1 - t1) * (t1 - 2 * t2) * jac[1][1]) *
               0.25 +
           n4.v * (t1 * (t2 - 1) * jac[0][1] +
                   (t1 - 1) * (t1 + 1) * jac[1][1] * 0.5) +
           n5.v * ((1 - t2) * (1 + t2) * jac[0][1] * 0.5 -
                   (1 + t1) * t2 * jac[1][1]) +
           n6.v * (-t1 * (1 + t2) * jac[0][1] +
                   (1 - t1) * (1 + t1) * jac[1][1] * 0.5) +
           n7.v * ((t2 - 1) * (t2 + 1) * jac[0][1] * 0.5 +
                   (t1 - 1) * t2 * jac[1][1])) *
         invdet;
  }

  // Transform field to global coordinates
  UnmapFields(ex, ey, ez, x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  // Drift medium?
  if (m_debug) {
    std::cout << m_className << "::ElectricField:\n";
    std::cout << "    Material " << element.matmap << ", drift flag "
              << m_materials[element.matmap].driftmedium << ".\n";
  }
  m = m_materials[element.matmap].medium;
  status = -5;
  if (m_materials[element.matmap].driftmedium) {
    if (m && m->IsDriftable()) status = 0;
  }
}

void ComponentAnsys121::WeightingField(const double xin, const double yin,
                                       const double zin, double& wx, double& wy,
                                       double& wz, const std::string& label) {
  // Initial values
  wx = wy = wz = 0;

  // Do not proceed if not properly initialised.
  if (!m_ready) return;

  // Look for the label.
  const int iw = GetWeightingFieldIndex(label);
  // Do not proceed if the requested weighting field does not exist.
  if (iw < 0) return;
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
  const int imap = FindElement5(x, y, z, t1, t2, t3, t4, jac, det);
  // Check if the point is in the mesh.
  if (imap < 0) return;

  const Element& element = m_elements[imap];
  if (m_debug) {
    PrintElement("WeightingField", x, y, z, t1, t2, t3, t4, element, 8, iw);
  }
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];
  // Calculate quadrilateral field, which can degenerate to a triangular field
  const double invdet = 1. / det;
  if (m_elements[imap].degenerate) {
    wx = -(n0.w[iw] * (4 * t1 - 1) * jac[0][1] +
           n1.w[iw] * (4 * t2 - 1) * jac[1][1] +
           n2.w[iw] * (4 * t3 - 1) * jac[2][1] +
           n3.w[iw] * (4 * t2 * jac[0][1] + 4 * t1 * jac[1][1]) +
           n4.w[iw] * (4 * t3 * jac[0][1] + 4 * t1 * jac[2][1]) +
           n5.w[iw] * (4 * t3 * jac[1][1] + 4 * t2 * jac[2][1])) *
         invdet;
    wy = -(n0.w[iw] * (4 * t1 - 1) * jac[0][2] +
           n1.w[iw] * (4 * t2 - 1) * jac[1][2] +
           n2.w[iw] * (4 * t3 - 1) * jac[2][2] +
           n3.w[iw] * (4 * t2 * jac[0][2] + 4 * t1 * jac[1][2]) +
           n4.w[iw] * (4 * t3 * jac[0][2] + 4 * t1 * jac[2][2]) +
           n5.w[iw] * (4 * t3 * jac[1][2] + 4 * t2 * jac[2][2])) *
         invdet;
  } else {
    const Node& n6 = m_nodes[element.emap[6]];
    const Node& n7 = m_nodes[element.emap[7]];
    wx = -(n0.w[iw] * ((1 - t2) * (2 * t1 + t2) * jac[0][0] +
                       (1 - t1) * (t1 + 2 * t2) * jac[1][0]) *
               0.25 +
           n1.w[iw] * ((1 - t2) * (2 * t1 - t2) * jac[0][0] -
                       (1 + t1) * (t1 - 2 * t2) * jac[1][0]) *
               0.25 +
           n2.w[iw] * ((1 + t2) * (2 * t1 + t2) * jac[0][0] +
                       (1 + t1) * (t1 + 2 * t2) * jac[1][0]) *
               0.25 +
           n3.w[iw] * ((1 + t2) * (2 * t1 - t2) * jac[0][0] -
                       (1 - t1) * (t1 - 2 * t2) * jac[1][0]) *
               0.25 +
           n4.w[iw] * (t1 * (t2 - 1) * jac[0][0] +
                       (t1 - 1) * (t1 + 1) * jac[1][0] * 0.5) +
           n5.w[iw] * ((1 - t2) * (1 + t2) * jac[0][0] * 0.5 -
                       (1 + t1) * t2 * jac[1][0]) +
           n6.w[iw] * (-t1 * (1 + t2) * jac[0][0] +
                       (1 - t1) * (1 + t1) * jac[1][0] * 0.5) +
           n7.w[iw] * ((t2 - 1) * (1 + t2) * jac[0][0] * 0.5 +
                       (t1 - 1) * t2 * jac[1][0])) *
         invdet;
    wy = -(n0.w[iw] * ((1 - t2) * (2 * t1 + t2) * jac[0][1] +
                       (1 - t1) * (t1 + 2 * t2) * jac[1][1]) *
               0.25 +
           n1.w[iw] * ((1 - t2) * (2 * t1 - t2) * jac[0][1] -
                       (1 + t1) * (t1 - 2 * t2) * jac[1][1]) *
               0.25 +
           n2.w[iw] * ((1 + t2) * (2 * t1 + t2) * jac[0][1] +
                       (1 + t1) * (t1 + 2 * t2) * jac[1][1]) *
               0.25 +
           n3.w[iw] * ((1 + t2) * (2 * t1 - t2) * jac[0][1] -
                       (1 - t1) * (t1 - 2 * t2) * jac[1][1]) *
               0.25 +
           n4.w[iw] * (t1 * (t2 - 1) * jac[0][1] +
                       (t1 - 1) * (t1 + 1) * jac[1][1] * 0.5) +
           n5.w[iw] * ((1 - t2) * (1 + t2) * jac[0][1] * 0.5 -
                       (1 + t1) * t2 * jac[1][1]) +
           n6.w[iw] * (-t1 * (1 + t2) * jac[0][1] +
                       (1 - t1) * (1 + t1) * jac[1][1] * 0.5) +
           n7.w[iw] * ((t2 - 1) * (t2 + 1) * jac[0][1] * 0.5 +
                       (t1 - 1) * t2 * jac[1][1])) *
         invdet;
  }

  // Transform field to global coordinates
  UnmapFields(wx, wy, wz, x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);
}

double ComponentAnsys121::WeightingPotential(const double xin, const double yin,
                                             const double zin,
                                             const std::string& label) {
  // Do not proceed if not properly initialised.
  if (!m_ready) return 0.;

  // Look for the label.
  const int iw = GetWeightingFieldIndex(label);
  // Do not proceed if the requested weighting field does not exist.
  if (iw < 0) return 0.;
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
  const int imap = FindElement5(x, y, z, t1, t2, t3, t4, jac, det);
  // Check if the point is in the mesh
  if (imap < 0) return 0.;

  const Element& element = m_elements[imap];
  if (m_debug) {
    PrintElement("WeightingPotential", x, y, z, t1, t2, t3, t4, element, 8, iw);
  }
  // Calculate quadrilateral field, which can degenerate to a triangular field
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const Node& n4 = m_nodes[element.emap[4]];
  const Node& n5 = m_nodes[element.emap[5]];
  if (element.degenerate) {
    return n0.w[iw] * t1 * (2 * t1 - 1) + n1.w[iw] * t2 * (2 * t2 - 1) +
           n2.w[iw] * t3 * (2 * t3 - 1) + 4 * n3.w[iw] * t1 * t2 +
           4 * n4.w[iw] * t1 * t3 + 4 * n5.w[iw] * t2 * t3;
  }

  const Node& n6 = m_nodes[element.emap[6]];
  const Node& n7 = m_nodes[element.emap[7]];
  return -n0.w[iw] * (1 - t1) * (1 - t2) * (1 + t1 + t2) * 0.25 -
         n1.w[iw] * (1 + t1) * (1 - t2) * (1 - t1 + t2) * 0.25 -
         n2.w[iw] * (1 + t1) * (1 + t2) * (1 - t1 - t2) * 0.25 -
         n3.w[iw] * (1 - t1) * (1 + t2) * (1 + t1 - t2) * 0.25 +
         n4.w[iw] * (1 - t1) * (1 + t1) * (1 - t2) * 0.5 +
         n5.w[iw] * (1 + t1) * (1 + t2) * (1 - t2) * 0.5 +
         n6.w[iw] * (1 - t1) * (1 + t1) * (1 + t2) * 0.5 +
         n7.w[iw] * (1 - t1) * (1 + t2) * (1 - t2) * 0.5;
}

Medium* ComponentAnsys121::GetMedium(const double xin, const double yin,
                                     const double zin) {
  // Copy the coordinates.
  double x = xin, y = yin, z = 0.;

  // Map the coordinates onto field map coordinates.
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  if (zin < m_minBoundingBox[2] || z > m_maxBoundingBox[2]) {
    return nullptr;
  }

  // Do not proceed if not properly initialised.
  if (!m_ready) {
    PrintNotReady("GetMedium");
    return nullptr;
  }
  if (m_warning) PrintWarning("GetMedium");

  // Find the element that contains this point.
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement5(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) {
    if (m_debug) {
      std::cerr << m_className << "::GetMedium:\n";
      std::cerr << "    Point (" << x << ", " << y << ") not in the mesh.\n";
    }
    return nullptr;
  }
  const Element& element = m_elements[imap];
  if (element.matmap >= m_materials.size()) {
    if (m_debug) {
      std::cerr << m_className << "::GetMedium:\n";
      std::cerr << "    Point (" << x << ", " << y << ")"
                << " has out of range material number " << imap << ".\n";
    }
    return nullptr;
  }

  if (m_debug) {
    PrintElement("GetMedium", x, y, z, t1, t2, t3, t4, element, 8);
  }

  // Assign a medium.
  return m_materials[element.matmap].medium;
}

void ComponentAnsys121::SetRangeZ(const double zmin, const double zmax) {
  if (fabs(zmax - zmin) <= 0.) {
    std::cerr << m_className << "::SetRangeZ: Zero range is not permitted.\n";
    return;
  }
  m_minBoundingBox[2] = m_mapmin[2] = std::min(zmin, zmax);
  m_maxBoundingBox[2] = m_mapmax[2] = std::max(zmin, zmax);
}

void ComponentAnsys121::UpdatePeriodicity() {
  UpdatePeriodicity2d();
  UpdatePeriodicityCommon();
}

double ComponentAnsys121::GetElementVolume(const unsigned int i) {
  if (i >= m_elements.size()) return 0.;
  const Element& element = m_elements[i];
  const Node& n0 = m_nodes[element.emap[0]];
  const Node& n1 = m_nodes[element.emap[1]];
  const Node& n2 = m_nodes[element.emap[2]];
  const Node& n3 = m_nodes[element.emap[3]];
  const double surf =
      0.5 *
      (fabs((n1.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n1.y - n0.y)) +
       fabs((n3.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n3.y - n0.y)));
  return surf;
}

void ComponentAnsys121::GetAspectRatio(const unsigned int i, double& dmin,
                                       double& dmax) {
  if (i >= m_elements.size()) {
    dmin = dmax = 0.;
    return;
  }

  const Element& element = m_elements[i];
  const int np = 8;
  // Loop over all pairs of vertices.
  for (int j = 0; j < np - 1; ++j) {
    const Node& nj = m_nodes[element.emap[j]];
    for (int k = j + 1; k < np; ++k) {
      const Node& nk = m_nodes[element.emap[k]];
      // Compute distance.
      const double dx = nj.x - nk.x;
      const double dy = nj.y - nk.y;
      const double dist = sqrt(dx * dx + dy * dy);
      if (k == 1) {
        dmin = dmax = dist;
      } else {
        if (dist < dmin) dmin = dist;
        if (dist > dmax) dmax = dist;
      }
    }
  }
}
}
