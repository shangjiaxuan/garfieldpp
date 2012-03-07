// Copied and modified ComponentAnsys123.cc
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>

#include "ComponentCST.hh"

namespace Garfield {

ComponentCST::ComponentCST() : ComponentFieldMap() {

  className = "ComponentCST";
  ready = false;
  // Default bounding box
  zMinBoundingBox = -50.;
  zMaxBoundingBox =  50.;
  m_xlines.clear();
  m_ylines.clear();
  m_zlines.clear();

}

bool
ComponentCST::Initialise(std::string elist, std::string nlist,
                         std::string mplist, std::string prnsol,
                         std::string unit) {
  ready = false;

  // Keep track of the success
  bool ok = true;

  // Buffer for reading
  const int size = 200;
  char line[size];
  // Open the material list
  std::ifstream fmplist;
  fmplist.open(mplist.c_str(), std::ios::in);
  if (fmplist.fail()) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Could not open material file " << mplist << " for reading.\n",
    std::cerr << "    The file perhaps does not exist.\n";
    return false;
  }

  // Read the material list
  nMaterials = 0;
  int il = 0;
  bool readerror = false;
  while (fmplist.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
         int(token[0]) == 10 || int(token[0]) == 13) continue;
    // Read number of materials,
    // ensure it does not exceed the maximum and initialise the list
    if (strcmp(token, "Materials") == 0) {
      token = strtok(NULL, " ");
      nMaterials = ReadInteger(token, -1, readerror);
      if (readerror) {
        std::cerr << className << "::Initialise:\n";
        std::cerr << "    Error reading file " <<  mplist 
                  << " (line " << il << ").\n";
        fmplist.close();
        ok = false;
        return false;
      }
      materials.resize(nMaterials);
      for (int i = 0; i < nMaterials; ++i) {
        materials[i].ohm = -1;
        materials[i].eps = -1;
        materials[i].medium = NULL;
      }
      if (debug) {
    	std::cout << className << "::Initialise:\n";
        std::cout << "    Number of materials: " << nMaterials << "\n";
      }
    } else if (strcmp(token, "Material") == 0) {
      token = strtok(NULL, " ");
      int imat = ReadInteger(token, -1, readerror);
      if (readerror) {
    	std::cerr << className << "::Initialise:\n";
    	std::cerr << "     Error reading file " << mplist
                  << " (line " << il << ".\n";
        fmplist.close();
        ok = false;
        return false;
      } else if (imat < 1 || imat > nMaterials) {
    	std::cerr << className << "::Initialise:\n";
    	std::cerr << "    Found out-of-range material index " << imat << "in\n";
        std::cerr << "    material properties file " << mplist << ".\n";
        ok = false;
      } else {
        token = strtok(NULL, " ");
        int itype = 0;
        if (strcmp(token,"PERX") == 0) {
          itype = 1;
        } else if (strcmp(token,"RSVX") == 0) {
          itype = 2;
        } else {
          std::cerr << className << "::Initialise:\n";
          std::cerr << "    Found unknown material property flag " << token << "\n";
          std::cerr << "    on material properties file " << mplist
                    << "(line " << il << ").\n";
          ok = false;
        }
        token = strtok(NULL, " ");
        if (itype == 1) {
          materials[imat - 1].eps = ReadDouble(token, -1, readerror);
        } else if (itype == 2) {
          materials[imat - 1].ohm = ReadDouble(token, -1, readerror);
          token = strtok(NULL, " ");
          if (strcmp(token,"PERX") != 0){
            std::cerr << className << "::Initialise:\n";
            std::cerr << "   Found unknown material property falg "<< token << "\n";
            std::cerr << "   on material file " << mplist
                      << " (material " << imat << ").\n)";
            ok = false;
          } else {
            token = strtok(NULL, " ");
            materials[imat - 1].eps = ReadDouble(token, -1, readerror);
          }
        }
        if (readerror) {
          std::cerr << className << "::Initialise:\n";
          std::cerr << "     Error reading file " << mplist << "(line " << il << ").\n";
          fmplist.close();
          ok = false;
          return false;
        }
        if (debug){
          std::cout << className << "::Initialise:\n";
          std::cout << "    Read material properties for material " << (imat-1) << "\n";
          if (itype == 2) {
            std::cout << "    eps = " << materials[imat - 1].eps 
                      << " ohm = " << materials[imat - 1].ohm << "\n";
          } else {
            std::cout << "    eps = " << materials[imat - 1].eps << "\n";
          }
        }
      }
    }
  }
  // Close the file
  fmplist.close();
  // Find the lowest epsilon, check for eps = 0, set default drift media
  double epsmin = -1; int iepsmin = -1;
  for (int imat = 0; imat < nMaterials; ++imat) {
    if (materials[imat].eps < 0) continue;
    if (materials[imat].eps == 0) {
      std::cout << className << "::Initialise:\n";
      std::cout << "    Material " << imat << " has been assigned a permittivity\n";
      std::cout << "    equal to zero in " << mplist << ".\n";
      ok = false;
    } else if (iepsmin < 0 || epsmin > materials[imat].eps) {
      epsmin = materials[imat].eps;
      iepsmin = imat;
    }
  }
  if (iepsmin < 0) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "     No material with positive permittivity found in\n";
    std::cerr << "     material list " << mplist.c_str() << ".\n";
    ok = false;
  } else {
    for (int imat = 0; imat < nMaterials; ++imat) {
      if (imat == iepsmin) {
        materials[imat].driftmedium = true;
      } else {
        materials[imat].driftmedium = false;
      }
    }
  }
  // Tell how many lines read
  std::cout << className << "::Initialise:\n";
  std::cout << "    Read properties of " << nMaterials << " materials\n";
  std::cout << "    from file " << mplist << ".\n";
  if (debug) PrintMaterials();

  // Check the value of the unit
  double funit;
  if (strcmp(unit.c_str(),"mum") == 0 || strcmp(unit.c_str(),"micron") == 0 ||
      strcmp(unit.c_str(),"micrometer") == 0) {
    funit = 0.0001;
  } else if (strcmp(unit.c_str(),"mm") == 0 ||
             strcmp(unit.c_str(),"millimeter") == 0) {
    funit = 0.1;
  } else if (strcmp(unit.c_str(),"cm") == 0 ||
             strcmp(unit.c_str(),"centimeter") == 0) {
    funit = 1.0;
  } else if (strcmp(unit.c_str(),"m") == 0 ||
             strcmp(unit.c_str(),"meter") == 0) {
    funit = 100.0;
  } else {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Unknown length unit " << unit << ".\n";
    ok = false;
    funit = 1.0;
  }
  if (debug) {
    std::cout << className << "::Initialise:\n";
    std::cout << "    Unit scaling factor = " << funit << ".\n";
  }

  // Open the node list
  std::ifstream fnlist;
  fnlist.open(nlist.c_str(), std::ios::in);
  if (fnlist.fail()) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Could not open nodes file " << nlist << " for reading.\n";
    std::cerr << "    The file perhaps does not exist.\n";
    return false;
  }
  // Read the node list
  nodes.clear();
  nNodes = 0;
  il = 0;
  int xlines = 0, ylines = 0, zlines = 0;
  int lines_type = -1;
  double line_tmp;
  while (fnlist.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    // Split into tokens
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13) continue;
    // Read max sizes
    if (strcmp(token, "xmax") == 0) {
      token = strtok(NULL, " "); xlines = ReadInteger(token, -1, readerror);
      token = strtok(NULL, " ");
      token = strtok(NULL, " "); ylines = ReadInteger(token, -1, readerror);
      token = strtok(NULL, " ");
      token = strtok(NULL, " "); zlines = ReadInteger(token, -1, readerror);
      if (readerror) break;
      continue;
    }
    if (strcmp(token, "x-lines\n") == 0 || strcmp(token, "x-lines") == 0) {
      lines_type = 1;
      if (debug) {
        std::cout << className << "::Initialise:\n";
        std::cout << "    Reading x-lines from file  " << nlist << ".\n";
      }
      continue;
    }
    if (strcmp(token, "y-lines\n") == 0 || strcmp(token, "y-lines") == 0) {
      lines_type = 2;
      if (debug) {
        std::cout << className << "::Initialise:\n";
        std::cout << "    Reading y-lines from file  " << nlist << ".\n";
      }
      continue;
    }
    if (strcmp(token, "z-lines\n") == 0 || strcmp(token, "z-lines") == 0) {
      lines_type = 3;
      if (debug) {
        std::cout << className << "::Initialise:\n";
        std::cout << "    Reading z-lines from file  " << nlist << ".\n";
      }
      continue;
    }
    line_tmp = ReadDouble(token, -1, readerror);
    if (lines_type == 1) m_xlines.push_back(line_tmp * funit);
    else if (lines_type == 2) m_ylines.push_back(line_tmp * funit);
    else if (lines_type == 3) m_zlines.push_back(line_tmp * funit);
    else {
      std::cerr << className << "::Initialise:\n";
      std::cerr << "    Line type was not set in  " << nlist
                << " (line " << il << ", token = " << token << ").\n";
      std::cerr << "    Maybe it is in the wrong format\n";
      std::cerr << "    e.g. missing tailing space after x-lines.\n";
      ok = false;
      break;
    }
    if (readerror) break;
  }
  // Check syntax
  if (readerror) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Error reading file " << nlist << " (line " << il << ").\n";
    fnlist.close();
    ok = false;
    return false;
  }
  // Close the file
  fnlist.close();
  // Calculate the node positions
  for (unsigned int z = 0; z < m_zlines.size(); z++) {
    for (unsigned int y = 0; y < m_ylines.size(); y++){
      for (unsigned int x = 0; x < m_xlines.size(); x++){
        node newNode;
        // Store the point coordinates
        newNode.x = m_xlines.at(x);
        newNode.y = m_ylines.at(y);
        newNode.z = m_zlines.at(z);
        nodes.push_back(newNode);
        ++nNodes;
      }
    }
  }
  if ((unsigned)xlines == m_xlines.size() &&
      (unsigned)ylines == m_ylines.size() &&
      (unsigned)zlines == m_zlines.size()) {
    std::cout << className << "::Initialise:\n";
    std::cout << "    Found in file " << nlist << "\n    "
                              << xlines << " x-lines\n    "
                              << ylines << " y-lines\n    "
                              << zlines << " z-lines\n";
  } else {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    There should be " << xlines << " x-lines, "
                                        << ylines << " y-lines and "
                                        << zlines << " z-lines in file "
                                        << nlist << " but I found :\n    "
                                        << m_xlines.size() << " x-lines, "
                                        << m_ylines.size() << " x-lines, "
                                        << m_zlines.size() << " z-lines.\n";
  }
  // Check synchronisation
  if ((xlines * ylines * zlines) != nNodes) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Synchronisation lost on file " << nlist << ".\n";
    std::cerr << "    Nodes: " << nNodes << " (expected "
              << (xlines * ylines * zlines) << ")\n";
    ok = false;
  }

  // Tell how many lines read
  std::cout << className << "::Initialise:\n";
  std::cout << "    Read " << nNodes << " nodes from file " << nlist << ".\n";
  // Check number of nodes

   // Open the element list
  std::ifstream felist;
  felist.open(elist.c_str(), std::ios::in);
  if (felist.fail()) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Could not open element file " << elist << " for reading.\n";
    std::cerr << "    The file perhaps does not exist.\n";
    return false;
  }
  // Read the element list
  elements.clear();
  nElements = 0;
  element newElement;
  int ndegenerate = 0;
  int nbackground = 0;
  il = 0;
  int highestnode = 0;
  while (felist.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    // Split into tokens
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token," ") == 0 || strcmp(token,"\n") == 0 ||
         int(token[0]) == 10 || int(token[0]) == 13 ||
         strcmp(token,"LIST") == 0 || strcmp(token,"ELEM") == 0) continue;
    // Read the element
    int ielem = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " "); int imat = ReadInteger(token, -1, readerror);
    // construct node numbers
    std::vector<int> nodes;
    GetNodesForElement(ielem,nodes);
    int in0 = nodes.at(0); int in1 = nodes.at(1); int in2 = nodes.at(2);
    int in3 = nodes.at(3); int in4 = nodes.at(4); int in5 = nodes.at(5);
    int in6 = nodes.at(6); int in7 = nodes.at(7);

    // Check synchronisation
    if (readerror) {
      std::cerr << className << "::Initialise:\n";
      std::cerr << "    Error reading file " << elist << " (line " << il << ").\n";
      felist.close();
      ok = false;
      return false;
    } else if (ielem != nElements + nbackground) {
      std::cerr << className << "::Initialise:\n";
      std::cerr << "    Synchronisation lost on file " << elist << " (line " << il << ").\n";
      std::cerr << "    Element: " << ielem << " (expected " << nElements << "), material: " << imat
                << ", nodes: ("
                << in0 << " "
                << in1 << " "
                << in2 << " "
                << in3 << " "
                << in4 << " "
                << in5 << " "
                << in6 << " "
                << in7 << ")\n";
      ok = false;
    }
    // Check the material number and ensure that epsilon is non-negative
    if (imat < 1 || imat > nMaterials) {
      std::cerr << className << "::Initialise:\n";
      std::cerr << "   Out-of-range material number on file " << elist 
                << " (line " << il << ").\n";
      std::cerr << "    Element: " << ielem << ", material: " << imat
                << ", nodes: ("
                << in0 << " "
                << in1 << " "
                << in2 << " "
                << in3 << " "
                << in4 << " "
                << in5 << " "
                << in6 << " "
                << in7 << ")\n";
      ok = false;
    }
    if (materials[imat - 1].eps < 0) {
      std::cerr << className << "::Initialise:\n";
      std::cerr << "    Element " << ielem << " in element list " << elist 
                << " uses material " << imat << " which\n";
      std::cerr << "    has not been assigned a positive permittivity\n";
      std::cerr << "    in material list " << mplist << ".\n";
      ok = false;
    }
     // Check the node numbers
    if (in0 < 0 || in1 < 0 || in2 < 0 || in3 < 0 ||
        in4 < 0 || in5 < 0 || in6 < 0 || in7 < 0) {
      std::cerr << className << "::Initialise:\n";
      std::cerr << "    Found a node number < 0 on file " << elist 
                << " (line " << il << ").\n";
      std::cerr << "    Element: " << ielem << ", material: " << imat << ",\n";
      std::cerr << "    nodes: ("
                << in0 << " "
                << in1 << " "
                << in2 << " "
                << in3 << " "
                << in4 << " "
                << in5 << " "
                << in6 << " "
                << in7 <<")\n";
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
    if (deleteBackground && materials[imat - 1].ohm == 0) {
      nbackground++;
      continue;
    }
     // No degenerated elements in CST
    newElement.degenerate = false;
    // Store the material reference
    newElement.matmap = imat - 1;
     // Node references
    newElement.emap[0] = in0;
    newElement.emap[1] = in1;
    newElement.emap[2] = in2;
    newElement.emap[3] = in3;
    newElement.emap[4] = in4;
    newElement.emap[5] = in5;
    newElement.emap[6] = in6;
    newElement.emap[7] = in7;
    elements.push_back(newElement);
    // if (nElements > 934100) {
    //   std::cout << className << "::Initialise\n";
    //   std::cout << "    Element: " << nElements << "\n";
    //   std::cout << "    Nodes  : " 
    //             << in0 << ", " << in1 << ", " << in2 << ", " << in3 << ", "
    //             << in4 << ", " << in5 << ", " << in6 << ", " << in7 << "\n";
     // }
    ++nElements;
  }
  // Close the file
  felist.close();
  // Tell how many lines read
  std::cout << className << "::Initialise:\n";
  std::cout << "    Read " << nElements << " elements from file " << elist << ",\n";
  std::cout << "    highest node number: " << highestnode << ",\n";
  std::cout << "    degenerate elements: " << ndegenerate << ",\n";
  std::cout << "    background elements skipped: " << nbackground << ".\n";

  if (nNodes != (highestnode+1)) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Number of nodes read (" << nNodes << ") on " << nlist << " \n";
    std::cerr << "    does not match element list (" << highestnode << ").\n";
    std::cerr << "    Maybe the line size exceeded 200 characters.\n";
    ok = false;
  }
  // Open the voltage list
  std::ifstream fprnsol;
  fprnsol.open(prnsol.c_str(), std::ios::in);
  if (fprnsol.fail()) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Could not open potential file " << prnsol << " for reading.\n";
    std::cerr << "    The file perhaps does not exist.\n";
    return false;
  }
  // Read the voltage list
  il = 0;
  int nread = 0;
  readerror = false;
  while (fprnsol.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token," ") == 0 || strcmp(token,"\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token,"Max") == 0) continue;
    // Read the element
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " "); double volt = ReadDouble(token, -1, readerror);
    // Check syntax
    if (readerror) {
      std::cerr << className << "::Initialise:\n";
      std::cerr << "    Error reading file " << prnsol << " (line " << il << ").\n";
      fprnsol.close();
      ok = false;
      return false;
    }
    // Check node number and store if OK
    if (inode < 1 || inode > nNodes) {
      std::cerr << className << "::Initialise:\n";
      std::cerr << "    Node number " << inode << " out of range\n";
      std::cerr << "    on potential file " << prnsol << " (line " << il << ").\n";
      ok = false;
    } else {
      nodes[inode - 1].v = volt;
      nread++;
    }
  }
  // Close the file
  fprnsol.close();
  // Tell how many lines read
  std::cout << className << "::Initialise:\n";
  std::cout << "    Read " << nread << " potentials from file " << prnsol << ".\n";
  // Check number of nodes
  if (nread != nNodes) {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Number of nodes read (" << nread 
              << ") on potential file " << prnsol << " does not\n";
    std::cerr << "    match the node list (" << nNodes << ").\n";
    ok = false;
  }
  // Set the ready flag
  if (ok) {
    ready = true;
  } else {
    std::cerr << className << "::Initialise:\n";
    std::cerr << "    Field map could not be read and cannot be interpolated.\n";
    return false;
  }

  // Establish the ranges
  SetRange();
  UpdatePeriodicity();
  return true;

}

bool
ComponentCST::SetWeightingField(std::string prnsol, std::string label) {

  if (!ready) {
    std::cerr << className << "::SetWeightingField:\n";
    std::cerr << "    No valid field map is present.\n";
    std::cerr << "    Weighting field cannot be added.\n";
    return false;
  }

  // Open the voltage list
  std::ifstream fprnsol;
  fprnsol.open(prnsol.c_str(), std::ios::in);
  if (fprnsol.fail()) {
    std::cerr << className << "::SetWeightingField:\n";
    std::cerr << "    Could not open potential file " << prnsol << " for reading.\n";
    std::cerr << "    The file perhaps does not exist.\n";
    return false;
  }
  // Check if a weighting field with the same label already exists.
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
    for (int j = nNodes; j--;) {
      nodes[j].w.resize(nWeightingFields);
    }
  } else {
    std::cout << className << "::SetWeightingField:\n";
    std::cout << "    Replacing existing weighting field " << label << ".\n";
  }
  wfields[iw] = label;
  wfieldsOk[iw] = false;

  // Buffer for reading
  const int size = 100;
  char line[size];

  bool ok = true;
  // Read the voltage list
  int il = 0;
  int nread = 0;
  bool readerror = false;
  while (fprnsol.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token," ") == 0 || strcmp(token,"\n") == 0 ||
         int(token[0]) == 10 || int(token[0]) == 13 ||
         strcmp(token,"PRINT")   == 0 || strcmp(token,"*****")   == 0 ||
         strcmp(token,"LOAD")    == 0 || strcmp(token,"TIME=")   == 0 ||
         strcmp(token,"MAXIMUM") == 0 || strcmp(token,"VALUE")   == 0 ||
         strcmp(token,"NODE")    == 0) continue;
    // Read the element
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    double volt = ReadDouble(token, -1, readerror);
    // Check syntax
    if (readerror) {
      std::cerr << className << "::SetWeightingField:\n";
      std::cerr << "    Error reading file " << prnsol
                << " (line " << il << ").\n";
      fprnsol.close();
      return false;
    }
    // Check node number and store if OK
    if (inode < 1 || inode > nNodes) {
      std::cerr << className << "::SetWeightingField:\n";
      std::cerr << "    Node number " << inode << " out of range.\n";
      std::cerr << "    on potential file " << prnsol
                << " (line " << il << ").\n";
      ok = false;
    } else {
      nodes[inode - 1].w[iw] = volt;
      nread++;
    }
  }
  // Close the file
  fprnsol.close();
  // Tell how many lines read
  std::cout << className << "::SetWeightingField:\n";
  std::cout << "    Read " << nread << " potentials from file " << prnsol << ".\n";
  // Check number of nodes
  if (nread != nNodes) {
    std::cerr << className << "::SetWeightingField:\n";
    std::cerr << "    Number of nodes read (" << nread << ")"
              << " on potential file "  << prnsol << "\n";
    std::cerr << "     does not match the node list (" << nNodes << ".\n";
    ok = false;
  }

  // Set the ready flag.
  wfieldsOk[iw] = ok;
  if (!ok) {
    std::cerr << className << "::SetWeightingField:\n";
    std::cerr << "    Field map could not be read "
              << "and cannot be interpolated.\n";
    return false;
  }

  return true;

}

void
ComponentCST::ElectricField(const double x, const double y, const double z,
double& ex, double& ey, double& ez,Medium*& m, int& status) {

  double v;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);

}

void
ComponentCST::ElectricField(
        const double xin, const double yin, const double zin,
        double& ex, double& ey, double& ez, double& volt,
        Medium*& m, int& status) {

  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirrored, ymirrored, zmirrored;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z,
  xmirrored, ymirrored, zmirrored,
  rcoordinate, rotation);

  // Initial values
  ex = ey = ez = volt = 0;
  status = 0;

  // Do not proceed if not properly initialised.
  if (!ready) {
    status = -10;
    std::cerr << className << "::ElectricField:\n";
    std::cerr << "    Field map not available for interpolation.\n";
    return;
  }

  if (warning) {
    std::cout << className << "::ElectricField:\n";
    std::cout << "    Warnings have been issued for this field map.\n";
  }
  double t1, t2, t3, jac[3][3], det;

  int imap = FindElementCube(x, y, z, t1, t2, t3, jac, det);

  if (imap < 0) {
    if (debug) {
      std::cout << className << "::ElectricField:\n";
      std::cout << "    Point (" << x << "," << y << "," << z 
                << ") not in the mesh,\n"; 
      std::cout << "    it is background or PEC.\n";
    }
    status = -6;
    return;
  }
  // Save element number of last element
  lastElement = imap;

  double inv_jac[3][3];
  det = fabs(det);
  inv_jac[0][0] = (1./det) * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]);
  inv_jac[0][1] = (1./det) * (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]);
  inv_jac[0][2] = (1./det) * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);
  inv_jac[1][0] = (1./det) * (jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]);
  inv_jac[1][1] = (1./det) * (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]);
  inv_jac[1][2] = (1./det) * (jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]);
  inv_jac[2][0] = (1./det) * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
  inv_jac[2][1] = (1./det) * (jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1]);
  inv_jac[2][2] = (1./det) * (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]);
  // Field calculation
  volt = (nodes[elements[imap].emap[0]].v * (1 - t1) * (1 - t2) * (1 - t3) +
          nodes[elements[imap].emap[1]].v * (1 + t1) * (1 - t2) * (1 - t3) +
          nodes[elements[imap].emap[2]].v * (1 + t1) * (1 + t2) * (1 - t3) +
          nodes[elements[imap].emap[3]].v * (1 - t1) * (1 + t2) * (1 - t3) +
          nodes[elements[imap].emap[4]].v * (1 - t1) * (1 - t2) * (1 + t3) +
          nodes[elements[imap].emap[5]].v * (1 + t1) * (1 - t2) * (1 + t3) +
          nodes[elements[imap].emap[6]].v * (1 + t1) * (1 + t2) * (1 + t3) +
          nodes[elements[imap].emap[7]].v * (1 - t1) * (1 + t2) * (1 + t3)) / 8.;
  ex = (1. / 8.) * (
          nodes[elements[imap].emap[0]].v * (-1 * (1 - t2) * (1 - t3) * inv_jac[0][0] +
                                             (1 - t1) * -1 * (1 - t3) * inv_jac[1][0] +
                                             (1 - t1) * (1 - t2) * -1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[1]].v * (+1 * (1 - t2) * (1 - t3) * inv_jac[0][0] +
                                             (1 + t1) * -1 * (1 - t3) * inv_jac[1][0] +
                                             (1 + t1) * (1 - t2) * -1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[2]].v * (+1 * (1 + t2) * (1 - t3) * inv_jac[0][0] +
                                             (1 + t1) * +1 * (1 - t3) * inv_jac[1][0] +
                                             (1 + t1) * (1 + t2) * -1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[3]].v * (-1 * (1 + t2) * (1 - t3) * inv_jac[0][0] +
                                             (1 - t1) * +1 * (1 - t3) * inv_jac[1][0] +
                                             (1 - t1) * (1 + t2) * -1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[4]].v * (-1 * (1 - t2) * (1 + t3) * inv_jac[0][0] +
                                             (1 - t1) * -1 * (1 + t3) * inv_jac[1][0] +
                                             (1 - t1) * (1 - t2) * +1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[5]].v * (+1 * (1 - t2) * (1 + t3) * inv_jac[0][0] +
                                             (1 + t1) * -1 * (1 + t3) * inv_jac[1][0] +
                                             (1 + t1) * (1 - t2) * +1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[6]].v * (+1 * (1 + t2) * (1 + t3) * inv_jac[0][0] +
                                             (1 + t1) * +1 * (1 + t3) * inv_jac[1][0] +
                                             (1 + t1) * (1 + t2) * +1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[7]].v * (-1 * (1 + t2) * (1 + t3) * inv_jac[0][0] +
                                             (1 - t1) * +1 * (1 + t3) * inv_jac[1][0] +
                                             (1 - t1) * (1 + t2) * +1 * inv_jac[2][0]));
  ey = (1. / 8.) * (
          nodes[elements[imap].emap[0]].v * (-1 * (1 - t2) * (1 - t3) * inv_jac[0][1] +
                                             (1 - t1) * -1 * (1 - t3) * inv_jac[1][1] +
                                             (1 - t1) * (1 - t2) * -1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[1]].v * (+1 * (1 - t2) * (1 - t3) * inv_jac[0][1] +
                                             (1 + t1) * -1 * (1 - t3) * inv_jac[1][1] +
                                             (1 + t1) * (1 - t2) * -1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[2]].v * (+1 * (1 + t2) * (1 - t3) * inv_jac[0][1] +
                                             (1 + t1) * +1 * (1 - t3) * inv_jac[1][1] +
                                             (1 + t1) * (1 + t2) * -1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[3]].v * (-1 * (1 + t2) * (1 - t3) * inv_jac[0][1] +
                                             (1 - t1) * +1 * (1 - t3) * inv_jac[1][1] +
                                             (1 - t1) * (1 + t2) * -1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[4]].v * (-1 * (1 - t2) * (1 + t3) * inv_jac[0][1] +
                                             (1 - t1) * -1 * (1 + t3) * inv_jac[1][1] +
                                             (1 - t1) * (1 - t2) * +1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[5]].v * (+1 * (1 - t2) * (1 + t3) * inv_jac[0][1] +
                                             (1 + t1) * -1 * (1 + t3) * inv_jac[1][1] +
                                             (1 + t1) * (1 - t2) * +1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[6]].v * (+1 * (1 + t2) * (1 + t3) * inv_jac[0][1] +
                                             (1 + t1) * +1 * (1 + t3) * inv_jac[1][1] +
                                             (1 + t1) * (1 + t2) * +1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[7]].v * (-1 * (1 + t2) * (1 + t3) * inv_jac[0][1] +
                                             (1 - t1) * +1 * (1 + t3) * inv_jac[1][1] +
                                             (1 - t1) * (1 + t2) * +1 * inv_jac[2][1]));
  ez = (-1. / 8.) * (
          nodes[elements[imap].emap[0]].v * (-1 * (1 - t2) * (1 - t3) * inv_jac[0][2] +
                                             (1 - t1) * -1 * (1 - t3) * inv_jac[1][2] +
                                             (1 - t1) * (1 - t2) * -1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[1]].v * (+1 * (1 - t2) * (1 - t3) * inv_jac[0][2] +
                                             (1 + t1) * -1 * (1 - t3) * inv_jac[1][2] +
                                             (1 + t1) * (1 - t2) * -1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[2]].v * (+1 * (1 + t2) * (1 - t3) * inv_jac[0][2] +
                                             (1 + t1) * +1 * (1 - t3) * inv_jac[1][2] +
                                             (1 + t1) * (1 + t2) * -1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[3]].v * (-1 * (1 + t2) * (1 - t3) * inv_jac[0][2] +
                                             (1 - t1) * +1 * (1 - t3) * inv_jac[1][2] +
                                             (1 - t1) * (1 + t2) * -1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[4]].v * (-1 * (1 - t2) * (1 + t3) * inv_jac[0][2] +
                                             (1 - t1) * -1 * (1 + t3) * inv_jac[1][2] +
                                             (1 - t1) * (1 - t2) * +1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[5]].v * (+1 * (1 - t2) * (1 + t3) * inv_jac[0][2] +
                                             (1 + t1) * -1 * (1 + t3) * inv_jac[1][2] +
                                             (1 + t1) * (1 - t2) * +1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[6]].v * (+1 * (1 + t2) * (1 + t3) * inv_jac[0][2] +
                                             (1 + t1) * +1 * (1 + t3) * inv_jac[1][2] +
                                             (1 + t1) * (1 + t2) * +1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[7]].v * (-1 * (1 + t2) * (1 + t3) * inv_jac[0][2] +
                                             (1 - t1) * +1 * (1 + t3) * inv_jac[1][2] +
                                             (1 - t1) * (1 + t2) * +1 * inv_jac[2][2]));
  // Transform field to global coordinates
  UnmapFields(ex, ey, ez, x, y, z,
              xmirrored, ymirrored, zmirrored,
              rcoordinate, rotation);

  if (debug) {
    std::cout << className << "::ElectricField:\n";
    std::cout << "    Element number: " << imap << ".\n";
    std::cout << "    Material " << elements[imap].matmap << ", drift flag " 
              << materials[elements[imap].matmap].driftmedium << ".\n";
    std::cout << "    Local Coordinates (" << t1 << "," << t2 << "," << t3 
              << ") Voltage: " << volt << "\n";
    std::cout << "    E-Field (" << ex << "," << ey << "," << ez << ")\n";
    std::cout << "*******End of ComponentCST::ElectricField********\n\n";
  }
  // Drift medium?
  m = materials[elements[imap].matmap].medium;
  status = -5;
  if (materials[elements[imap].matmap].driftmedium) {
    if (m != 0) {
      if (m->IsDriftable()) status = 0;
    }
  }
}

void
ComponentCST::WeightingField(
        const double xin, const double yin, const double zin,
        double& wx, double& wy, double& wz,
        const std::string label) {

  // Initial values
  wx = wy = wz = 0;

  // Do not proceed if not properly initialised.
  if (!ready) return;

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

  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirrored, ymirrored, zmirrored;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z,
                 xmirrored, ymirrored, zmirrored,
                 rcoordinate, rotation);

  if (warning) {
    std::cout << className << "::WeightingField:\n";
    std::cout << "    Warnings have been issued for this field map.\n";
  }

  // Find the element that contains this point
  double t1, t2, t3, jac[3][3], det;

  int imap = FindElementCube(x, y, z, t1, t2, t3, jac, det);

  // Check if the point is in the mesh
  if (imap < 0) return;

  if (debug) {
    std::cout << className << "::WeightingField:\n";
    std::cout << "    Global: (" << x << "," << y << "," << z << "),\n";
    std::cout << "    Local: (" << t1 << "," << t2 << "," << t3 
              << ") in element " << imap << " \n";
    std::cout << "    Node xyzV\n";
    for (int i = 0; i < 8; i++) {
      std::cout << "  " << elements[imap].emap[i]
                << " " << nodes[elements[imap].emap[i]].x
                << " " << nodes[elements[imap].emap[i]].y
                << " " << nodes[elements[imap].emap[i]].z
                << " " << nodes[elements[imap].emap[i]].w[iw] << "\n";
    }
  }
  double inv_jac[3][3];
  det = fabs(det);
  inv_jac[0][0] = (1. / det) * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]);
  inv_jac[0][1] = (1. / det) * (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]);
  inv_jac[0][2] = (1. / det) * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);
  inv_jac[1][0] = (1. / det) * (jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]);
  inv_jac[1][1] = (1. / det) * (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]);
  inv_jac[1][2] = (1. / det) * (jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]);
  inv_jac[2][0] = (1. / det) * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
  inv_jac[2][1] = (1. / det) * (jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1]);
  inv_jac[2][2] = (1. / det) * (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]);

  wx = (1. / 8.) * (
          nodes[elements[imap].emap[0]].w[iw] * (-1 * (1 - t2) * (1 - t3) * inv_jac[0][0] +
                                                 (1 - t1) * -1 * (1 - t3) * inv_jac[1][0] +
                                                 (1 - t1) * (1 - t2) * -1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[1]].w[iw] * (+1 * (1 - t2) * (1 - t3) * inv_jac[0][0] +
                                                 (1 + t1) * -1 * (1 - t3) * inv_jac[1][0] +
                                                 (1 + t1) * (1 - t2) * -1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[2]].w[iw] * (+1 * (1 + t2) * (1 - t3) * inv_jac[0][0] +
                                                 (1 + t1) * +1 * (1 - t3) * inv_jac[1][0] +
                                                 (1 + t1) * (1 + t2) * -1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[3]].w[iw] * (-1 * (1 + t2) * (1 - t3) * inv_jac[0][0] +
                                                 (1 - t1) * +1 * (1 - t3) * inv_jac[1][0] +
                                                 (1 - t1) * (1 + t2) * -1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[4]].w[iw] * (-1 * (1 - t2) * (1 + t3) * inv_jac[0][0] +
                                                 (1 - t1) * -1 * (1 + t3) * inv_jac[1][0] +
                                                 (1 - t1) * (1 - t2) * +1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[5]].w[iw] * (+1 * (1 - t2) * (1 + t3) * inv_jac[0][0] +
                                                 (1 + t1) * -1 * (1 + t3) * inv_jac[1][0] +
                                                 (1 + t1) * (1 - t2) * +1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[6]].w[iw] * (+1 * (1 + t2) * (1 + t3) * inv_jac[0][0] +
                                                 (1 + t1) * +1 * (1 + t3) * inv_jac[1][0] +
                                                 (1 + t1) * (1 + t2) * +1 * inv_jac[2][0]) +
          nodes[elements[imap].emap[7]].w[iw] * (-1 * (1 + t2) * (1 + t3) * inv_jac[0][0] +
                                                 (1 - t1) * +1 * (1 + t3) * inv_jac[1][0] +
                                                 (1 - t1) * (1 + t2) * +1 * inv_jac[2][0]));
  wy = (1. / 8.) * (
          nodes[elements[imap].emap[0]].w[iw] * (-1 * (1 - t2) * (1 - t3) * inv_jac[0][1] +
                                                 (1 - t1) * -1 * (1 - t3) * inv_jac[1][1] +
                                                 (1 - t1) * (1 - t2) * -1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[1]].w[iw] * (+1 * (1 - t2) * (1 - t3) * inv_jac[0][1] +
                                                 (1 + t1) * -1 * (1 - t3) * inv_jac[1][1] +
                                                 (1 + t1) * (1 - t2) * -1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[2]].w[iw] * (+1 * (1 + t2) * (1 - t3) * inv_jac[0][1] +
                                                 (1 + t1) * +1 * (1 - t3) * inv_jac[1][1] +
                                                 (1 + t1) * (1 + t2) * -1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[3]].w[iw] * (-1 * (1 + t2) * (1 - t3) * inv_jac[0][1] +
                                                 (1 - t1) * +1 * (1 - t3) * inv_jac[1][1] +
                                                 (1 - t1) * (1 + t2) * -1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[4]].w[iw] * (-1 * (1 - t2) * (1 + t3) * inv_jac[0][1] +
                                                 (1 - t1) * -1 * (1 + t3) * inv_jac[1][1] +
                                                 (1 - t1) * (1 - t2) * +1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[5]].w[iw] * (+1 * (1 - t2) * (1 + t3) * inv_jac[0][1] +
                                                 (1 + t1) * -1 * (1 + t3) * inv_jac[1][1] +
                                                 (1 + t1) * (1 - t2) * +1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[6]].w[iw] * (+1 * (1 + t2) * (1 + t3) * inv_jac[0][1] +
                                                 (1 + t1) * +1 * (1 + t3) * inv_jac[1][1] +
                                                 (1 + t1) * (1 + t2) * +1 * inv_jac[2][1]) +
          nodes[elements[imap].emap[7]].w[iw] * (-1 * (1 + t2) * (1 + t3) * inv_jac[0][1] +
                                                 (1 - t1) * +1 * (1 + t3) * inv_jac[1][1] +
                                                 (1 - t1) * (1 + t2) * +1 * inv_jac[2][1]));
  wz = (-1. / 8.) * (
          nodes[elements[imap].emap[0]].w[iw] * (-1 * (1 - t2) * (1 - t3) * inv_jac[0][2] +
                                                 (1 - t1) * -1 * (1 - t3) * inv_jac[1][2] +
                                                 (1 - t1) * (1 - t2) * -1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[1]].w[iw] * (+1 * (1 - t2) * (1 - t3) * inv_jac[0][2] +
                                                 (1 + t1) * -1 * (1 - t3) * inv_jac[1][2] +
                                                 (1 + t1) * (1 - t2) * -1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[2]].w[iw] * (+1 * (1 + t2) * (1 - t3) * inv_jac[0][2] +
                                                 (1 + t1) * +1 * (1 - t3) * inv_jac[1][2] +
                                                 (1 + t1) * (1 + t2) * -1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[3]].w[iw] * (-1 * (1 + t2) * (1 - t3) * inv_jac[0][2] +
                                                 (1 - t1) * +1 * (1 - t3) * inv_jac[1][2] +
                                                 (1 - t1) * (1 + t2) * -1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[4]].w[iw] * (-1 * (1 - t2) * (1 + t3) * inv_jac[0][2] +
                                                 (1 - t1) * -1 * (1 + t3) * inv_jac[1][2] +
                                                 (1 - t1) * (1 - t2) * +1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[5]].w[iw] * (+1 * (1 - t2) * (1 + t3) * inv_jac[0][2] +
                                                 (1 + t1) * -1 * (1 + t3) * inv_jac[1][2] +
                                                 (1 + t1) * (1 - t2) * +1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[6]].w[iw] * (+1 * (1 + t2) * (1 + t3) * inv_jac[0][2] +
                                                 (1 + t1) * +1 * (1 + t3) * inv_jac[1][2] +
                                                 (1 + t1) * (1 + t2) * +1 * inv_jac[2][2]) +
          nodes[elements[imap].emap[7]].w[iw] * (-1 * (1 + t2) * (1 + t3) * inv_jac[0][2] +
                                                 (1 - t1) * +1 * (1 + t3) * inv_jac[1][2] +
                                                 (1 - t1) * (1 + t2) * +1 * inv_jac[2][2]));

  // Transform field to global coordinates
  UnmapFields(wx, wy, wz, x, y, z,
              xmirrored, ymirrored, zmirrored,
              rcoordinate, rotation);

}


double
ComponentCST::WeightingPotential(
  const double xin, const double yin, const double zin,
  const std::string label) {

  // Do not proceed if not properly initialised.
  if (!ready) return 0.;

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

  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirrored, ymirrored, zmirrored;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z,
                 xmirrored, ymirrored, zmirrored,
                 rcoordinate, rotation);

  if (warning) {
    std::cout << className << "::WeightingPotential:\n";
    std::cout << "Warnings have been issued for this field map.\n";
  }

  // Find the element that contains this point
  double t1, t2, t3, jac[3][3], det;

  int imap = FindElementCube(x, y, z, t1, t2, t3, jac, det);

  // Check if the point is in the mesh
  if (imap < 0) return 0.;

  if (debug) {
    std::cout << className << "::WeightingPotential:\n";
    std::cout << "    Global: (" << x << "," << y << "," << z << "),\n";
    std::cout << "    Local: (" << t1 << "," << t2 << "," << t3 
              << ") in element " << imap << "\n";
    std::cout << "  Node xyzV\n";
    for (int i = 0; i < 8; ++i) {
      std::cout << "  " << elements[imap].emap[i]
                << " " <<  nodes[elements[imap].emap[i]].x
                << " " <<  nodes[elements[imap].emap[i]].y
                << " " <<  nodes[elements[imap].emap[i]].z
                << " " <<  nodes[elements[imap].emap[i]].v << "\n";
    }
  }

  return (1. / 8.) * ( 
          nodes[elements[imap].emap[0]].w[iw] * (1 - t1) * (1 - t2) * (1 - t3) +
          nodes[elements[imap].emap[1]].w[iw] * (1 + t1) * (1 - t2) * (1 - t3) +
          nodes[elements[imap].emap[2]].w[iw] * (1 + t1) * (1 + t2) * (1 - t3) +
          nodes[elements[imap].emap[3]].w[iw] * (1 - t1) * (1 + t2) * (1 - t3) +
          nodes[elements[imap].emap[4]].w[iw] * (1 - t1) * (1 - t2) * (1 + t3) +
          nodes[elements[imap].emap[5]].w[iw] * (1 + t1) * (1 - t2) * (1 + t3) +
          nodes[elements[imap].emap[6]].w[iw] * (1 + t1) * (1 + t2) * (1 + t3) +
          nodes[elements[imap].emap[7]].w[iw] * (1 - t1) * (1 + t2) * (1 + t3));

}

bool
ComponentCST::GetMedium(const double xin, const double yin, const double zin, 
                        Medium*& m) {

  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirrored, ymirrored, zmirrored;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z,
                 xmirrored, ymirrored, zmirrored,
                 rcoordinate, rotation);

  // Initial value
  m = 0;

  // Do not proceed if not properly initialised.
  if (!ready) {
    std::cerr << className << "::GetMedium:\n";
    std::cerr << "Field map not available for interpolation.\n";
    return false;
  }
  if (warning) {
    std::cout << className << "::GetMedium:\n";
    std::cout << "Warnings have been issued for this field map.\n";
  }

  // Find the element that contains this point.
  double t1, t2, t3, jac[3][3], det;

  int imap = FindElementCube(x, y, z, t1, t2, t3, jac, det);

  if (imap < 0) {
    if (debug) {
      std::cerr << className << "::GetMedium:\n";
      std::cerr << "Point (" << x << "," << y << "," << z << ") not in the mesh.\n";
    }
    return false;
  }
  if (elements[imap].matmap < 0 || elements[imap].matmap >= nMaterials ) {
    if (debug) {
      std::cerr << className << "::GetMedium:\n";
      std::cerr << "Point (" << x << "," << y << "," << z 
                << ") has out of range material number " << imap << ".\n";
    }
    return false;
  }

  if (debug) {
    std::cout << className << "::GetMedium:\n";
    std::cout << "    Global: (" << x << "," << y << "," << z << "),\n";
    std::cout << "    Local: (" << t1 << "," << t2 << "," << t3 
              << ") in element " << imap 
              << " (degenerate: " << elements[imap].degenerate << ")\n";
    std::cout << "    Node xyzV\n";
    for (int i = 0; i < 8; i++) {
      std::cout << "     " <<  elements[imap].emap[i]
                << " " <<  nodes[elements[imap].emap[i]].x
                << " " <<  nodes[elements[imap].emap[i]].y
                << " " <<  nodes[elements[imap].emap[i]].z
                << " " <<  nodes[elements[imap].emap[i]].v << "\n";
    }
  }

  // Assign a medium.
  m = materials[elements[imap].matmap].medium;
  if (m == 0) return false;
  return true;

}

void
ComponentCST::SetRangeZ(const double zmin, const double zmax) {

  if (fabs(zmax - zmin) <= 0.) {
    std::cerr << className << "::SetRangeZ:\n";
    std::cerr << "    Zero range is not permitted.\n";
    return;
  }
  zMinBoundingBox = std::min(zmin, zmax);
  zMaxBoundingBox = std::min(zmin, zmax);
}

void
ComponentCST::UpdatePeriodicity() {

  UpdatePeriodicity2d();
  UpdatePeriodicityCommon();

}

double
  ComponentCST::GetElementVolume(const int i) {

  if (i < 0 || i >= nElements) return 0.;
  const double volume =
      fabs((nodes[elements[i].emap[1]].x - nodes[elements[i].emap[0]].x) *
           (nodes[elements[i].emap[4]].y - nodes[elements[i].emap[0]].y) *
           (nodes[elements[i].emap[3]].z - nodes[elements[i].emap[0]].z));
  return volume;
}

void
ComponentCST::GetAspectRatio(const int i, double& dmin, double& dmax) {

  if (i < 0 || i >= nElements) {
    dmin = dmax = 0.;
    return;
  }

  const int np = 8;
  // Loop over all pairs of vertices.
  for (int j = 0; j < np - 1; ++j) {
    for (int k = j + 1; k < np; ++k) {
      // Compute distance.
      const double dist = sqrt(
      pow(nodes[elements[i].emap[j]].x - nodes[elements[i].emap[k]].x, 2) +
      pow(nodes[elements[i].emap[j]].y - nodes[elements[i].emap[k]].y, 2));
      if (k == 1) {
        dmin = dist;
        dmax = dist;
      } else {
        if (dist < dmin) dmin = dist;
        if (dist > dmax) dmax = dist;
      }
    }
  }
}


void
ComponentCST::Element2Index(int element, int &i, int &j, int &k) {

  /* Here the index along x,y,z direction of the given element
   * is calculated (i,j,k).
   * i,j,k start at 0 and reach at maximum
   * m_xlines-1,m_ylines-1,m_zlines-1
   */

  int tmp = element;
  k = element / ((m_xlines.size() - 1) * (m_ylines.size() - 1));
  tmp -= k * (m_xlines.size() - 1) * (m_ylines.size() - 1);
  j = tmp / (m_xlines.size() - 1);
  i = element - j * (m_xlines.size() - 1) - k * (m_xlines.size() - 1) * (m_ylines.size() - 1);
}

void
ComponentCST::GetNodesForElement(int element, std::vector<int> &nodes){
   /*
   global coordinates   8 _ _ _ _7    t3    t2
                       /       /|     ^   /|
     ^ z              /       / |     |   /
     |             5 /_______/6 |     |  /
     |              |        |  |     | /
     |              |  4     |  /3    |/     t1
      ------->      |        | /       ------->
     /      y       |        |/       local coordinates
    /               1--------2
   /
  v x
  */

  int i,j,k;
  Element2Index(element,i,j,k);
  nodes.clear();
  nodes.push_back((i+1)+    j*m_xlines.size()+    k*m_xlines.size()*m_ylines.size());
  nodes.push_back((i+1)+(j+1)*m_xlines.size()+    k*m_xlines.size()*m_ylines.size());
  nodes.push_back(    i+(j+1)*m_xlines.size()+    k*m_xlines.size()*m_ylines.size());
  nodes.push_back(    i+    j*m_xlines.size()+    k*m_xlines.size()*m_ylines.size());
  nodes.push_back((i+1)+    j*m_xlines.size()+(k+1)*m_xlines.size()*m_ylines.size());
  nodes.push_back((i+1)+(j+1)*m_xlines.size()+(k+1)*m_xlines.size()*m_ylines.size());
  nodes.push_back(    i+(j+1)*m_xlines.size()+(k+1)*m_xlines.size()*m_ylines.size());
  nodes.push_back(    i+    j*m_xlines.size()+(k+1)*m_xlines.size()*m_ylines.size());
}

int
ComponentCST::FindElementCube(const double x, const double y, const double z,
                      double& t1, double& t2, double& t3,
                      double jac[3][3], double& det){

  int imap = -1;
  // check if the last element matches
  // if yes than leave directly
  if (lastElement >= 0 &&
      x >= nodes[elements[lastElement].emap[3]].x &&
      y >= nodes[elements[lastElement].emap[3]].y &&
      z >= nodes[elements[lastElement].emap[3]].z &&
      x < nodes[elements[lastElement].emap[0]].x &&
      y < nodes[elements[lastElement].emap[2]].y &&
      z < nodes[elements[lastElement].emap[7]].z) {
    imap = lastElement;
    CoordinatesCube(x,y,z,t1,t2,t3,jac,det,imap);
    return imap;
  }

  // CST specific element search without an element loop
  int index_x = 0, index_y = 0, index_z = 0;
  if (imap == -1) {
    std::vector<double>::iterator it;
    double my_x[] = {x};
    it = std::search(m_xlines.begin(),m_xlines.end(),my_x,my_x+1,Greater);
    index_x = std::distance(m_xlines.begin(),(it-1));
    double my_y[] = {y};
    it = std::search(m_ylines.begin(),m_ylines.end(),my_y,my_y+1,Greater);
    index_y = std::distance(m_ylines.begin(),(it-1));
    double my_z[] = {z};
    it = std::search(m_zlines.begin(),m_zlines.end(),my_z,my_z+1,Greater);
    index_z = std::distance(m_zlines.begin(),(it-1));
    /* Behavior at borders:
    * x < x_min -> index_x = -1
     * x = x_min -> index_x = 0
     * x = x_max -> index_x = m_xlines.size() - 1
     * x > x_max -> index_x = m_xlines.size() - 1
    */

    // check if the point is out of the mesh
    if(index_x < 0 || index_y < 0 || index_z < 0 ||
       index_x == (m_xlines.size()-1) || index_y == (m_ylines.size()-1) || index_z == (m_zlines.size()-1))
      return -1;
    else
      imap = index_x + (m_xlines.size() - 1) * index_y + (m_xlines.size() - 1) * (m_ylines.size() - 1) * index_z;

    if(debug && imap != -1) {
      if( x < nodes[elements[imap].emap[3]].x || x > nodes[elements[imap].emap[0]].x ||
          y < nodes[elements[imap].emap[3]].y || y > nodes[elements[imap].emap[2]].y ||
          z < nodes[elements[imap].emap[3]].z || z > nodes[elements[imap].emap[7]].z) {
        std::cout << "Element: " << imap << "\tPoint: (" << x << "," << y << "," << z << ")\n"
                << "x: " << nodes[elements[imap].emap[3]].x << " - " << nodes[elements[imap].emap[0]].x << "\n"
                << "y: " << nodes[elements[imap].emap[3]].y << " - " << nodes[elements[imap].emap[2]].y << "\n"
                << "z: " << nodes[elements[imap].emap[3]].z << " - " << nodes[elements[imap].emap[7]].z << "\n\n";
      }
    }
  }
  // final test - should never fail
  if (imap < 0 || imap > nElements) {
    std::cerr << className << "::FindElementCube:\n";
    std::cerr << "    Index of the element (imap is " << imap << ") is to large!"
              << "    Number of Elements: " << nElements
              << "\n    index_x: "   << index_x
              << "\n    index_y: " << index_y
              << "\n    index_z: " << index_z << std::endl;
    if (debug) {
      std::cout << className << "::FindElementCube:\n";
      std::cout << "    Point (" << x << "," << y << "," << z
                << ") not in the mesh, it is background or PEC.\n";
      std::cout << "    First node ("
                << nodes[elements[0].emap[3]].x << ","
                << nodes[elements[0].emap[3]].y << ","
                << nodes[elements[0].emap[3]].z
                << ") in the mesh.\n";
      std::cout << "    dx= " << (nodes[elements[0].emap[0]].x-nodes[elements[0].emap[3]].x)
                << ", dy= " << (nodes[elements[0].emap[2]].y-nodes[elements[0].emap[3]].y)
                << ", dz= " << (nodes[elements[0].emap[7]].z-nodes[elements[0].emap[3]].z)
                << "\n";
      std::cout << "    Last node (" << nodes[elements[nElements-1].emap[5]].x
                << "," << nodes[elements[nElements-1].emap[5]].y
                << "," << nodes[elements[nElements-1].emap[5]].z
                << ") in the mesh.\n";
      std::cout << "  dx= " << (nodes[elements[nElements-1].emap[0]].x-nodes[elements[nElements-1].emap[3]].x)
                << ", dy= " << (nodes[elements[nElements-1].emap[2]].y-nodes[elements[nElements-1].emap[3]].y)
                << ", dz= " << (nodes[elements[nElements-1].emap[7]].z-nodes[elements[nElements-1].emap[3]].z)
                << "\n";
    }
    return -1;
  }
  CoordinatesCube(x,y,z,t1,t2,t3,jac,det,imap);
  if (debug) {
    std::cout << className << "::FindElementCube:\n";
    std::cout << "Global: (" << x << "," << y << "," << z << ") in element "
              << imap << " (degenerate: "
              << elements[imap].degenerate << ")\n";
   std::cout << "      Node xyzV\n";
    for (int i = 0; i < 8; i++) {
      std::cout << "  " << elements[imap].emap[i]
                << " " << nodes[elements[imap].emap[i]].x
                << " " << nodes[elements[imap].emap[i]].y
                << " " << nodes[elements[imap].emap[i]].z
                << " " << nodes[elements[imap].emap[i]].v
                << "\n";
    }
  }
  return imap;
}
}