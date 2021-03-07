#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Garfield/ComponentTcadBase.hh"
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

template<size_t N>
void ComponentTcadBase<N>::WeightingField(const double x, const double y,
                                          const double z, double& wx, 
                                          double& wy, double& wz, 
                                          const std::string& /*label*/) {
  wx = wy = wz = 0.;
  if (m_wf.empty()) {
    std::cerr << m_className << "::WeightingField: Not available.\n";
    return;
  }
  Interpolate(x, y, z, m_wf, wx, wy, wz);
}

template<size_t N>
double ComponentTcadBase<N>::WeightingPotential(const double x, const double y,
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

template<size_t N>
bool ComponentTcadBase<N>::GetVoltageRange(double& vmin, double& vmax) {
  if (!m_ready) return false;
  vmin = m_pMin;
  vmax = m_pMax;
  return true;
}

template<size_t N>
bool ComponentTcadBase<N>::SetWeightingField(const std::string& datfile1,
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
  std::vector<std::array<double, N> > wf1;
  std::vector<double> wp1;
  if (!LoadWeightingField(datfile1, wf1, wp1)) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Could not import data from " << datfile1 << ".\n";
    return false;
  }

  // Then load the field/potential for the configuration with the potential 
  // at the electrode to be read out increased by small voltage dv. 
  std::vector<std::array<double, N> > wf2;
  std::vector<double> wy2;
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
    m_wf.resize(nVertices);
    for (size_t i = 0; i < nVertices; ++i) {
      for (size_t j = 0; j < N; ++j) {
        m_wf[i][j] = (wf2[i][j] - wf1[i][j]) * s;
      } 
    }
  }
  if (foundPotential) {
    m_wp.assign(nVertices, 0.);
    for (size_t i = 0; i < nVertices; ++i) {
      m_wp[i] = (wp2[i] - wp1[i]) * s; 
    }
  }
  return true;
}

template<size_t N>
bool ComponentTcadBase<N>::LoadWeightingField(
    const std::string& datafilename,
    std::vector<std::array<double, N> >& wf, std::vector<double>& wp) {

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
    // Strip white space from the beginning of the line.
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
    if (dataset != "ElectrostaticPotential" && dataset != "ElectricField") {
      continue;
    }
    bool field = false;
    if (dataset == "ElectricField") {
      wf.clear();
      wf.resize(nVertices);
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
    if (field) nValues /= N;
    // Mark the vertices belonging to this region.
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
      std::array<double, N> val;
      if (field) {
        for (size_t k = 0; k < N; ++k) datafile >> val[k];
      } else {
        datafile >> val[0];
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
        wf[ivertex] = val;
      } else {
        wp[ivertex] = val[0];
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

template<size_t N>
void ComponentTcadBase<N>::PrintRegions() const {

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

template<size_t N>
void ComponentTcadBase<N>::GetRegion(const size_t i, std::string& name,
                                     bool& active) const {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::GetRegion: Index out of range.\n";
    return;
  }
  name = m_regions[i].name;
  active = m_regions[i].drift;
}

template<size_t N>
void ComponentTcadBase<N>::SetDriftRegion(const size_t i) {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::SetDriftRegion: Index out of range.\n";
    return;
  }
  m_regions[i].drift = true;
}

template<size_t N>
void ComponentTcadBase<N>::UnsetDriftRegion(const size_t i) {
  if (i >= m_regions.size()) {
    std::cerr << m_className << "::UnsetDriftRegion: Index out of range.\n";
    return;
  }
  m_regions[i].drift = false;
}

template<size_t N>
void ComponentTcadBase<N>::SetMedium(const size_t i, Medium* medium) {
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

template<size_t N>
bool ComponentTcadBase<N>::SetDonor(const size_t donorNumber,
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

template<size_t N>
bool ComponentTcadBase<N>::SetAcceptor(const size_t acceptorNumber,
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

template<size_t N>
bool ComponentTcadBase<N>::ElectronAttachment(const double x, const double y,
                                              const double z, double& eta) {
  Interpolate(x, y, z, m_eAttachment, eta);
  return true;
}

template<size_t N>
bool ComponentTcadBase<N>::HoleAttachment(const double x, const double y,
                                          const double z, double& eta) {
  Interpolate(x, y, z, m_hAttachment, eta);
  return true;
}

template<size_t N>
bool ComponentTcadBase<N>::ElectronVelocity(const double x, const double y,
                                            const double z, double& vx, 
                                            double& vy, double& vz) {
  return Interpolate(x, y, z, m_eVelocity, vx, vy, vz);
}

template<size_t N>
bool ComponentTcadBase<N>::HoleVelocity(const double x, const double y,
                                        const double z, double& vx, 
                                        double& vy, double& vz) {
  return Interpolate(x, y, z, m_hVelocity, vx, vy, vz);
}

template<size_t N>
bool ComponentTcadBase<N>::GetElectronLifetime(const double x, const double y,
                                               const double z, double& tau) {
  return Interpolate(x, y, z, m_eLifetime, tau);
}

template<size_t N>
bool ComponentTcadBase<N>::GetHoleLifetime(const double x, const double y,
                                           const double z, double& tau) {
  return Interpolate(x, y, z, m_hLifetime, tau);
}

template<size_t N>
bool ComponentTcadBase<N>::GetElectronMobility(const double x, const double y,
                                               const double z, double& mob) {
  return Interpolate(x, y, z, m_eMobility, mob);
}

template<size_t N>
bool ComponentTcadBase<N>::GetHoleMobility(const double x, const double y,
                                           const double z, double& mob) {
  return Interpolate(x, y, z, m_hMobility, mob);
}

template<size_t N>
void ComponentTcadBase<N>::UpdatePeriodicity() {
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

template<size_t N>
void ComponentTcadBase<N>::Cleanup() {
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

template<size_t N>
void ComponentTcadBase<N>::MapCoordinates(std::array<double, N>& x, 
                                          std::array<bool, N>& mirr) const {
  mirr.fill(false);
  for (size_t i = 0; i < N; ++i) {
    // In case of periodicity, reduce to the cell volume.
    const double cellsx = m_bbMax[i] - m_bbMin[i];
    if (m_periodic[i]) {
      x[i] = m_bbMin[i] + fmod(x[i] - m_bbMin[i], cellsx);
      if (x[i] < m_bbMin[i]) x[i] += cellsx;
    } else if (m_mirrorPeriodic[i]) {
      double xNew = m_bbMin[i] + fmod(x[i] - m_bbMin[i], cellsx);
      if (xNew < m_bbMin[i]) xNew += cellsx;
      const int nx = int(floor(0.5 + (xNew - x[i]) / cellsx));
      if (nx != 2 * (nx / 2)) {
        xNew = m_bbMin[i] + m_bbMax[i] - xNew;
        mirr[i] = true;
      }
      x[i] = xNew;
    }
  }
}

template<size_t N>
size_t ComponentTcadBase<N>::FindRegion(const std::string& name) const {
  const auto nRegions = m_regions.size();
  for (size_t j = 0; j < nRegions; ++j) {
    if (name == m_regions[j].name) return j;
  }
  return m_regions.size();
}

template<size_t N>
void ComponentTcadBase<N>::UpdateAttachment() {

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

template class ComponentTcadBase<2>;
template class ComponentTcadBase<3>;
}
