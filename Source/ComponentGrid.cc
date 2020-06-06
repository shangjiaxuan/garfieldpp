#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>

#include "Garfield/ComponentGrid.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Utilities.hh"

namespace {

unsigned int GetFormat(std::string format) {
  std::transform(format.begin(), format.end(), format.begin(), toupper);
  unsigned int fmt = 0;
  if (format == "XY") {
    fmt = 1;
  } else if (format == "XYZ") {
    fmt = 2;
  } else if (format == "IJ") {
    fmt = 3;
  } else if (format == "IJK") {
    fmt = 4;
  } else if (format == "YXZ") {
    fmt = 5;
  }
  return fmt;
}

void PrintError(const std::string& fcn, const unsigned int line,
                const std::string& par) {
  std::cerr << fcn << ": Error reading line " << line << ".\n"
            << "    Could not read " << par << ".\n";
}

void PrintNotReady(const std::string& fcn) {
  std::cerr << fcn << ": Map not available.\n";
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

}  // namespace

namespace Garfield {

ComponentGrid::ComponentGrid() : ComponentBase() {
  m_className = "ComponentGrid";
}

void ComponentGrid::ElectricField(const double x, const double y,
                                  const double z, double& ex, double& ey,
                                  double& ez, double& p, Medium*& m,
                                  int& status) {
  m = nullptr;
  status = 0;

  // Make sure the field map has been loaded.
  if (!m_ready) {
    PrintNotReady(m_className + "::ElectricField");
    status = -10;
    return;
  }

  status = 0;
  bool active = true;
  if (!GetField(x, y, z, m_efields, ex, ey, ez, p, active)) {
    status = -11;
    return;
  }
  if (!active) {
    status = -5;
    return;
  }
  m = m_medium;
  if (!m) status = -5;
}

void ComponentGrid::ElectricField(const double x, const double y,
                                  const double z, double& ex, double& ey,
                                  double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

void ComponentGrid::WeightingField(const double x, const double y,
                                   const double z, double& wx, double& wy,
                                   double& wz, const std::string& /*label*/) {
  wx = wy = wz = 0.;
  if (!m_hasWfield) return;
  const double xx = x - m_wField_xOffset;
  const double yy = y - m_wField_yOffset;
  const double zz = z - m_wField_zOffset;
  double wp = 0.;
  bool active = true;
  GetField(xx, yy, zz, m_wfields, wx, wy, wz, wp, active);
}

double ComponentGrid::WeightingPotential(const double x, const double y,
                                         const double z,
                                         const std::string& /*label*/) {
  if (!m_hasWfield) return 0.;
  const double xx = x - m_wField_xOffset;
  const double yy = y - m_wField_yOffset;
  const double zz = z - m_wField_zOffset;
  double wx = 0., wy = 0., wz = 0.;
  double wp = 0.;
  bool active = true;
  if (!GetField(xx, yy, zz, m_wfields, wx, wy, wz, wp, active)) return 0.;
  return wp;
}

void ComponentGrid::DelayedWeightingField(const double x, const double y,
                                          const double z, const double t,
                                          double& wx, double& wy, double& wz,
                                          const std::string& /*label*/) {
  wx = wy = wz = 0.;
  if (m_wdtimes.empty()) return;
  // Assume no weighting field for times outside the range of available maps.
  if (t < m_wdtimes.front() || t > m_wdtimes.back()) return;

  const double xx = x - m_wField_xOffset;
  const double yy = y - m_wField_yOffset;
  const double zz = z - m_wField_zOffset;

  const auto it1 = std::upper_bound(m_wdtimes.cbegin(), m_wdtimes.cend(), t);
  const auto it0 = std::prev(it1);

  const double dt = t - *it0;
  double wp = 0.;
  const unsigned int i0 = it0 - m_wdtimes.cbegin();
  double wx0 = 0., wy0 = 0., wz0 = 0.;
  bool active = true;
  if (!GetField(xx, yy, zz, m_wdfields[i0], wx0, wy0, wz0, wp, active)) return;

  if (dt < Small || it1 == m_wdtimes.cend()) {
    wx = wx0;
    wy = wy0;
    wz = wz0;
    return;
  }
  const unsigned int i1 = it1 - m_wdtimes.cbegin();
  double wx1 = 0., wy1 = 0., wz1 = 0.;
  if (!GetField(xx, yy, zz, m_wdfields[i1], wx1, wy1, wz1, wp, active)) return;

  const double f1 = dt / (*it1 - *it0);
  const double f0 = 1. - f1;
  wx = f0 * wx0 + f1 * wx1;
  wy = f0 * wy0 + f1 * wy1;
  wz = f0 * wz0 + f1 * wz1;
}

void ComponentGrid::SetWeightingFieldOffset(const double x, const double y,
                                            const double z) {
  m_wField_xOffset = x;
  m_wField_yOffset = y;
  m_wField_zOffset = z;
}

void ComponentGrid::MagneticField(const double x, const double y,
                                  const double z, double& bx, double& by,
                                  double& bz, int& status) {
  status = 0;
  if (!m_hasBfield) {
    return ComponentBase::MagneticField(x, y, z, bx, by, bz, status);
  }

  double p = 0.;
  bool active = true;
  if (!GetField(x, y, z, m_bfields, bx, by, bz, p, active)) {
    status = -11;
  }
}

Medium* ComponentGrid::GetMedium(const double x, const double y,
                                 const double z) {
  // Make sure the field map has been loaded.
  if (!m_ready) {
    PrintNotReady(m_className + "::GetMedium");
    return nullptr;
  }

  if (!m_periodic[0] && !m_mirrorPeriodic[0] && (x < m_xMin || x > m_xMax)) {
    return nullptr;
  }
  if (!m_periodic[1] && !m_mirrorPeriodic[1] && (y < m_yMin || x > m_yMax)) {
    return nullptr;
  }
  if (!m_periodic[2] && !m_mirrorPeriodic[2] && (z < m_zMin || x > m_zMax)) {
    return nullptr;
  }
  if (m_active.empty()) return m_medium;

  bool mirrored = false;
  const double xx =
      Reduce(x, m_xMin, m_xMax, m_periodic[0], m_mirrorPeriodic[0], mirrored);
  const double yy =
      Reduce(y, m_yMin, m_yMax, m_periodic[1], m_mirrorPeriodic[1], mirrored);
  const double zz =
      Reduce(z, m_zMin, m_zMax, m_periodic[2], m_mirrorPeriodic[2], mirrored);
  // Get the indices.
  const double sx = (xx - m_xMin) / m_dx;
  const double sy = (yy - m_yMin) / m_dy;
  const double sz = (zz - m_zMin) / m_dz;
  const unsigned int i0 = static_cast<unsigned int>(std::floor(sx));
  const unsigned int j0 = static_cast<unsigned int>(std::floor(sy));
  const unsigned int k0 = static_cast<unsigned int>(std::floor(sz));
  const unsigned int i1 = std::min(i0 + 1, m_nX - 1);
  const unsigned int j1 = std::min(j0 + 1, m_nY - 1);
  const unsigned int k1 = std::min(k0 + 1, m_nZ - 1);
  if (m_active[i0][j0][k0] && m_active[i0][j0][k1] && m_active[i0][j1][k0] &&
      m_active[i0][j1][k1] && m_active[i1][j0][k0] && m_active[i1][j0][k1] &&
      m_active[i1][j1][k0] && m_active[i1][j1][k1]) {
    return m_medium;
  }
  return nullptr;
}

bool ComponentGrid::SetMesh(const unsigned int nx, const unsigned int ny,
                            const unsigned int nz, const double xmin,
                            const double xmax, const double ymin,
                            const double ymax, const double zmin,
                            const double zmax) {
  Reset();
  if (nx == 0 || ny == 0 || nz == 0) {
    std::cerr << m_className << "::SetMesh:\n"
              << "    Number of mesh elements must be positive.\n";
    return false;
  }
  if (xmin >= xmax) {
    std::cerr << m_className << "::SetMesh: Invalid x range.\n";
    return false;
  } else if (ymin >= ymax) {
    std::cerr << m_className << "::SetMesh: Invalid y range.\n";
    return false;
  } else if (zmin >= zmax) {
    std::cerr << m_className << "::SetMesh: Invalid z range.\n";
    return false;
  }
  m_nX = nx;
  m_nY = ny;
  m_nZ = nz;
  m_xMin = xmin;
  m_yMin = ymin;
  m_zMin = zmin;
  m_xMax = xmax;
  m_yMax = ymax;
  m_zMax = zmax;
  m_dx = m_nX > 1 ? (m_xMax - m_xMin) / (m_nX - 1) : (m_xMax - m_xMin);
  m_dy = m_nY > 1 ? (m_yMax - m_yMin) / (m_nY - 1) : (m_yMax - m_yMin);
  m_dz = m_nZ > 1 ? (m_zMax - m_zMin) / (m_nZ - 1) : (m_zMax - m_zMin);
  m_hasMesh = true;
  return true;
}

bool ComponentGrid::GetMesh(unsigned int& nx, unsigned int& ny,
                            unsigned int& nz, double& xmin, double& xmax,
                            double& ymin, double& ymax, double& zmin,
                            double& zmax) const {
  if (!m_hasMesh) return false;
  nx = m_nX;
  ny = m_nY;
  nz = m_nZ;
  xmin = m_xMin;
  ymin = m_yMin;
  zmin = m_zMin;
  xmax = m_xMax;
  ymax = m_yMax;
  zmax = m_zMax;
  return true;
}

bool ComponentGrid::LoadElectricField(const std::string& fname,
                                      const std::string& fmt, const bool withP,
                                      const bool withFlag, const double scaleX,
                                      const double scaleE,
                                      const double scaleP) {
  m_ready = false;
  m_hasPotential = m_hasEfield = false;
  m_active.assign(m_nX, std::vector<std::vector<bool> >(
                            m_nY, std::vector<bool>(m_nZ, true)));
  // Read the file.
  m_pMin = withP ? +1. : 0.;
  m_pMax = withP ? -1. : 0.;
  if (!LoadField(fname, fmt, withP, withFlag, scaleX, scaleE, scaleP,
                 m_efields)) {
    return false;
  }
  m_hasEfield = true;
  m_ready = true;
  if (withP) m_hasPotential = true;
  return true;
}

bool ComponentGrid::LoadWeightingField(const std::string& fname,
                                       const std::string& fmt, const bool withP,
                                       const double scaleX, const double scaleE,
                                       const double scaleP) {
  m_hasWfield = false;
  // Read the file.
  if (!LoadField(fname, fmt, withP, false, scaleX, scaleE, scaleP, m_wfields)) {
    return false;
  }
  m_hasWfield = true;
  return true;
}

bool ComponentGrid::LoadWeightingField(const std::string& fname,
                                       const std::string& fmt, const double t,
                                       const bool withP, const double scaleX,
                                       const double scaleE,
                                       const double scaleP) {
  std::vector<std::vector<std::vector<Node> > > wfield;
  // Read the file.
  if (!LoadField(fname, fmt, withP, false, scaleX, scaleE, scaleP, wfield)) {
    return false;
  }
  if (m_wdtimes.empty() || t > m_wdtimes.back()) {
    m_wdtimes.push_back(t);
    m_wdfields.push_back(std::move(wfield));
  } else {
    const auto it = std::upper_bound(m_wdtimes.cbegin(), m_wdtimes.cend(), t);
    const auto n = std::distance(m_wdtimes.cbegin(), it);
    m_wdtimes.insert(it, t);
    m_wdfields.insert(m_wdfields.cbegin() + n, std::move(wfield));
  }
  return true;
}

bool ComponentGrid::LoadMagneticField(const std::string& fname,
                                      const std::string& fmt,
                                      const double scaleX,
                                      const double scaleB) {
  m_hasBfield = false;
  // Read the file.
  if (!LoadField(fname, fmt, false, false, scaleX, scaleB, 1., m_bfields)) {
    return false;
  }
  m_hasBfield = true;
  return true;
}

bool ComponentGrid::SaveElectricField(ComponentBase* cmp,
                                      const std::string& filename,
                                      const std::string& format) {
  if (!cmp) {
    std::cerr << m_className << "::SaveElectricField: Null pointer.\n";
    return false;
  }
  if (!m_hasMesh) {
    std::cerr << m_className << "::SaveElectricField: Mesh not set.\n";
    return false;
  }
  const unsigned int fmt = GetFormat(format);
  if (fmt == 0) {
    std::cerr << m_className << "::SaveElectricField:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out);
  if (!outfile) {
    std::cerr << m_className << "::SaveElectricField:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }
  std::cout << m_className << "::SaveElectricField:\n"
            << "    Exporting field/potential to " << filename << ".\n"
            << "    Be patient...\n";
  PrintProgress(0.);
  outfile << "# XMIN = " << m_xMin << ", XMAX = " << m_xMax << ", NX = " << m_nX
          << "\n";
  outfile << "# YMIN = " << m_yMin << ", YMAX = " << m_yMax << ", NY = " << m_nY
          << "\n";
  outfile << "# ZMIN = " << m_zMin << ", ZMAX = " << m_zMax << ", NZ = " << m_nZ
          << "\n";

  const unsigned int nValues = m_nX * m_nY * m_nZ;
  const unsigned int nPrint =
      std::pow(10, static_cast<unsigned int>(
                       std::max(std::floor(std::log10(nValues)) - 1, 1.)));
  unsigned int nLines = 0;
  Medium* medium = nullptr;
  int status = 0;
  for (unsigned int i = 0; i < m_nX; ++i) {
    const double x = m_xMin + i * m_dx;
    for (unsigned int j = 0; j < m_nY; ++j) {
      const double y = m_yMin + j * m_dy;
      for (unsigned int k = 0; k < m_nZ; ++k) {
        const double z = m_zMin + k * m_dz;
        if (fmt == 1) {
          outfile << x << "  " << y << "  ";
        } else if (fmt == 2) {
          outfile << x << "  " << y << "  " << z << "  ";
        } else if (fmt == 3) {
          outfile << i << "  " << j << "  ";
        } else if (fmt == 4) {
          outfile << i << "  " << j << "  " << k << "  ";
        } else if (fmt == 5) {
          outfile << y << "  " << x << "  " << z << "  ";
        }
        double ex = 0., ey = 0., ez = 0., v = 0.;
        cmp->ElectricField(x, y, z, ex, ey, ez, v, medium, status);
        outfile << ex << "  " << ey << "  " << ez << "  " << v << "\n";
        ++nLines;
        if (nLines % nPrint == 0) PrintProgress(double(nLines) / nValues);
      }
    }
  }
  outfile.close();
  std::cout << std::endl << m_className << "::SaveElectricField: Done.\n";
  return true;
}

bool ComponentGrid::SaveWeightingField(ComponentBase* cmp,
                                       const std::string& id,
                                       const std::string& filename,
                                       const std::string& format) {
  if (!cmp) {
    std::cerr << m_className << "::SaveWeightingField: Null pointer.\n";
    return false;
  }
  if (!m_hasMesh) {
    std::cerr << m_className << "::SaveWeightingField: Mesh not set.\n";
    return false;
  }
  const unsigned int fmt = GetFormat(format);
  if (fmt == 0) {
    std::cerr << m_className << "::SaveWeightingField:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out);
  if (!outfile) {
    std::cerr << m_className << "::SaveWeightingField:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }
  std::cout << m_className << "::SaveWeightingField:\n"
            << "    Exporting field/potential to " << filename << ".\n"
            << "    Be patient...\n";
  PrintProgress(0.);
  outfile << "# XMIN = " << m_xMin << ", XMAX = " << m_xMax << ", NX = " << m_nX
          << "\n";
  outfile << "# YMIN = " << m_yMin << ", YMAX = " << m_yMax << ", NY = " << m_nY
          << "\n";
  outfile << "# ZMIN = " << m_zMin << ", ZMAX = " << m_zMax << ", NZ = " << m_nZ
          << "\n";
  const unsigned int nValues = m_nX * m_nY * m_nZ;
  const unsigned int nPrint =
      std::pow(10, static_cast<unsigned int>(
                       std::max(std::floor(std::log10(nValues)) - 1, 1.)));
  unsigned int nLines = 0;
  for (unsigned int i = 0; i < m_nX; ++i) {
    const double x = m_xMin + i * m_dx;
    for (unsigned int j = 0; j < m_nY; ++j) {
      const double y = m_yMin + j * m_dy;
      for (unsigned int k = 0; k < m_nZ; ++k) {
        const double z = m_zMin + k * m_dz;
        if (fmt == 1) {
          outfile << x << "  " << y << "  ";
        } else if (fmt == 2) {
          outfile << x << "  " << y << "  " << z << "  ";
        } else if (fmt == 3) {
          outfile << i << "  " << j << "  ";
        } else if (fmt == 4) {
          outfile << i << "  " << j << "  " << k << "  ";
        } else if (fmt == 5) {
          outfile << y << "  " << x << "  " << z << "  ";
        }
        double wx = 0., wy = 0., wz = 0.;
        cmp->WeightingField(x, y, z, wx, wy, wz, id);
        const double v = cmp->WeightingPotential(x, y, z, id);
        outfile << wx << "  " << wy << "  " << wz << "  " << v << "\n";
        ++nLines;
        if (nLines % nPrint == 0) PrintProgress(double(nLines) / nValues);
      }
    }
  }
  outfile.close();
  std::cout << std::endl << m_className << "::SaveWeightingField: Done.\n";
  return true;
}

bool ComponentGrid::LoadMesh(const std::string& filename, std::string format,
                             const double scaleX) {
  const unsigned int fmt = GetFormat(format);
  if (fmt == 0) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }

  // Keep track of which mesh parameters we have found.
  std::bitset<9> found;
  found.reset();
  double xmin = 0., ymin = 0., zmin = 0.;
  double xmax = 0., ymax = 0., zmax = 0.;
  unsigned int nx = 0, ny = 0, nz = 0;

  // Parse the comment lines in the file.
  std::ifstream infile;
  infile.open(filename.c_str(), std::ios::in);
  if (!infile) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }
  std::string line;
  unsigned int nLines = 0;
  while (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;
    // Strip white space from the beginning of the line.
    ltrim(line);
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip lines that are not comments.
    if (line[0] != '#' && !(line[0] == '/' && line[1] == '/')) {
      continue;
    }
    std::size_t pos0 = 0;
    std::size_t pos1 = line.find("=", pos0);
    while (pos1 != std::string::npos) {
      std::string key = line.substr(pos0, pos1 - pos0);
      std::transform(key.begin(), key.end(), key.begin(), toupper);
      const std::size_t pos2 = line.find_first_of(",;", pos1 + 1);
      std::istringstream val(line.substr(pos1 + 1, pos2 - pos1 - 1));
      if (key.find("XMIN") != std::string::npos) {
        val >> xmin;
        found.set(0);
      } else if (key.find("YMIN") != std::string::npos) {
        val >> ymin;
        found.set(1);
      } else if (key.find("ZMIN") != std::string::npos) {
        val >> zmin;
        found.set(2);
      } else if (key.find("XMAX") != std::string::npos) {
        val >> xmax;
        found.set(3);
      } else if (key.find("YMAX") != std::string::npos) {
        val >> ymax;
        found.set(4);
      } else if (key.find("ZMAX") != std::string::npos) {
        val >> zmax;
        found.set(5);
      } else if (key.find("NX") != std::string::npos) {
        val >> nx;
        found.set(6);
      } else if (key.find("NY") != std::string::npos) {
        val >> ny;
        found.set(7);
      } else if (key.find("NZ") != std::string::npos) {
        val >> nz;
        found.set(8);
      }
      if (pos2 == std::string::npos) break;
      pos0 = pos2 + 1;
      pos1 = line.find("=", pos0);
    }
  }
  infile.close();

  if (fmt == 1 || fmt == 3) {
    // Try to complement missing information on the z-range.
    if (!found[8]) {
      nz = 1;
      found.set(8);
    }
    if (!found[2]) {
      if (found[0] || found[1] || found[3] || found[4] || found[5]) {
        zmin = -std::max(
            {fabs(xmin), fabs(xmax), fabs(ymin), fabs(ymax), fabs(zmax)});
      } else {
        zmin = -100.;
      }
      found.set(2);
    }
    if (!found[5]) {
      zmax = std::max(
          {fabs(xmin), fabs(xmax), fabs(ymin), fabs(ymax), fabs(zmin)});
      found.set(5);
    }
  }
  if (found.all()) {
    std::cout << m_className << "::LoadMesh:\n";
    std::printf("%12.6f < x [cm] < %12.6f, %5u points\n", xmin, xmax, nx);
    std::printf("%12.6f < y [cm] < %12.6f, %5u points\n", ymin, ymax, ny);
    std::printf("%12.6f < z [cm] < %12.6f, %5u points\n", zmin, zmax, nz);
    return SetMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax);
  }

  if ((fmt == 3 || fmt == 4) && !(found[0] && found[3])) {
    std::cerr << m_className << "::LoadMesh: x-limits not found.\n";
    return false;
  } else if ((fmt == 3 || fmt == 4) && !(found[1] && found[4])) {
    std::cerr << m_className << "::LoadMesh: y-limits not found.\n";
    return false;
  } else if (fmt == 4 && !(found[2] && found[5])) {
    std::cerr << m_className << "::LoadMesh: z-limits not found.\n";
    return false;
  }

  unsigned int nValues = 0;
  infile.open(filename.c_str(), std::ios::in);
  if (!infile) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }

  if (!found[0]) xmin = std::numeric_limits<double>::max();
  if (!found[1]) ymin = std::numeric_limits<double>::max();
  if (!found[2]) zmin = std::numeric_limits<double>::max();
  if (!found[3]) xmax = std::numeric_limits<double>::min();
  if (!found[4]) ymax = std::numeric_limits<double>::min();
  if (!found[5]) zmax = std::numeric_limits<double>::min();
  constexpr double tol = 1.e-10;
  auto cmp = [](double x, double y) {
    return x < y - tol * (std::fabs(x) + std::fabs(y));
  };
  std::set<double, decltype(cmp)> xLines(cmp);
  std::set<double, decltype(cmp)> yLines(cmp);
  std::set<double, decltype(cmp)> zLines(cmp);
  nLines = 0;
  bool bad = false;
  while (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;
    // Strip white space from the beginning of the line.
    ltrim(line);
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip comments.
    if (line[0] == '#') continue;
    if (line[0] == '/' && line[1] == '/') continue;
    std::istringstream data;
    data.str(line);
    if (fmt == 1) {
      // "XY"
      double x, y;
      data >> x >> y;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      if (!found[0]) xmin = std::min(x, xmin);
      if (!found[1]) ymin = std::min(y, ymin);
      if (!found[3]) xmax = std::max(x, xmin);
      if (!found[4]) ymax = std::max(y, ymin);
      xLines.insert(x);
      yLines.insert(y);
    } else if (fmt == 2) {
      // "XYZ"
      double x, y, z;
      data >> x >> y >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      if (!found[0]) xmin = std::min(x, xmin);
      if (!found[1]) ymin = std::min(y, ymin);
      if (!found[2]) zmin = std::min(z, zmin);
      if (!found[3]) xmax = std::max(x, xmax);
      if (!found[4]) ymax = std::max(y, ymax);
      if (!found[5]) zmax = std::max(z, zmax);
      xLines.insert(x);
      yLines.insert(y);
      zLines.insert(z);
    } else if (fmt == 3) {
      // "IJ"
      unsigned int i = 0, j = 0;
      data >> i >> j;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "indices");
        bad = true;
        break;
      }
      if (!found[6]) nx = std::max(nx, i);
      if (!found[7]) ny = std::max(ny, j);
    } else if (fmt == 4) {
      // "IJK"
      unsigned int i = 0, j = 0, k = 0;
      data >> i >> j >> k;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "indices");
        bad = true;
        break;
      }
      if (!found[6]) nx = std::max(nx, i);
      if (!found[7]) ny = std::max(ny, j);
      if (!found[8]) nz = std::max(nz, k);
    } else if (fmt == 5) {
      // "YXZ"
      double x, y, z;
      data >> y >> x >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      if (!found[0]) xmin = std::min(x, xmin);
      if (!found[1]) ymin = std::min(y, ymin);
      if (!found[2]) zmin = std::min(z, zmin);
      if (!found[3]) xmax = std::max(x, xmax);
      if (!found[4]) ymax = std::max(y, ymax);
      if (!found[5]) zmax = std::max(z, zmax);
      xLines.insert(x);
      yLines.insert(y);
      zLines.insert(z);
    }
    ++nValues;
  }
  infile.close();
  if (bad) return false;

  if (fmt == 1 || fmt == 2 || fmt == 5) {
    if (!found[6]) nx = xLines.size();
    if (!found[7]) ny = yLines.size();
    if (!found[8]) nz = zLines.size();
  }

  std::cout << m_className << "::LoadMesh:\n";
  std::printf("%12.6f < x [cm] < %12.6f, %5u points\n", xmin, xmax, nx);
  std::printf("%12.6f < y [cm] < %12.6f, %5u points\n", ymin, ymax, ny);
  std::printf("%12.6f < z [cm] < %12.6f, %5u points\n", zmin, zmax, nz);
  unsigned int nExpected = nx * ny;
  if (fmt == 2 || fmt == 4 || fmt == 5) nExpected *= nz;
  if (nExpected != nValues) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "   Warning: Expected " << nExpected << " lines, read "
              << nValues << ".\n";
  }
  return SetMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax);
}

bool ComponentGrid::LoadField(
    const std::string& filename, std::string format, const bool withPotential,
    const bool withFlag, const double scaleX, const double scaleF,
    const double scaleP,
    std::vector<std::vector<std::vector<Node> > >& fields) {
  if (!m_hasMesh) {
    if (!LoadMesh(filename, format, scaleX)) {
      std::cerr << m_className << "::LoadField: Mesh not set.\n";
      return false;
    }
  }

  const unsigned int fmt = GetFormat(format);
  if (fmt == 0) {
    std::cerr << m_className << "::LoadField:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }

  // Set up the grid.
  Initialise(fields);

  unsigned int nValues = 0;
  // Keep track of which elements have been read.
  std::vector<std::vector<std::vector<bool> > > isSet(
      m_nX,
      std::vector<std::vector<bool> >(m_nY, std::vector<bool>(m_nZ, false)));

  std::ifstream infile;
  infile.open(filename.c_str(), std::ios::in);
  if (!infile) {
    std::cerr << m_className << "::LoadField:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }

  std::string line;
  unsigned int nLines = 0;
  bool bad = false;
  while (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;
    // Strip white space from beginning of line.
    ltrim(line);
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip comments.
    if (line[0] == '#') continue;
    if (line[0] == '/' && line[1] == '/') continue;
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;
    double fx = 0.;
    double fy = 0.;
    double fz = 0.;
    double p = 0.;
    std::istringstream data;
    data.str(line);
    if (fmt == 1) {
      // "XY"
      double x, y;
      data >> x >> y;
      if (data.fail()) {
        PrintError(m_className + "::LoadField", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      const double u = std::round((x - m_xMin) / m_dx);
      const double v = std::round((y - m_yMin) / m_dy);
      i = u < 0. ? 0 : static_cast<unsigned int>(u);
      j = v < 0. ? 0 : static_cast<unsigned int>(v);
      if (i >= m_nX) i = m_nX - 1;
      if (j >= m_nY) j = m_nY - 1;
    } else if (fmt == 2) {
      // "XYZ"
      double x, y, z;
      data >> x >> y >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadField", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      const double u = std::round((x - m_xMin) / m_dx);
      const double v = std::round((y - m_yMin) / m_dy);
      const double w = std::round((z - m_zMin) / m_dz);
      i = u < 0. ? 0 : static_cast<unsigned int>(u);
      j = v < 0. ? 0 : static_cast<unsigned int>(v);
      j = w < 0. ? 0 : static_cast<unsigned int>(w);
      if (i >= m_nX) i = m_nX - 1;
      if (j >= m_nY) j = m_nY - 1;
      if (k >= m_nZ) k = m_nZ - 1;
    } else if (fmt == 3) {
      // "IJ"
      data >> i >> j;
      if (data.fail()) {
        PrintError(m_className + "::LoadField", nLines, "indices");
        bad = true;
        break;
      }
    } else if (fmt == 4) {
      // "IJK"
      data >> i >> j >> k;
      if (data.fail()) {
        PrintError(m_className + "::LoadField", nLines, "indices");
        bad = true;
        break;
      }
    } else if (fmt == 5) {
      // "YXZ"
      double x, y, z;
      data >> y >> x >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadField", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      const double u = std::round((x - m_xMin) / m_dx);
      const double v = std::round((y - m_yMin) / m_dy);
      const double w = std::round((z - m_zMin) / m_dz);
      i = u < 0. ? 0 : static_cast<unsigned int>(u);
      j = v < 0. ? 0 : static_cast<unsigned int>(v);
      j = w < 0. ? 0 : static_cast<unsigned int>(w);
      if (i >= m_nX) i = m_nX - 1;
      if (j >= m_nY) j = m_nY - 1;
      if (k >= m_nZ) k = m_nZ - 1;
    }
    // Check the indices.
    if (i >= m_nX || j >= m_nY || k >= m_nZ) {
      std::cerr << m_className << "::LoadField:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Index (" << i << ", " << j << ", " << k
                << ") out of range.\n";
      continue;
    }
    if (isSet[i][j][k]) {
      std::cerr << m_className << "::LoadField:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Node (" << i << ", " << j << ", " << k
                << ") has already been set.\n";
      continue;
    }
    // Get the field values.
    if (fmt == 1 || fmt == 3) {
      // Two-dimensional field-map
      fz = 0.;
      data >> fx >> fy;
    } else if (fmt == 5) {
      data >> fy >> fx >> fz;
    } else {
      data >> fx >> fy >> fz;
    }
    if (data.fail()) {
      PrintError(m_className + "::LoadField", nLines, "field components");
      bad = true;
      break;
    }
    fx *= scaleF;
    fy *= scaleF;
    fz *= scaleF;
    if (withPotential) {
      data >> p;
      if (data.fail()) {
        PrintError(m_className + "::LoadField", nLines, "potential");
        bad = true;
        break;
      }
      p *= scaleP;
      if (m_pMin > m_pMax) {
        // First value.
        m_pMin = p;
        m_pMax = p;
      } else {
        if (p < m_pMin) m_pMin = p;
        if (p > m_pMax) m_pMax = p;
      }
    }
    int flag = 0;
    if (withFlag) {
      data >> flag;
      if (data.fail()) {
        PrintError(m_className + "::LoadField", nLines, "region");
        bad = true;
        break;
      }
    }
    const bool isActive = flag == 0 ? false : true;
    if (fmt == 1 || fmt == 3) {
      // Two-dimensional field-map
      for (unsigned int kk = 0; kk < m_nZ; ++kk) {
        fields[i][j][kk].fx = fx;
        fields[i][j][kk].fy = fy;
        fields[i][j][kk].fz = fz;
        fields[i][j][kk].v = p;
        if (withFlag) m_active[i][j][kk] = isActive;
        isSet[i][j][kk] = true;
      }
    } else {
      fields[i][j][k].fx = fx;
      fields[i][j][k].fy = fy;
      fields[i][j][k].fz = fz;
      fields[i][j][k].v = p;
      isSet[i][j][k] = true;
    }
    ++nValues;
  }
  infile.close();
  if (bad) return false;
  std::cout << m_className << "::LoadField:\n"
            << "    Read " << nValues << " values from " << filename << ".\n";
  unsigned int nExpected = m_nX * m_nY;
  if (fmt == 2 || fmt == 4 || fmt == 5) nExpected *= m_nZ;
  if (nExpected != nValues) {
    std::cerr << m_className << "::LoadField:\n"
              << "   Expected " << nExpected << " values.\n";
  }
  return true;
}

bool ComponentGrid::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                   double& xmax, double& ymax, double& zmax) {
  if (!m_ready) return false;
  if (m_periodic[0] || m_mirrorPeriodic[0]) {
    xmin = -INFINITY;
    xmax = +INFINITY;
  } else {
    xmin = m_xMin;
    xmax = m_xMax;
  }

  if (m_periodic[1] || m_mirrorPeriodic[1]) {
    ymin = -INFINITY;
    ymax = +INFINITY;
  } else {
    ymin = m_yMin;
    ymax = m_yMax;
  }

  if (m_periodic[2] || m_mirrorPeriodic[2]) {
    zmin = -INFINITY;
    zmax = +INFINITY;
  } else {
    zmin = m_zMin;
    zmax = m_zMax;
  }
  return true;
}

bool ComponentGrid::GetVoltageRange(double& vmin, double& vmax) {
  if (!m_ready) return false;
  vmin = m_pMin;
  vmax = m_pMax;
  return true;
}

bool ComponentGrid::GetElectricFieldRange(double& exmin, double& exmax,
                                          double& eymin, double& eymax,
                                          double& ezmin, double& ezmax) {
  if (!m_ready) {
    PrintNotReady(m_className + "::GetElectricFieldRange");
    return false;
  }

  exmin = exmax = m_efields[0][0][0].fx;
  eymin = eymax = m_efields[0][0][0].fy;
  ezmin = ezmax = m_efields[0][0][0].fz;
  for (unsigned int i = 0; i < m_nX; ++i) {
    for (unsigned int j = 0; j < m_nY; ++j) {
      for (unsigned int k = 0; k < m_nZ; ++k) {
        const Node& node = m_efields[i][j][k];
        if (node.fx < exmin) exmin = node.fx;
        if (node.fx > exmax) exmax = node.fx;
        if (node.fy < eymin) eymin = node.fy;
        if (node.fy > eymax) eymax = node.fy;
        if (node.fz < ezmin) ezmin = node.fz;
        if (node.fz > ezmax) ezmax = node.fz;
      }
    }
  }
  return true;
}

void ComponentGrid::SetMedium(Medium* m) {
  if (!m) {
    std::cerr << m_className << "::SetMedium: Null pointer.\n";
  }
  m_medium = m;
}

bool ComponentGrid::GetField(
    const double xi, const double yi, const double zi,
    const std::vector<std::vector<std::vector<Node> > >& field, double& fx,
    double& fy, double& fz, double& p, bool& active) {
  if (!m_hasMesh) {
    std::cerr << m_className << "::GetField: Mesh is not set.\n";
    return false;
  }

  // Reduce the point to the basic cell (in case of periodicity) and
  // check if it is inside the mesh.
  bool xMirrored = false;
  const double x =
      Reduce(xi, m_xMin, m_xMax, m_periodic[0], m_mirrorPeriodic[0], xMirrored);
  if (x < m_xMin || x > m_xMax) return false;
  bool yMirrored = false;
  const double y =
      Reduce(yi, m_yMin, m_yMax, m_periodic[1], m_mirrorPeriodic[1], yMirrored);
  if (y < m_yMin || y > m_yMax) return false;
  bool zMirrored = false;
  const double z =
      Reduce(zi, m_zMin, m_zMax, m_periodic[2], m_mirrorPeriodic[2], zMirrored);
  if (z < m_zMin || z > m_zMax) return false;

  // Get the indices.
  const double sx = (x - m_xMin) / m_dx;
  const double sy = (y - m_yMin) / m_dy;
  const double sz = (z - m_zMin) / m_dz;
  const unsigned int i0 = static_cast<unsigned int>(std::floor(sx));
  const unsigned int j0 = static_cast<unsigned int>(std::floor(sy));
  const unsigned int k0 = static_cast<unsigned int>(std::floor(sz));
  const double ux = sx - i0;
  const double uy = sy - j0;
  const double uz = sz - k0;
  const unsigned int i1 = std::min(i0 + 1, m_nX - 1);
  const unsigned int j1 = std::min(j0 + 1, m_nY - 1);
  const unsigned int k1 = std::min(k0 + 1, m_nZ - 1);
  const double vx = 1. - ux;
  const double vy = 1. - uy;
  const double vz = 1. - uz;
  if (!m_active.empty()) {
    active = m_active[i0][j0][k0] && m_active[i0][j0][k1] &&
             m_active[i0][j1][k0] && m_active[i0][j1][k1] &&
             m_active[i1][j0][k0] && m_active[i1][j0][k1] &&
             m_active[i1][j1][k0] && m_active[i1][j1][k1];
  }
  const Node& n000 = field[i0][j0][k0];
  const Node& n100 = field[i1][j0][k0];
  const Node& n010 = field[i0][j1][k0];
  const Node& n110 = field[i1][j1][k0];
  const Node& n001 = field[i0][j0][k1];
  const Node& n101 = field[i1][j0][k1];
  const Node& n011 = field[i0][j1][k1];
  const Node& n111 = field[i1][j1][k1];

  if (m_debug) {
    std::cout << m_className << "::GetField: Determining field at (" << xi
              << ", " << yi << ", " << zi << ").\n"
              << "    X: " << i0 << " (" << ux << ") - " << i1 << " (" << vx
              << ").\n"
              << "    Y: " << j0 << " (" << uy << ") - " << j1 << " (" << vy
              << ").\n"
              << "    Z: " << k0 << " (" << uz << ") - " << k1 << " (" << vz
              << ").\n";
  }
  fx = ((n000.fx * vx + n100.fx * ux) * vy +
        (n010.fx * vx + n110.fx * ux) * uy) *
           vz +
       ((n001.fx * vx + n101.fx * ux) * vy +
        (n011.fx * vx + n111.fx * ux) * uy) *
           uz;
  fy = ((n000.fy * vx + n100.fy * ux) * vy +
        (n010.fy * vx + n110.fy * ux) * uy) *
           vz +
       ((n001.fy * vx + n101.fy * ux) * vy +
        (n011.fy * vx + n111.fy * ux) * uy) *
           uz;
  fz = ((n000.fz * vx + n100.fz * ux) * vy +
        (n010.fz * vx + n110.fz * ux) * uy) *
           vz +
       ((n001.fz * vx + n101.fz * ux) * vy +
        (n011.fz * vx + n111.fz * ux) * uy) *
           uz;
  p = ((n000.v * vx + n100.v * ux) * vy + (n010.v * vx + n110.v * ux) * uy) *
          vz +
      ((n001.v * vx + n101.v * ux) * vy + (n011.v * vx + n111.v * ux) * uy) *
          uz;
  if (xMirrored) fx = -fx;
  if (yMirrored) fy = -fy;
  if (zMirrored) fz = -fz;
  return true;
}

bool ComponentGrid::GetElectricField(const unsigned int i, const unsigned int j,
                                     const unsigned int k, double& v,
                                     double& ex, double& ey, double& ez) const {
  v = ex = ey = ez = 0.;
  if (!m_ready) {
    if (!m_hasMesh) {
      std::cerr << m_className << "::GetElectricField: Mesh not set.\n";
      return false;
    }
    PrintNotReady(m_className + "::GetElectricField");
    return false;
  }
  if (i >= m_nX || j >= m_nY || k >= m_nZ) {
    std::cerr << m_className << "::GetElectricField: Index out of range.\n";
    return false;
  }
  const Node& node = m_efields[i][j][k];
  v = node.v;
  ex = node.fx;
  ey = node.fy;
  ez = node.fz;
  return true;
}

void ComponentGrid::Reset() {
  m_efields.clear();
  m_bfields.clear();
  m_wfields.clear();
  m_eattachment.clear();
  m_hattachment.clear();

  m_wdfields.clear();
  m_wdtimes.clear();

  m_active.clear();

  m_nX = m_nY = m_nZ = 0;
  m_xMin = m_yMin = m_zMin = 0.;
  m_xMax = m_yMax = m_zMax = 0.;
  m_pMin = m_pMax = 0.;
  m_medium = nullptr;

  m_hasMesh = false;
  m_hasPotential = false;
  m_hasEfield = false;
  m_hasBfield = false;
  m_hasWfield = false;
  m_ready = false;

  m_wField_xOffset = 0.;
  m_wField_yOffset = 0.;
  m_wField_zOffset = 0.;
}

void ComponentGrid::UpdatePeriodicity() {
  if (!m_ready) {
    PrintNotReady(m_className + "::UpdatePeriodicity");
    return;
  }

  // Check for conflicts.
  for (unsigned int i = 0; i < 3; ++i) {
    if (m_periodic[i] && m_mirrorPeriodic[i]) {
      std::cerr << m_className << "::UpdatePeriodicity:\n"
                << "    Both simple and mirror periodicity requested. Reset.\n";
      m_periodic[i] = m_mirrorPeriodic[i] = false;
    }
  }

  if (m_axiallyPeriodic[0] || m_axiallyPeriodic[1] || m_axiallyPeriodic[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Axial symmetry is not supported. Reset.\n";
    m_axiallyPeriodic.fill(false);
  }

  if (m_rotationSymmetric[0] || m_rotationSymmetric[1] ||
      m_rotationSymmetric[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Rotation symmetry is not supported. Reset.\n";
    m_rotationSymmetric.fill(false);
  }
}

double ComponentGrid::Reduce(const double xin, const double xmin,
                             const double xmax, const bool simplePeriodic,
                             const bool mirrorPeriodic, bool& mirrored) const {
  // In case of periodicity, reduce the coordinate to the basic cell.
  double x = xin;
  const double lx = xmax - xmin;
  if (simplePeriodic) {
    x = xmin + fmod(x - xmin, lx);
    if (x < xmin) x += lx;
  } else if (mirrorPeriodic) {
    double xNew = xmin + fmod(x - xmin, lx);
    if (xNew < xmin) xNew += lx;
    const int nx = int(floor(0.5 + (xNew - x) / lx));
    if (nx != 2 * (nx / 2)) {
      xNew = xmin + xmax - xNew;
      mirrored = true;
    }
    x = xNew;
  }
  return x;
}

void ComponentGrid::Initialise(
    std::vector<std::vector<std::vector<Node> > >& fields) {
  fields.resize(m_nX);
  for (unsigned int i = 0; i < m_nX; ++i) {
    fields[i].resize(m_nY);
    for (unsigned int j = 0; j < m_nY; ++j) {
      fields[i][j].resize(m_nZ);
      for (unsigned int k = 0; k < m_nZ; ++k) {
        fields[i][j][k].fx = 0.;
        fields[i][j][k].fy = 0.;
        fields[i][j][k].fz = 0.;
        fields[i][j][k].v = 0.;
      }
    }
  }
}

bool ComponentGrid::LoadElectronAttachment(const std::string& fname,
                                           const std::string& fmt, 
                                           const unsigned int col,
                                           const double scaleX) {
  // Read the file.
  return LoadData(fname, fmt, scaleX, m_eattachment, col);
}

bool ComponentGrid::LoadHoleAttachment(const std::string& fname,
                                       const std::string& fmt, 
                                       const unsigned int col,
                                       const double scaleX) {
  // Read the file.
  return LoadData(fname, fmt, scaleX, m_hattachment, col);
}

bool ComponentGrid::LoadData(
    const std::string& filename, std::string format, const double scaleX,
    std::vector<std::vector<std::vector<double> > >& tab, 
    const unsigned int col) {
  if (!m_hasMesh) {
    if (!LoadMesh(filename, format, scaleX)) {
      std::cerr << m_className << "::LoadData: Mesh not set.\n";
      return false;
    }
  }

  const unsigned int fmt = GetFormat(format);
  if (fmt == 0) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }
  // Check the column index.
  unsigned int offset = 0;
  if (fmt == 1 || fmt == 3) {
    if (col < 2) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Unexpected column index (" << col << ").\n";
      return false; 
    }
    offset = 2;
  } else {
    if (col < 3) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Unexpected column index (" << col << ").\n";
      return false; 
    }
    offset = 3;
  } 

  // Set up the grid.
  tab.assign(
      m_nX, 
      std::vector<std::vector<double> >(m_nY, std::vector<double>(m_nZ, 0.)));

  unsigned int nValues = 0;
  // Keep track of which elements have been read.
  std::vector<std::vector<std::vector<bool> > > isSet(
      m_nX,
      std::vector<std::vector<bool> >(m_nY, std::vector<bool>(m_nZ, false)));

  std::ifstream infile;
  infile.open(filename.c_str(), std::ios::in);
  if (!infile) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }

  std::string line;
  unsigned int nLines = 0;
  bool bad = false;
  while (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;
    // Strip white space from beginning of line.
    ltrim(line);
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip comments.
    if (line[0] == '#') continue;
    if (line[0] == '/' && line[1] == '/') continue;
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;
    double val = 0;
    std::istringstream data;
    data.str(line);
    if (fmt == 1) {
      // "XY"
      double x, y;
      data >> x >> y;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      const double u = std::round((x - m_xMin) / m_dx);
      const double v = std::round((y - m_yMin) / m_dy);
      i = u < 0. ? 0 : static_cast<unsigned int>(u);
      j = v < 0. ? 0 : static_cast<unsigned int>(v);
      if (i >= m_nX) i = m_nX - 1;
      if (j >= m_nY) j = m_nY - 1;
    } else if (fmt == 2) {
      // "XYZ"
      double x, y, z;
      data >> x >> y >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      const double u = std::round((x - m_xMin) / m_dx);
      const double v = std::round((y - m_yMin) / m_dy);
      const double w = std::round((z - m_zMin) / m_dz);
      i = u < 0. ? 0 : static_cast<unsigned int>(u);
      j = v < 0. ? 0 : static_cast<unsigned int>(v);
      j = w < 0. ? 0 : static_cast<unsigned int>(w);
      if (i >= m_nX) i = m_nX - 1;
      if (j >= m_nY) j = m_nY - 1;
      if (k >= m_nZ) k = m_nZ - 1;
    } else if (fmt == 3) {
      // "IJ"
      data >> i >> j;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "indices");
        bad = true;
        break;
      }
    } else if (fmt == 4) {
      // "IJK"
      data >> i >> j >> k;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "indices");
        bad = true;
        break;
      }
    } else if (fmt == 5) {
      // "YXZ"
      double x, y, z;
      data >> y >> x >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      const double u = std::round((x - m_xMin) / m_dx);
      const double v = std::round((y - m_yMin) / m_dy);
      const double w = std::round((z - m_zMin) / m_dz);
      i = u < 0. ? 0 : static_cast<unsigned int>(u);
      j = v < 0. ? 0 : static_cast<unsigned int>(v);
      j = w < 0. ? 0 : static_cast<unsigned int>(w);
      if (i >= m_nX) i = m_nX - 1;
      if (j >= m_nY) j = m_nY - 1;
      if (k >= m_nZ) k = m_nZ - 1;
    }
    // Check the indices.
    if (i >= m_nX || j >= m_nY || k >= m_nZ) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Index (" << i << ", " << j << ", " << k
                << ") out of range.\n";
      continue;
    }
    if (isSet[i][j][k]) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Node (" << i << ", " << j << ", " << k
                << ") has already been set.\n";
      continue;
    }

    // Skip to the requested column.
    for (unsigned int i = 0; i < col - offset; ++i) {
      double dummy = 0.;
      data >> dummy;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, 
                   "column " + std::to_string(offset + i));
        break;
      }
    }
    if (data.fail()) {
      bad = true;
      break;
    }
    data >> val;

    if (data.fail()) {
      PrintError(m_className + "::LoadData", nLines, 
                 "column " + std::to_string(col));
      bad = true;
      break;
    }

    if (fmt == 1 || fmt == 3) {
      // Two-dimensional map
      for (unsigned int kk = 0; kk < m_nZ; ++kk) {
        tab[i][j][kk] = val;
        isSet[i][j][kk] = true;
      }
    } else {
      tab[i][j][k] = val;
      isSet[i][j][k] = true;
    }
    ++nValues;
  }
  infile.close();
  if (bad) return false;
  std::cout << m_className << "::LoadData:\n"
            << "    Read " << nValues << " values from " << filename << ".\n";
  unsigned int nExpected = m_nX * m_nY;
  if (fmt == 2 || fmt == 4 || fmt == 5) nExpected *= m_nZ;
  if (nExpected != nValues) {
    std::cerr << m_className << "::LoadData:\n"
              << "   Expected " << nExpected << " values.\n";
  }
  return true;
}

bool ComponentGrid::GetData(
    const double xi, const double yi, const double zi,
    const std::vector<std::vector<std::vector<double> > >& tab, double& val) {
  if (!m_hasMesh) {
    std::cerr << m_className << "::GetData: Mesh is not set.\n";
    return false;
  }

  // Reduce the point to the basic cell (in case of periodicity) and
  // check if it is inside the mesh.
  bool xMirrored = false;
  const double x =
      Reduce(xi, m_xMin, m_xMax, m_periodic[0], m_mirrorPeriodic[0], xMirrored);
  if (x < m_xMin || x > m_xMax) return false;
  bool yMirrored = false;
  const double y =
      Reduce(yi, m_yMin, m_yMax, m_periodic[1], m_mirrorPeriodic[1], yMirrored);
  if (y < m_yMin || y > m_yMax) return false;
  bool zMirrored = false;
  const double z =
      Reduce(zi, m_zMin, m_zMax, m_periodic[2], m_mirrorPeriodic[2], zMirrored);
  if (z < m_zMin || z > m_zMax) return false;

  // Get the indices.
  const double sx = (x - m_xMin) / m_dx;
  const double sy = (y - m_yMin) / m_dy;
  const double sz = (z - m_zMin) / m_dz;
  const unsigned int i0 = static_cast<unsigned int>(std::floor(sx));
  const unsigned int j0 = static_cast<unsigned int>(std::floor(sy));
  const unsigned int k0 = static_cast<unsigned int>(std::floor(sz));
  const double ux = sx - i0;
  const double uy = sy - j0;
  const double uz = sz - k0;
  const unsigned int i1 = std::min(i0 + 1, m_nX - 1);
  const unsigned int j1 = std::min(j0 + 1, m_nY - 1);
  const unsigned int k1 = std::min(k0 + 1, m_nZ - 1);
  const double vx = 1. - ux;
  const double vy = 1. - uy;
  const double vz = 1. - uz;
  const double n000 = tab[i0][j0][k0];
  const double n100 = tab[i1][j0][k0];
  const double n010 = tab[i0][j1][k0];
  const double n110 = tab[i1][j1][k0];
  const double n001 = tab[i0][j0][k1];
  const double n101 = tab[i1][j0][k1];
  const double n011 = tab[i0][j1][k1];
  const double n111 = tab[i1][j1][k1];

  if (m_debug) {
    std::cout << m_className << "::GetData: Interpolating at (" << xi
              << ", " << yi << ", " << zi << ").\n"
              << "    X: " << i0 << " (" << ux << ") - " << i1 << " (" << vx
              << ").\n"
              << "    Y: " << j0 << " (" << uy << ") - " << j1 << " (" << vy
              << ").\n"
              << "    Z: " << k0 << " (" << uz << ") - " << k1 << " (" << vz
              << ").\n";
  }
  val = ((n000 * vx + n100 * ux) * vy + (n010 * vx + n110 * ux) * uy) * vz +
        ((n001 * vx + n101 * ux) * vy + (n011 * vx + n111 * ux) * uy) * uz;

  return true;
}

bool ComponentGrid::ElectronAttachment(const double x, const double y,
                                       const double z, double& att) {
  // Make sure the map has been loaded.
  if (m_eattachment.empty()) {
    PrintNotReady(m_className + "::ElectronAttachment");
    return false;
  }
  return GetData(x, y, z, m_eattachment, att);
}

bool ComponentGrid::HoleAttachment(const double x, const double y,
                                   const double z, double& att) {
  // Make sure the map has been loaded.
  if (m_hattachment.empty()) {
    PrintNotReady(m_className + "::HoleAttachment");
    return false;
  }
  return GetData(x, y, z, m_hattachment, att);
}
}  // namespace Garfield
