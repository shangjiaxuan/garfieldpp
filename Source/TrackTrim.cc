#include <fstream>
#include <iostream>
#include <algorithm>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Numerics.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewBase.hh"
#include "Garfield/TrackTrim.hh"

namespace {

}

namespace Garfield {

TrackTrim::TrackTrim() : Track() { m_className = "TrackTrim"; }

bool TrackTrim::ReadRangeFile(const std::string& filename) {

  // TRMRER
  std::ifstream file;
  file.open(filename.c_str(), std::ios::in);
  if (file.fail()) {
    std::cerr << m_className << "::ReadRangeFile:\n"
              << "    Could not open " << filename << ".\n";
    return false;
  }

  const size_t size = 100;
  char line[size];
  while (file.getline(line, size, '\n')) {
    if (strstr(line, "TARGET")) {
      break;
    }
  }

  double xmin = 0.;
  double xmax = 0.;
  double rho = 0.;
  while (file.getline(line, size, '\n')) { // 20
    if (strstr(line, "Ion       Depth")) break;
    if (!strstr(line, "Layer #")) continue;
    char* token = std::strtok(line, "Layer #");
    token = std::strtok(nullptr, " ");
    if (!token) continue;
    unsigned int layer = std::atoi(token);
    if (layer == m_layer) {
      token = std::strtok(nullptr, "Depth=");
      if (token) {
        token = std::strtok(nullptr, " ");
        if (token) xmax = std::atof(token);
        continue;
      }
      token = std::strtok(nullptr, "Density =");
      if (token) {
        token = std::strtok(nullptr, " ");
        if (token) rho = std::atof(token);
      }
    } else if (m_layer > 0 && layer == m_layer - 1) {
      token = std::strtok(nullptr, "Depth=");
      if (token) {
        token = std::strtok(nullptr, " ");
        if (token) xmin = std::atof(token);
        continue;
      }
    }
  }
  if (rho <= 0.) {
    std::cerr << m_className << "::ReadRangeFile:\n"
              << "    WARNING: Density in file is zero.\n";
  } else {
    m_density = rho;
  }
  // TODO
  m_xmin = xmin;
  m_xmax = xmax;
  file.close();
  if (m_debug) {
    std::cout << m_className << "::ReadRangeFile:\n";
    std::printf("    Layer selected: %12d\n", m_layer);
    std::printf("    Layer starts at %12.5f A\n", m_xmin);
    std::printf("    Layer ends at   %12.5f A\n", m_xmax);
    std::printf("    Layer density:  %12.5f g/cm3\n", m_density);
  }
  return true;
}

bool TrackTrim::ReadXyzFile(const std::string& filename) {

  // TRMREE
  std::ifstream file(filename.c_str(), std::ios::in);
  if (file.fail()) {
    std::cerr << m_className << "::ReadXyzFile:\n"
              << "    Could not open " << filename << ".\n";
    return false;
  }

  bool first = true;
  const size_t size = 100;
  char line[size];
  while (file.getline(line, size, '\n')) { // 10
    if (!strstr(line, "Energy")) continue;
    if (first) {
      first = false;
      continue;
    }
    // Get the ion energy.

    // Skip heading lines.
    // Reading through until the desired ion is reached
    // Reading in the ion data
  }

  // First value of dE/dX may be too high - residual from previous layer.
  // if (TRMEMI[0] > 10 * TRMEMI[1]) {
  //   TRMEMI[0] = TRMEMI[1];
  // }
  file.close();
  return true;
}

bool TrackTrim::NewTrack(const double x0, const double y0, const double z0,
                         const double t0, const double dx0, const double dy0,
                         const double dz0) {
  // Verify that a sensor has been set.
  if (!m_sensor) {
    std::cerr << m_className << "::NewTrack: Sensor is not defined.\n";
    return false;
  }

  // Get the bounding box.
  double xmin = 0., ymin = 0., zmin = 0.;
  double xmax = 0., ymax = 0., zmax = 0.;
  if (!m_sensor->GetArea(xmin, ymin, zmin, xmax, ymax, zmax)) {
    std::cerr << m_className << "::NewTrack: Drift area is not set.\n";
    return false;
  } else if (x0 < xmin || x0 > xmax || y0 < ymin || y0 > ymax || z0 < zmin ||
             z0 > zmax) {
    std::cerr << m_className << "::NewTrack:\n"
              << "    Initial position outside bounding box.\n";
    return false;
  }

  // Make sure the initial position is inside an ionisable medium.
  Medium* medium = nullptr;
  if (!m_sensor->GetMedium(x0, y0, z0, medium)) {
    std::cerr << m_className << "::NewTrack: No medium at initial position.\n";
    return false;
  } else if (!medium->IsIonisable()) {
    std::cerr << m_className << "::NewTrack:\n"
              << "    Medium at initial position is not ionisable.\n";
    return false;
  }

  // Normalise and store the direction.
  const double normdir = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  double xdir = dx0;
  double ydir = dy0;
  double zdir = dz0;
  if (normdir < Small) {
    // Null vector. Sample the direction isotropically.
    RndmDirection(xdir, ydir, zdir);
  } else {
    // Normalise the direction vector.
    xdir /= normdir;
    ydir /= normdir;
    zdir /= normdir;
  }

  // Make sure all necessary parameters have been set.

  // Reset the cluster count
  m_currcluster = 0;
  m_clusters.clear();

  // Initial situation: starting position
  double x = x0;
  double y = y0;
  double z = z0;

  return true;
}

bool TrackTrim::GetCluster(double& xcls, double& ycls, double& zcls,
                           double& tcls, int& n, double& e, double& extra) {
  if (m_debug) {
    std::cout << m_className << "::GetCluster: Cluster " << m_currcluster
              << " of " << m_clusters.size() << "\n";
  }
  // Stop if we have exhausted the list of clusters.
  if (m_currcluster >= m_clusters.size()) return false;

  const auto& cluster = m_clusters[m_currcluster];
  xcls = cluster.x;
  ycls = cluster.y;
  zcls = cluster.z;
  tcls = cluster.t;

  n = cluster.electrons;
  e = cluster.ec;
  extra = cluster.kinetic;
  // Move to next cluster
  ++m_currcluster;
  return true;
}
}
