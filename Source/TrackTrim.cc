#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackTrim.hh"

namespace {

std::vector<std::string> tokenize(const std::string& line) {

  std::vector<std::string> words;
  std::istringstream ss(line);
	for (std::string word; ss >> word;) {
		words.push_back(word);
  }
	return words;
}

}

namespace Garfield {

TrackTrim::TrackTrim() : Track() { m_className = "TrackTrim"; }

bool TrackTrim::ReadFile(const std::string& filename, 
                         const unsigned int nIons, const unsigned int nSkip) {

  // TRMREE - Reads the TRIM EXYZ file.

  // Reset.
  m_ions.clear();
  m_ion = 0;
  m_clusters.clear();
  m_cluster = 0;

  std::ifstream infile;
  infile.open(filename.c_str(), std::ios::in);
  if (infile.fail()) {
    std::cerr << m_className << "::ReadFile:\n"
              << "    Unable to open the EXYZ file (" << filename << ").\n";
    return false;
  }

  constexpr double Angstrom = 1.e-8;
  unsigned int nRead = 0;
  unsigned int ionNumber = 0;
  bool header = true;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> dedx;
  std::vector<float> ekin;
  for (std::string line; std::getline(infile, line);) {
    if (line.find("------- ") != std::string::npos) {
      // Reached the end of the header.
      header = false;
      continue;
    } else if (header) {
      if (line.find("Ion Data: ") != std::string::npos) {
        // Read the next line.
        std::getline(infile, line);
        auto words = tokenize(line);
        if (words.size() >= 3) {
          m_particleName = words[0];
          auto pos = words[2].find("keV");
          if (pos != std::string::npos) {
            m_ekin = 1.e3 * std::stod(words[2].substr(0, pos));
          }
        } 
      } 
      // Otherwise, skip the header.
      continue;
    } 
    auto words = tokenize(line);
    if (words.size() < 6) {
      std::cerr << m_className << "::ReadFile: Unexpected line:\n"
                << line << "\n";
      continue;
    }
    if (ionNumber != std::stoul(words[0])) {
      // New ion.
      if (ionNumber > 0) {
        if (nRead >= nSkip) AddIon(x, y, z, dedx, ekin);
        x.clear();
        y.clear();
        z.clear();
        dedx.clear(); 
        ekin.clear();
        ++nRead;
        // Stop if we are done reading the requested number of ions.
        if (nIons > 0 && nRead >= nIons) break;
      }
      ionNumber = std::stoi(words[0]);
    }
    if (nRead < nSkip) continue;
    // Convert coordinates to cm.
    x.push_back(std::stof(words[2]) * Angstrom);
    y.push_back(std::stof(words[3]) * Angstrom);
    z.push_back(std::stof(words[4]) * Angstrom);
    // Convert stopping power from eV/A to eV/cm.
    dedx.push_back(std::stof(words[5]) / Angstrom);
    // Convert ion energy from keV to eV.
    ekin.push_back(std::stof(words[1]) * 1.e3);
  }
  infile.close();
  AddIon(x, y, z, dedx, ekin);
  std::cout << m_className << "::ReadFile: Read energy vs position for " 
            << m_ions.size() << " ions.\n";
  return true;
}

void TrackTrim::AddIon(const std::vector<float>& x,
                       const std::vector<float>& y,
                       const std::vector<float>& z,
                       const std::vector<float>& dedx, 
                       const std::vector<float>& ekin) {

  const size_t nPoints = x.size();
  if (nPoints < 2) return;
  std::vector<std::array<float, 5> > path;
  for (size_t i = 0; i < nPoints; ++i) {
    float eloss = 0.;
    if (i < nPoints - 1) {
      const float dx = x[i + 1] - x[i];
      const float dy = y[i + 1] - y[i];
      const float dz = z[i + 1] - z[i];
      const float step = sqrt(dx * dx + dy * dy + dz * dz);
      if (i == 0 && dedx[i] > 10. * dedx[i + 1]) {
        eloss = step * dedx[i + 1];
      } else { 
        eloss = step * dedx[i];
      }
      const float dekin = ekin[i] - ekin[i + 1];
      if (dekin > 0.) eloss = std::min(eloss, dekin); 
    }
    path.push_back({x[i], y[i], z[i], eloss, ekin[i]});
  }
  m_ions.push_back(path);
}

void TrackTrim::Print() {
  std::cout << m_className << "::Print:\n";
  if (m_ions.empty()) {
    std::cerr << "    No TRIM data present.\n";
    return;
  }
  std::cout << "    Projectile: " << m_particleName << ", "
            << m_ekin * 1.e-3 << " keV\n"
            << "    Number of tracks: " << m_ions.size() << "\n"
            << "    Work function: " << m_work << " eV\n"
            << "    Fano factor: " << m_fano << "\n";
}

bool TrackTrim::NewTrack(const double x0, const double y0, const double z0,
    const double t0, const double dx0, const double dy0, const double dz0) {

  // TRMGEN - Generates TRIM clusters

  if (m_ions.empty()) {
    std::cerr << m_className << "::NewTrack: No TRIM data present.\n";
    return false;
  }

  if (m_ion >= m_ions.size()) {
    // Rewind.
    std::cout << m_className << "::NewTrack: Rewinding.\n";
    m_ion = 0;
  }

  // Verify that a sensor has been set.
  if (!m_sensor) {
    std::cerr << m_className << "::NewTrack: Sensor is not defined.\n";
    return false;
  }

  // Normalise and store the direction.
  const double d0 = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  double dx = dx0;
  double dy = dy0;
  double dz = dz0;
  if (d0 < Small) {
    if (m_debug) {
      std::cout << m_className << "::NewTrack: Randomizing initial direction.\n";
    }
    // Null vector. Sample the direction isotropically.
    RndmDirection(dx, dy, dz);
  } else {
    // Normalise the direction vector.
    dx /= d0;
    dy /= d0;
    dz /= d0;
  }
  const double dt = sqrt(dx * dx + dy * dy);
  double phi = 0.;
  double theta = 0.; 
  if (dt > 0.) {
    phi = atan2(dy, dx);
    theta = atan2(dz, dt);
  } else {
    theta = dz < 0. ? -HalfPi : HalfPi;
  }
  const double ctheta = cos(theta);
  const double stheta = sin(theta);
  const double cphi = cos(phi);
  const double sphi = sin(phi);

  // Make sure all necessary parameters have been set.
  if (m_work < Small) {
    std::cerr << m_className << "::NewTrack: Work function not set.\n";
    return false;
  }
 
  // Plot.
  if (m_viewer) PlotNewTrack(x0, y0, z0);
 
  // Reset the cluster count.
  m_cluster = 0;
  m_clusters.clear();

  // Pool of unused energy
  double epool = 0.0;

  const auto& path = m_ions[m_ion];
  const size_t nPoints = path.size();
  for (size_t i = 1; i < nPoints; ++i) {
    const double x = path[i][0];
    const double y = path[i][1];
    const double z = path[i][2];
    Cluster cluster;
    cluster.x = x0 + cphi * ctheta * x - sphi * y - cphi * stheta * z;
    cluster.y = y0 + sphi * ctheta * x + cphi * y - sphi * stheta * z;
    cluster.z = z0 + stheta * x + ctheta * z;
    // Is this point inside an ionisable medium?
    Medium* medium = nullptr;
    if (!m_sensor->GetMedium(cluster.x, cluster.y, cluster.z, medium)) {
      continue;
    } 
    if (!medium || !medium->IsIonisable()) continue;
    cluster.t = t0;
    double eloss = path[i - 1][3];
    if (m_fano < Small) {
      // No fluctuations.
      cluster.electrons = int((eloss + epool) / m_work);
      cluster.ec = m_work * cluster.electrons;
    } else {
      double ec = eloss + epool;
      cluster.electrons = 0;
      cluster.ec = 0.0;
      while (true) {
        const double er = RndmHeedWF(m_work, m_fano);
        if (er > ec) break;
        cluster.electrons++;
        cluster.ec += er;
        ec -= er;
      }
    }
    cluster.ekin = path[i][4];
    epool += eloss - cluster.ec;
    m_clusters.push_back(std::move(cluster));
    if (m_viewer) PlotCluster(cluster.x, cluster.y, cluster.z);
  }
  // Move to the next ion in the list.
  ++m_ion;
  return true;
}

bool TrackTrim::GetCluster(double& xcls, double& ycls, double& zcls,
                           double& tcls, int& n, double& e, double& extra) {
  if (m_debug) {
    std::cout << m_className << "::GetCluster: Cluster " << m_cluster
              << " of " << m_clusters.size() << "\n";
  }
  // Stop if we have exhausted the list of clusters.
  if (m_cluster >= m_clusters.size()) return false;

  const auto& cluster = m_clusters[m_cluster];
  xcls = cluster.x;
  ycls = cluster.y;
  zcls = cluster.z;
  tcls = cluster.t;

  n = cluster.electrons;
  e = cluster.ec;
  extra = cluster.ekin;
  // Move to the next cluster.
  ++m_cluster;
  return true;
}
}
