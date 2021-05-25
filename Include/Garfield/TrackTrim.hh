#ifndef G_TRACK_TRIM_H
#define G_TRACK_TRIM_H

#include <string>
#include <vector>

#include "Track.hh"

namespace Garfield {

/// Generate tracks based on TRIM output files.
///  - http://www.srim.org

class TrackTrim : public Track {
 public:
  /// Constructor
  TrackTrim();
  /// Destructor
  virtual ~TrackTrim() {}

  /// Load data from an EXYZ.txt file.
  bool ReadFile(const std::string& file, const unsigned int nIons = 0, 
                const unsigned int nSkip = 0);
  /// Print a summary of the available TRIM data.
  void Print();

  /// Set the W value [eV].
  void SetWorkFunction(const double w) { m_work = w; }
  /// Get the W value [eV].
  double GetWorkFunction() const { return m_work; }
  /// Set the Fano factor.
  void SetFanoFactor(const double f) { m_fano = f; }
  /// Get the Fano factor.
  double GetFanoFactor() const { return m_fano; }

  bool NewTrack(const double x0, const double y0, const double z0,
                const double t0, const double dx0, const double dy0,
                const double dz0) override;
  bool GetCluster(double& xcls, double& ycls, double& zcls,
                  double& tcls, int& n, double& e, double& extra) override;

 protected:
  /// Work function [eV].
  double m_work = -1.;
  /// Fano factor [-].
  double m_fano = -1.;

  /// Projectile energy [eV].
  double m_ekin = 0.;
  /// List of tracks. 
  std::vector<std::vector<std::array<float, 5> > > m_ions;
  /// Index of the current track.
  size_t m_ion = 0;

  struct Cluster {
    double x, y, z, t;  ///< Cluster location and time
    double ec;          ///< Energy spent to make the cluster
    double ekin;        ///< Ion energy when cluster was created
    int electrons;      ///< Number of electrons in this cluster
  };
  /// Clusters on the current track.
  std::vector<Cluster> m_clusters;
  /// Index of the next cluster to be returned.
  size_t m_cluster = 0;

  void AddIon(const std::vector<float>& x, const std::vector<float>& y,
              const std::vector<float>& z, const std::vector<float>& dedx,
              const std::vector<float>& ekin);
};
}

#endif
