#ifndef G_TRACK_TRIM_H
#define G_TRACK_TRIM_H
#include <vector>
#include "Track.hh"

namespace Garfield {

/// Generate tracks based on events simulated by the TRIM program.
/// Based on the Garfield TRIMCAT interface contributed by J. Butterworth. 

class TrackTrim : public Track {
 public:
  /// Constructor
  TrackTrim();
  /// Destructor
  virtual ~TrackTrim() {}

  /// Set the W value [eV].
  void SetWorkFunction(const double w) { m_work = w; }
  /// Get the W value [eV].
  double GetWorkFunction() const { return m_work; }
  /// Set the Fano factor.
  void SetFanoFactor(const double f) { m_fano = f; }
  /// Get the Fano factor.
  double GetFanoFactor() const { return m_fano; }

  /// Load data.
  bool ReadFile(const std::string& file);
  void Print();
  /// Make a plot of the electromagnetic, hadronic, and total energy loss
  /// as function of the projectile energy.
  void PlotEnergyLoss();
  /// Make a plot of the projected range as function of the projectile energy.
  void PlotRange();
  /// Make a plot of the transverse and longitudinal straggling as function
  /// of the projectile energy. 
  void PlotStraggling();

  bool NewTrack(const double x0, const double y0, const double z0,
                const double t0, const double dx0, const double dy0,
                const double dz0) override;
  bool GetCluster(double& xcls, double& ycls, double& zcls,
                  double& tcls, int& n, double& e, double& extra) override;

 protected:
  /// Work function [eV]
  double m_work = 0.;
  /// Fano factor [-]
  double m_fano = 0.17;
  /// Density
  double m_density = 0.;

  /// Layer number of the active volume
  unsigned int m_layer = 0;
  /// Minimum x-coordinate of the active layer.
  double m_xmin = 0.;
  /// Maximum x-coordinate of the active layer.
  double m_xmax = 0.;

  // ITRIM   : The number of the ion
  // NCTRIM  : The number of clusters

  // TRMEMI  : EM energy loss
  // TRMHDI  : HD energy loss
  // TRMTGD  : TRIM target depth
  // TRMIOE  : TRIM ion energy
  // TRMY  : Corresponding Y location of the ion
  // TRMZ  : Corresponding Z location of the ion

  std::vector<double> m_emloss;
  /// Hadronic energy loss [MeV cm2/g]
  std::vector<double> m_hdloss;

  struct Cluster {
    double x, y, z, t;  ///< Cluster location and time
    double ec;          ///< Energy spent to make the cluster
    double kinetic;     ///< Ion energy when cluster was created
    int electrons;      ///< Number of electrons in this cluster
  };
  std::vector<Cluster> m_clusters;
  /// Index of the next cluster to be returned.
  size_t m_currcluster = 0;

  bool ReadRangeFile(const std::string& filename);
  bool ReadXyzFile(const std::string& filename);
};
}

#endif
