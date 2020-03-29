#ifndef G_VIEW_DRIFT
#define G_VIEW_DRIFT

#include <memory>
#include <string>
#include <vector>
#include <array>
#include <utility>

#include <TView.h>

#include "ViewBase.hh"

namespace Garfield {

/// Visualize drift lines and tracks.

class ViewDrift : public ViewBase {
 public:
  /// Constructor.
  ViewDrift();
  /// Destructor.
  ~ViewDrift();

  /// Delete existing drift lines, tracks and markers.
  void Clear();

  /// Draw the drift lines.
  void Plot(const bool twod = false, const bool axis = true);

  /// Set the size of the cluster markers (see TAttMarker).
  void SetClusterMarkerSize(const double size);
  /// Set the size of the collision markers (see TAttMarker).
  void SetCollisionMarkerSize(const double size);

  // Functions used by the transport classes.
  void NewElectronDriftLine(const unsigned int np, int& id, const float x0,
                            const float y0, const float z0);
  void NewHoleDriftLine(const unsigned int np, int& id, const float x0,
                        const float y0, const float z0);
  void NewIonDriftLine(const unsigned int np, int& id, const float x0,
                       const float y0, const float z0);
  void NewChargedParticleTrack(const unsigned int np, int& id, const float x0,
                               const float y0, const float z0);

  void SetDriftLinePoint(const unsigned int iL, const unsigned int iP,
                         const float x, const float y, const float z);
  void AddDriftLinePoint(const unsigned int iL, const float x, const float y,
                         const float z);
  void SetTrackPoint(const unsigned int iL, const unsigned int iP,
                     const float x, const float y, const float z);
  void AddTrackPoint(const unsigned int iL, const float x, const float y,
                     const float z);
  void AddExcitation(const float x, const float y, const float z);
  void AddIonisation(const float x, const float y, const float z);
  void AddAttachment(const float x, const float y, const float z);

  void AddPhoton(const float x0, const float y0, const float z0,
                 const float x1, const float y1, const float z1);

  friend class ViewFEMesh;

 private:

  // View
  std::unique_ptr<TView> m_view;

  enum class Particle {
    Electron,
    Hole,
    Ion
  };
  std::vector<std::pair<std::vector<std::array<float, 3> >,
                        Particle> > m_driftLines;

  std::vector<std::vector<std::array<float, 3> > > m_tracks;
  std::vector<std::array<std::array<float, 3>, 2> > m_photons;

  std::vector<std::array<float, 3> > m_exc;
  std::vector<std::array<float, 3> > m_ion;
  std::vector<std::array<float, 3> > m_att;

  double m_markerSizeCluster = 1.;
  double m_markerSizeCollision = 0.5;

  void Plot2d(const bool axis);
  void Plot3d(const bool axis);
};
}
#endif
