#ifndef G_AVALANCHE_MICROSCOPIC_H
#define G_AVALANCHE_MICROSCOPIC_H

#include <string>
#include <vector>

#include <TH1.h>

#include "GarfieldConstants.hh"
#include "Sensor.hh"
#include "ViewDrift.hh"

namespace Garfield {

/// Calculate electron drift lines and avalanches using microscopic tracking.

class AvalancheMicroscopic {
 public:
  /// Constructor
  AvalancheMicroscopic();
  /// Destructor
  ~AvalancheMicroscopic() {}

  /// Set the sensor.
  void SetSensor(Sensor* sensor);

  /// Switch on drift line plotting.
  void EnablePlotting(ViewDrift* view);
  void DisablePlotting();
  /// Draw a marker at every excitation.
  void EnableExcitationMarkers() { m_plotExcitations = true; }
  void DisableExcitationMarkers() { m_plotExcitations = false; }
  /// Draw a marker at every ionising collision.
  void EnableIonisationMarkers() { m_plotIonisations = true; }
  void DisableIonisationMarkers() { m_plotIonisations = false; }
  /// Draw a marker at every attachment.
  void EnableAttachmentMarkers() { m_plotAttachments = true; }
  void DisableAttachmentMarkers() { m_plotAttachments = false; }

  /// Switch on calculation of induced currents
  void EnableSignalCalculation() { m_useSignal = true; }
  void DisableSignalCalculation() { m_useSignal = false; }

  /// Switch on calculation of the total induced charge
  void EnableInducedChargeCalculation() { m_useInducedCharge = true; }
  void DisableInducedChargeCalculation() { m_useInducedCharge = false; }

  /// Switch on filling histograms for electron energy distribution
  void EnableElectronEnergyHistogramming(TH1* histo);
  void DisableElectronEnergyHistogramming() { m_histElectronEnergy = nullptr; }
  void EnableHoleEnergyHistogramming(TH1* histo);
  void DisableHoleEnergyHistogramming() { m_histHoleEnergy = nullptr; }

  /** Fill histograms of the distance between successive collisions.
    * \param histo
             pointer to the histogram to be filled
    * \param opt
             direction ('x', 'y', 'z', 'r')
    */
  void SetDistanceHistogram(TH1* histo, const char opt = 'r');
  void EnableDistanceHistogramming(const int type);
  /// Switch on distance distribution histograms for a given collision type.
  void DisableDistanceHistogramming(const int type);
  void DisableDistanceHistogramming();
  /// Fill histograms of the energy of electrons emitted in ionising collisions.
  void EnableSecondaryEnergyHistogramming(TH1* histo);
  void DisableSecondaryEnergyHistogramming() { m_histSecondary = nullptr; }

  /// Switch on storage of drift lines.
  void EnableDriftLines() { m_useDriftLines = true; }
  void DisableDriftLines() { m_useDriftLines = false; }

  /** Switch on photon transport.
    * \remark This feature has not been tested thoroughly. */
  void EnablePhotonTransport() { m_usePhotons = true; }
  void DisablePhotonTransport() { m_usePhotons = false; }

  /// Switch on stepping according to band structure E(k), for semiconductors.
  void EnableBandStructure() { m_useBandStructureDefault = true; }
  void DisableBandStructure() { m_useBandStructureDefault = false; }

  /// Switch on update of coordinates for null-collision steps (default: off).
  void EnableNullCollisionSteps() { m_useNullCollisionSteps = true; }
  void DisableNullCollisionSteps() { m_useNullCollisionSteps = false; }

  /** Set a (lower) energy threshold for electron transport.
    * This can be useful for simulating delta electrons. */
  void SetElectronTransportCut(const double cut) { m_deltaCut = cut; }
  double GetElectronTransportCut() const { return m_deltaCut; }

  /// Set an energy threshold for photon transport.
  void SetPhotonTransportCut(const double cut) { m_gammaCut = cut; }
  double GetPhotonTransportCut() const { return m_gammaCut; }

  /** Set a max. avalanche size (i. e. ignore ionising collisions
      once this size has been reached). */
  void EnableAvalancheSizeLimit(const unsigned int size) { m_sizeCut = size; }
  void DisableAvalancheSizeLimit() { m_sizeCut = 0; }
  int GetAvalancheSizeLimit() const { return m_sizeCut; }

  /// Enable magnetic field in stepping algorithm (default: off).
  void EnableMagneticField() { m_useBfield = true; }
  void DisableMagneticField() { m_useBfield = false; }

  /// Set number of collisions to be skipped for plotting
  void SetCollisionSteps(const unsigned int n) { m_nCollSkip = n; }

  /// Define a time interval (only carriers inside the interval are simulated).
  void SetTimeWindow(const double t0, const double t1);
  void UnsetTimeWindow();

  /// Return the number of electrons and ions in the avalanche.
  void GetAvalancheSize(int& ne, int& ni) const {
    ne = m_nElectrons;
    ni = m_nIons;
  }
  void GetAvalancheSize(int& ne, int& nh, int& ni) const {
    ne = m_nElectrons;
    nh = m_nHoles;
    ni = m_nIons;
  }

  /** Return the number of electron trajectories in the last
    * simulated avalanche (including captured electrons). */
  unsigned int GetNumberOfElectronEndpoints() const {
    return m_endpointsElectrons.size();
  }
  /** Return the coordinates and time of start and end point of a given
    * electron drift line.
    * \param i index of the drift line
    * \param x0,y0,z0,t0 coordinates and time of the starting point
    * \param x1,y1,z1,t1 coordinates and time of the end point
    * \param e0,e1 initial and final energy
    * \param status status code (see GarfieldConstants.hh)
    */
  void GetElectronEndpoint(const unsigned int i, double& x0, double& y0,
                           double& z0, double& t0, double& e0, double& x1,
                           double& y1, double& z1, double& t1, double& e1,
                           int& status) const;
  void GetElectronEndpoint(const unsigned int i, double& x0, double& y0,
                           double& z0, double& t0, double& e0, double& x1,
                           double& y1, double& z1, double& t1, double& e1,
                           double& dx1, double& dy1, double& dz1,
                           int& status) const;
  unsigned int GetNumberOfElectronDriftLinePoints(
      const unsigned int i = 0) const;
  unsigned int GetNumberOfHoleDriftLinePoints(const unsigned int i = 0) const;
  void GetElectronDriftLinePoint(double& x, double& y, double& z, double& t,
                                 const int ip,
                                 const unsigned int iel = 0) const;
  void GetHoleDriftLinePoint(double& x, double& y, double& z, double& t,
                             const int ip, const unsigned int iel = 0) const;

  unsigned int GetNumberOfHoleEndpoints() const {
    return m_endpointsHoles.size();
  }
  void GetHoleEndpoint(const unsigned int i, double& x0, double& y0, double& z0,
                       double& t0, double& e0, double& x1, double& y1,
                       double& z1, double& t1, double& e1, int& status) const;

  unsigned int GetNumberOfPhotons() const { return m_photons.size(); }
  // Status codes:
  //   -2: photon absorbed by gas molecule
  void GetPhoton(const unsigned int i, double& e, double& x0, double& y0,
                 double& z0, double& t0, double& x1, double& y1, double& z1,
                 double& t1, int& status) const;

  /** Calculate an electron drift line.
    * \param x0,y0,z0,t0 starting point of the electron
    * \param e0 initial energy of the electron
    * \param dx0,dy0,dz0 initial direction of the electron
    * If the initial direction is not specified, it is sampled randomly.
    * Secondary electrons are not transported. */
  bool DriftElectron(const double x0, const double y0, const double z0,
                     const double t0, const double e0, const double dx0 = 0.,
                     const double dy0 = 0., const double dz0 = 0.);

  /// Calculate an avalanche initiated by a given electron.
  bool AvalancheElectron(const double x0, const double y0, const double z0,
                         const double t0, const double e0,
                         const double dx0 = 0., const double dy0 = 0.,
                         const double dz0 = 0.);

  /// Set a user handling procedure. This function is called at every step.
  void SetUserHandleStep(void (*f)(double x, double y, double z, double t,
                                   double e, double dx, double dy, double dz,
                                   bool hole));
  void UnsetUserHandleStep() { m_userHandleStep = nullptr; }
  /// Set a user handling procedure, to be called at every (real) collision.
  void SetUserHandleCollision(void (*f)(double x, double y, double z, double t,
                                        int type, int level, Medium* m,
                                        double e0, double e1,
                                        double dx0, double dy0, double dz0,
                                        double dx1, double dy1, double dz1));
  void UnsetUserHandleCollision() { m_userHandleCollision = nullptr; }
  /// Set a user handling procedure, to be called at every attachment.
  void SetUserHandleAttachment(void (*f)(double x, double y, double z, double t,
                                         int type, int level, Medium* m));
  void UnsetUserHandleAttachment() { m_userHandleAttachment = nullptr; }
  /// Set a user handling procedure, to be called at every inelastic collision.
  void SetUserHandleInelastic(void (*f)(double x, double y, double z, double t,
                                        int type, int level, Medium* m));
  void UnsetUserHandleInelastic() { m_userHandleInelastic = nullptr; }
  /// Set a user handling procedure, to be called at every ionising collision.
  void SetUserHandleIonisation(void (*f)(double x, double y, double z, double t,
                                         int type, int level, Medium* m));
  void UnsetUserHandleIonisation() { m_userHandleIonisation = nullptr; }

  /// Switch on debugging messages.
  void EnableDebugging() { m_debug = true; }
  void DisableDebugging() { m_debug = false; }

 private:
  std::string m_className = "AvalancheMicroscopic";

  Sensor* m_sensor = nullptr;

  struct point {
    double x, y, z, t;
  };

  struct Electron {
    int status;                    //< Status.
    bool hole;                     //< Electron or hole.
    double x0, y0, z0, t0;         //< Starting point and time.
    double e0;                     //< Initial kinetic energy.
    int band;                      //< Band.
    double x, y, z, t;             //< Current position and time.
    double kx, ky, kz;             //< Current direction/wave vector.
    double energy;                 //< Current kinetic energy.
    std::vector<point> driftLine;  //< Drift line.
    double xLast, yLast, zLast;    //< Previous position.
  };
  std::vector<Electron> m_endpointsElectrons;
  std::vector<Electron> m_endpointsHoles;

  struct photon {
    int status;             //< Status
    double energy;          //< Energy
    double x0, y0, z0, t0;  //< Starting point and time.
    double x1, y1, z1, t1;  //< End point and time.
  };
  std::vector<photon> m_photons;

  /// Number of electrons produced
  int m_nElectrons = 0;
  /// Number of holes produced
  int m_nHoles = 0;
  /// Number of ions produced
  int m_nIons = 0;

  bool m_usePlotting = false;
  ViewDrift* m_viewer = nullptr;
  bool m_plotExcitations = true;
  bool m_plotIonisations = true;
  bool m_plotAttachments = true;

  TH1* m_histElectronEnergy = nullptr;
  TH1* m_histHoleEnergy = nullptr;
  TH1* m_histDistance = nullptr;
  char m_distanceOption = 'r';
  std::vector<int> m_distanceHistogramType;

  TH1* m_histSecondary = nullptr;

  bool m_useSignal = false;
  bool m_useInducedCharge = false;
  bool m_useDriftLines = false;
  bool m_usePhotons = false;
  bool m_useBandStructureDefault = true;
  bool m_useNullCollisionSteps = false;
  bool m_useBfield = false;

  // Rotation matrices
  double m_rb11 = 1., m_rb12 = 0., m_rb13 = 0.;
  double m_rb21 = 0., m_rb22 = 1., m_rb23 = 0.;
  double m_rb31 = 0., m_rb32 = 0., m_rb33 = 1.;
  double m_rx22 = 1., m_rx23 = 0., m_rx32 = 0., m_rx33 = 1.;

  // Transport cuts
  double m_deltaCut = 0.;
  double m_gammaCut = 0.;

  // Max. avalanche size
  unsigned int m_sizeCut = 0;

  unsigned int m_nCollSkip = 100;

  bool m_hasTimeWindow = false;
  double m_tMin = 0.;
  double m_tMax = 0.;

  // User procedures
  void (*m_userHandleStep)(double x, double y, double z, double t, double e,
                           double dx, double dy, double dz,
                           bool hole) = nullptr;
  void (*m_userHandleCollision)(double x, double y, double z, double t,
                                int type, int level, Medium* m,
                                double e0, double e1,
                                double dx0, double dy0, double dz0,
                                double dx1, double dy1, double dz1) = nullptr;
  void (*m_userHandleAttachment)(double x, double y, double z, double t,
                                 int type, int level, Medium* m) = nullptr;
  void (*m_userHandleInelastic)(double x, double y, double z, double t,
                                int type, int level, Medium* m) = nullptr;
  void (*m_userHandleIonisation)(double x, double y, double z, double t,
                                 int type, int level, Medium* m) = nullptr;

  // Switch on/off debugging messages
  bool m_debug = false;

  // Electron transport
  bool TransportElectron(const double x0, const double y0, const double z0,
                         const double t0, double e0, const double dx0,
                         const double dy0, const double dz0, const bool aval,
                         bool hole);
  // Photon transport
  void TransportPhoton(const double x, const double y, const double z,
                       const double t, const double e,
                       std::vector<Electron>& stack);

  void ComputeRotationMatrix(const double bx, const double by, const double bz,
                             const double bmag, const double ex,
                             const double ey, const double ez);

  void RotateGlobal2Local(double& dx, double& dy, double& dz) const;
  void RotateLocal2Global(double& dx, double& dy, double& dz) const;

  static bool IsInactive(const Electron& item) {
    return item.status == StatusLeftDriftMedium ||
           item.status == StatusBelowTransportCut ||
           item.status == StatusOutsideTimeWindow ||
           item.status == StatusLeftDriftArea || item.status == StatusAttached;
  }
  void Update(std::vector<Electron>::iterator it, const double x,
              const double y, const double z, const double t,
              const double energy, const double kx, const double ky,
              const double kz, const int band);
  void AddToEndPoints(const Electron& item, const bool hole) {
    if (hole) {
      m_endpointsHoles.push_back(item);
    } else {
      m_endpointsElectrons.push_back(item);
    }
  }

  /// Add a new electron/hole (with random direction) to a container.
  void AddToStack(const double x, const double y, const double z,
                  const double t, const double energy, const bool hole,
                  std::vector<Electron>& container) const;
  /// Add a new electron/hole to a container.
  void AddToStack(const double x, const double y, const double z,
                  const double t, const double energy, const double dx,
                  const double dy, const double dz, const int band,
                  const bool hole, std::vector<Electron>& container) const;
  void Terminate(double x0, double y0, double z0, double t0, double& x1,
                 double& y1, double& z1, double& t1);
};
}

#endif
