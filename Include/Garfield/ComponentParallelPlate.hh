#ifndef G_COMPONENT_PP_H
#define G_COMPONENT_PP_H

#include <string>

#include "Component.hh"
#include "Medium.hh"

namespace Garfield {

/// Component for parallel-plate geometries.

class GARFIELD_IMPORTEXPROT ComponentParallelPlate : public Component {
 public:
  /// Constructor
  ComponentParallelPlate();
  /// Destructor
  ~ComponentParallelPlate() {}

  /** Define the geometry.
   * \param g size of the gap along positive \f$z\f$.
   * \param b thickness of the resistive layer along negative \f$z\f$.
   * \param eps relative permittivity of the resistive layer.
   * \param sigma conductivity of the resistive layer (must be larger then zero,
   * otherwise do not pass it in the function).
   * \param V applied potential difference between the parallel plates.
   */
  void Setup(double g, double b, double eps, double v, double sigma = 0.);

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;

  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  double DelayedWeightingPotential(const double x, const double y,
                                   const double z, const double t,
                                   const std::string& label) override;

  void DelayedWeightingField(const double x, const double y, const double z,
                             const double t, double& wx, double& wy, double& wz,
                             const std::string& label) override;

  bool GetVoltageRange(double& vmin, double& vmax) override;

  /** Add a pixel electrode.
   * \param x,z position of the center of the electrode in the xy-plane.
   * \param lx width in the along \f$x\f$.
   * \param ly width in the along \f$y\f$.
   * \param label give name using a string.
   */
  void AddPixel(double x, double y, double lx, double ly,
                const std::string& label);
  /// Add strip electrode.
  void AddStrip(double x, double lx, const std::string& label);

  /// Add plane electrode, if you want to read the signal from the cathode set
  /// the second argument to false.
  void AddPlane(const std::string& label, bool anode = true);

  // Setting the medium
  void SetMedium(Medium* medium) { m_medium = medium; }

  // This will calculate the electrode's time-dependent weighting potential on
  // the specified grid.
  void SetWeightingPotentialGrid(const std::string& label, const double xmin,
                                 const double xmax, const double xsteps,
                                 const double ymin, const double ymax,
                                 const double ysteps, const double zmin,
                                 const double zmax, const double zsteps,
                                 const double tmin, const double tmax,
                                 const double tsteps);

  // This will calculate all electrodes time-dependent weighting potential on
  // the specified grid.
  void SetWeightingPotentialGrids(const double xmin, const double xmax,
                                  const double xsteps, const double ymin,
                                  const double ymax, const double ysteps,
                                  const double zmin, const double zmax,
                                  const double zsteps, const double tmin,
                                  const double tmax, const double tsteps);

  Medium* GetMedium(const double x, const double y, const double z) override;

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
                      double& xmax, double& ymax, double& zmax) override;
 private:
  static constexpr double m_precision = 1.e-30;
  // Size of the gap.
  double m_g = 0.;
  // Thickness of the resistive element.
  double m_b = 0.;
  // Applied voltage on the electrode to
  // calculate the weighting potential.
  static constexpr double m_Vw = 1.;
  double m_eps = 1.;
  double m_eps0 = 8.85418782e-3;
  // Voltage difference between the parallel plates.
  double m_V = 0.;
  // Electric field in the gap.
  double m_ezg = 0.;
  // Electric field in the resistive layer.
  double m_ezb = 0.;
  // Conductivity of the resistive layer.
  double m_sigma = 0.;

  Medium* m_medium = nullptr;

  /// Structure that captures the information of the electrodes under study
  struct Electrode {
    std::string label;                     ///< Label.
    int ind = structureelectrode::NotSet;  ///< Readout group.
    double xpos, ypos;                     ///< Coordinates in x/y.
    double lx, ly;                         ///< Dimensions in the x-y plane.
    double flip = 1;                       ///< Dimensions in the x-y plane.

    bool m_usegrid = false;
    std::vector<std::vector<std::vector<double>>> gridPromptV;
    std::vector<std::vector<std::vector<std::vector<double>>>> gridDelayedV;

    double gridXSteps = 0;
    double gridYSteps = 0;
    double gridZSteps = 0;
    double gridTSteps = 0;

    double gridX0 = 0;
    double gridY0 = 0;
    double gridZ0 = 0;
    double gridT0 = 0;

    double gridXStepSize = 0;
    double gridYStepSize = 0;
    double gridZStepSize = 0;
    double gridTStepSize = 0;
  };

  enum fieldcomponent { xcomp = 0, ycomp, zcomp };

  /// Possible readout groups
  enum structureelectrode { NotSet = -1, Plane, Strip, Pixel };

  // Vectors storing the readout electrodes.
  std::vector<std::string> m_readout;
  std::vector<Electrode> m_readout_p;

  // Functions that calculate the electric field and potential
  double IntegrateField(const Electrode& el, int comp, const double x,
                        const double y, const double z);

  double IntegrateDelayedField(const Electrode& el, int comp, const double x,
                               const double y, const double z, const double t);

  double IntegratePromptPotential(const Electrode& el, const double x,
                                  const double y, const double z);
  double IntegrateDelayedPotential(const Electrode& el, const double x,
                                   const double y, const double z,
                                   const double t);

  void CalculateDynamicalWeightingPotential(const Electrode& el);

  double FindWeightingPotentialInGrid(const Electrode& el, const double x,
                                      const double y, const double z);

  double FindDelayedWeightingPotentialInGrid(const Electrode& el,
                                             const double x, const double y,
                                             const double z, const double t);

  double FindWeightFactor(const Electrode& el, const double dx, const double dy,
                          const double dz);

  double FindWeightFactor(const Electrode& el, const double dx, const double dy,
                          const double dz, const double dt);

  void UpdatePeriodicity() override;
  void Reset() override;
};
}  // namespace Garfield
#endif
