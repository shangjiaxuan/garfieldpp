#ifndef G_COMPONENT_PP_H
#define G_COMPONENT_PP_H

#include <string>

#include "Component.hh"
#include "Medium.hh"

namespace Garfield {

/// Component for parallel-plate geometries.

class ComponentParallelPlate : public Component {
 public:
  /// Constructor
  ComponentParallelPlate();
  /// Destructor
  ~ComponentParallelPlate() {}

  /// Set the parallel-plate geometries without resistive elements.
  void Setup(double g, double b, double eps, double V);

  /** Set the parallel-plate geometries with resistive elements.
   * \param g size of the gap along positive \f$z\f$.
   * \param b thickness of the resistive layer along negative \f$z\f$.
   * \param eps relative permittivity of the resistive layer.
   * \param sigma conductivity of the resistive layer (must be larger then zero,
   * otherwise do not pass it in the function). 
   * \param V applied potential difference between the parallel plates.
   */
  void Setup(double g, double b, double eps, double V, double sigma);

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

  void DelayedWeightingField(const double x, const double y,
                             const double z, const double t, double& wx,
                             double& wy, double& wz,
                             const std::string& label) override;

  bool GetVoltageRange(double& vmin, double& vmax) override;

  /** Set the parallel-plate geometries with resistive elements.
   * \param x,y position of the center of the electrode in the xy-plane.
   * \param lx_input width in the along \f$x\f$.
   * \param ly_inputwidth in the along \f$y\f$.
   * \param label give name using a string.
   */
  void AddPixel(double x, double y, double lx_input, double ly_input,
                const std::string& label);
  /// Add strip electrode.
  void AddStrip(double x, double lx_input, const std::string& label);

  /// Add plane electrode.
  void AddPlane(const std::string& label);

  // Setting the medium
  void SetMedium(Medium* medium) { m_medium = medium; }

  Medium* GetMedium(const double x, const double y, const double z) override;

 private:
  // Size of the gap.
  double m_g = 0.;
  // Thickness of the resistive element.
  double m_b = 0.;
  // Applied voltage on the electrode to
  // calculate the weighting potential.
  static constexpr double m_Vw = 1.;  
  double m_eps = 1.;
  // Voltage across the RPC
  double m_V = 0.;     
  // Electric field in the gap.
  double m_ezg = 0;
  // Electric field in the resistive layer.
  double m_ezb = 0;  
  double m_sigma = 0;

  std::string m_className = "ComponentParallelPlate";

  Medium* m_medium = nullptr;

  /// Structure that captures the information of the electrodes under study
  struct Electrode {
    std::string type;                      ///< Label.
    int ind = structureelectrode::NotSet;  ///< Readout group.
    double xpos, ypos;                     ///< Coordinates in x/y.
    double lx, ly;                         ///< Dimensions in the x-y plane.
  };

  enum fieldcomponent { xcomp = 0, ycomp, zcomp };

  /// Possible readout groups
  enum structureelectrode { NotSet = -1, Plane, Strip, Pixel };

  /// Vectors storing the readouts
  std::vector<std::string> m_readout;
  std::vector<Electrode> m_readout_p;

  /// Functions that calculate the electric field and potential
  double IntegrateField(const Electrode& el, int comp, const double X,
                        const double Y, const double Z);

  double IntegratePromptPotential(const Electrode& el, const double X,
                                  const double Y, const double Z);
  double IntegrateDelayedPotential(const Electrode& el, const double X,
                                   const double Y, const double Z,
                                   const double t);

  void UpdatePeriodicity() override;
  void Reset() override;
};
}  // namespace Garfield
#endif
