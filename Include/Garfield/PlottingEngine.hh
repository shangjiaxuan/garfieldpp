#ifndef G_PLOTTING_ENGINE_H
#define G_PLOTTING_ENGINE_H

namespace Garfield {

/// Abstract base class for plotting engines.

class PlottingEngine {
 public:
  /// Default constructor
  PlottingEngine() = delete;
  /// Constructor
  PlottingEngine(const std::string& name) : m_className(name) {}
  /// Destructor
  virtual ~PlottingEngine() {}

  /// Use serif font.
  void SetSerif() { m_serif = true; }
  /// Use sans-serif font.
  void SetSansSerif() { m_serif = false; }

  // Set/get colors.
  void SetLineColor1(const std::string& col) { m_colorLine1 = col; }
  void SetLineColor2(const std::string& col) { m_colorLine2 = col; }
  void SetElectronColor(const std::string& col) { m_colorElectron = col; }
  void SetHoleColor(const std::string& col) { m_colorHole = col; }
  void SetIonColor(const std::string& col) { m_colorIon = col; }
  void SetPhotonColor(const std::string& col) { m_colorPhoton = col; }
  void SetChargedParticleColor(const std::string& col) {
    m_colorChargedParticle = col;
  }

  std::string GetLineColor1() const { return m_colorLine1; }
  std::string GetLineColor2() const { return m_colorLine2; }
  std::string GetElectronColor() const { return m_colorElectron; }
  std::string GetHoleColor() const { return m_colorHole; }
  std::string GetIonColor() const { return m_colorIon; }
  std::string GetPhotonColor() const { return m_colorPhoton; }
  std::string GetChargedParticleColor() const { return m_colorChargedParticle; }

  /// Switch debugging messages on/off.
  void EnableDebugging(const bool on = true) { m_debug = on; }

 protected:
  std::string m_className = "PlottingEngine";

  bool m_serif = false;

  std::string m_colorLine1 = "dark-blue";
  std::string m_colorLine2 = "olive";
  std::string m_colorElectron = "orange";
  std::string m_colorHole = "red";
  std::string m_colorIon = "dark-red";
  std::string m_colorPhoton = "blue";
  std::string m_colorChargedParticle = "dark-green";

  bool m_debug = false;
};
}

#endif
