#ifndef G_VIEW_SIGNAL
#define G_VIEW_SIGNAL

#include <Rtypes.h>
#include <TH1D.h>

#include <array>
#include <memory>
#include <string>

#include "ViewBase.hh"

namespace Garfield {

class Sensor;

/// Plot the signal computed by a sensor as a ROOT histogram.

class ViewSignal : public ViewBase {
 public:
  /// Constructor
  ViewSignal();
  /// Destructor
  ~ViewSignal() = default;

  /// Set the sensor from which to retrieve the signal.
  void SetSensor(Sensor* s);

  /** Plot the signal.
   * \param label Identifier (weighting field) of the signal to be plotted.
   * \param total Flag whether to plot the total induced signal.
   * \param electron Flag whether to plot the electron-induced signal.
   * \param ion Flag whether to plot the ion/hole-induced signal.
   * \param delayed Flag whether to plot the delayed signal component.
   * \param same Flag to keep existing plots on the canvas or not.
   */
  void PlotSignal(const std::string& label, const bool total = true,
                  const bool electron = false, const bool ion = false,
                  const bool delayed = false, const bool same = false);
    
    /** Plot the signal.
     * \param label Identifier (weighting field) of the signal to be plotted.
     * \param setting String containing information about the components you want to plot. The syntax is the first letter of charge carrier signal you want to plot: "t" for total, "e" for electron and "i" for Ion/hole. This is followed by which signal components you want to be plotted: "T" for total, "P" for promt and "D" for delayed. Enabeling all posibilities will look like: "tTPDeTPDiTPD".
     * \param same Flag to keep existing plots on the canvas or not.
     */
    void PlotSignal(const std::string& label, const std::string setting,const bool same = false);

  /** Retrieve the histogram for the total, prompt and delayed induced charge or
   * signal.. \param label Identifier (weighting field) of the signal to be
   * plotted. \param h histogram to be returned \param electron Flag whether to
   * plot the electron-induced signal. \param ion Flag whether to plot the
   * ion/hole-induced signal. \param delayed Flag whether to plot the delayed
   * signal component. \param same Flag to keep existing plots on the canvas or
   * not. \param getsignal is true for plotting the induced signal and false for
   * plotting the induced charge.
   */

  void Plot(const std::string& label, const bool getsignal,
            const bool total = true, const bool delayed = true,
            const bool same = false);

  /** Retrieve the histogram for the total, prompt and delayed induced charge or
    signal.
    * \param h histogram to be returned
               ('t': total, 'e': electron-induced, 'i': ion/hole-induced).
    **/
  
  TH1D* GetHistogram(const char h = 't') {
    return h == 'e' ? m_hSignalElectrons.get()
                    : h == 'i' ? m_hSignalIons.get() : m_hSignal.get();
  }

  /// Set the x-axis limits explicitly.
  void SetRangeX(const double xmin, const double xmax);
  /// Remove the user-defined x-axis limits.
  void UnsetRangeX() { m_userRangeX = false; }

  /// Set the y-axis limits explicitly.
  void SetRangeY(const double ymin, const double ymax);
  /// Remove the user-defined y-axis limits.
  void UnsetRangeY() { m_userRangeY = false; }

  /// Override the default y-axis label.
  void SetLabelY(const std::string& label) { m_labelY = label; }

  /// Set the (ROOT) colour with which to draw the total signal.
  void SetColourTotal(const short col) { m_colTotal = col; }
  /// Set the (ROOT) colour with which to draw the electron component.
  void SetColourElectrons(const short col) { m_colElectrons = col; }
  /// Set the (ROOT) colour with hich to draw the hole/ion component.
  void SetColourIons(const short col) { m_colIons = col; }
  /// Set the (ROOT) colour with hich to draw the hole/ion component.
  void SetColourHoles(const short col) { m_colIons = col; }

  /// Set the (ROOT) colours with which to draw the delayed signal(s).
  void SetColourDelayed(const short colTotal,
                        const short colElectrons = kYellow - 7,
                        const short colIons = kRed - 9) {
    m_colDelayed = {colTotal, colElectrons, colIons};
  }

 private:
  // Sensor.
  Sensor* m_sensor = nullptr;

  // Axis range.
  double m_xmin = 0.;
  double m_xmax = 0.;
  bool m_userRangeX = false;
  double m_ymin = 0.;
  double m_ymax = 0.;
  bool m_userRangeY = false;

  // Axis label.
  std::string m_labelY = "";

  // Histograms.
  std::unique_ptr<TH1D> m_hSignal;
  std::unique_ptr<TH1D> m_hPromptSignal;
    std::unique_ptr<TH1D> m_hPromptElectrons;
    std::unique_ptr<TH1D> m_hPromptIons;

  std::unique_ptr<TH1D> m_hCharge;
  std::unique_ptr<TH1D> m_hPromptCharge;
  std::unique_ptr<TH1D> m_hDelayedCharge;

  std::unique_ptr<TH1D> m_hSignalElectrons;
  std::unique_ptr<TH1D> m_hSignalIons;
  std::unique_ptr<TH1D> m_hDelayedSignal;
  std::unique_ptr<TH1D> m_hDelayedSignalElectrons;
  std::unique_ptr<TH1D> m_hDelayedSignalIons;

  // Colours.
  short m_colTotal = kBlue + 3;
  short m_colElectrons = kOrange - 3;
  short m_colIons = kRed + 1;
  std::array<short, 6> m_colDelayed{{kCyan + 2, kYellow - 7, kRed - 9, kGreen+1, kYellow-4, kRed-9}};
    
    std::array<short, 3> m_colPrompt{{kAzure+10, kRed-4, kMagenta+2}};

  void DrawHistogram(TH1D* h, const std::string& opt,
                     const std::string& xlabel, const std::string& ylabel);
    
  bool ComponentPlotting(std::string entry, std::string prefix);
};
}  // namespace Garfield
#endif
