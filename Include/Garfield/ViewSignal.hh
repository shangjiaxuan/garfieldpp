#ifndef G_VIEW_SIGNAL
#define G_VIEW_SIGNAL

#include <memory>
#include <string>

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>

namespace Garfield {

class Sensor;

/// Plot the signal computed by a sensor as a ROOT histogram.

class ViewSignal {
 public:
  /// Constructor
  ViewSignal();
  /// Destructor
  ~ViewSignal();

  /// Set the sensor from which to retrieve the signal.
  void SetSensor(Sensor* s);
  /// Set the pad on which to draw the histogram.
  void SetCanvas(TCanvas* c);

  /** Plot the signal.
    * \param label Identifier (weighting field) of the signal to be plotted.
    * \param total Flag whether to plot the total induced signal.
    * \param electron Flag whether to plot the electron-induced signal.
    * \param ion Flag whether to plot the ion/hole-induced signal.
    */
  void PlotSignal(const std::string& label, const bool total = true,
                  const bool electron = false, const bool ion = false);
  /** Retrieve the histogram for the induced signal.
    * \param h histogram to be returned
               ('t': total, 'e': electron-induced, 'h': ion-induced).
    **/
  TH1D* GetHistogram(const char h = 't') {
    return h == 'e' ? m_hSignalElectrons.get()
                    : 'i' ? m_hSignalIons.get() : m_hSignal.get();
  }

  /// Set the x-axis limits explicitly.
  void SetRangeX(const double xmin, const double xmax);
  /// Remove the user-defined x-axis limits.
  void UnsetRangeX() { m_userRangeX = false; }

  /// Set the y-axis limits explicitly.
  void SetRangeY(const double ymin, const double ymax);
  /// Remove the user-defined y-axis limits.
  void UnsetRangeY() { m_userRangeY = false; }

  /// Enable/disable debugging output.
  void EnableDebugging(const bool on = true) { m_debug = on; }

 private:
  std::string m_className = "ViewSignal";

  // Options
  bool m_debug = false;

  // Sensor
  Sensor* m_sensor = nullptr;

  // Canvas
  TCanvas* m_canvas = nullptr;
  bool m_hasExternalCanvas = false;

  // Axis range.
  double m_xmin = 0.;
  double m_xmax = 0.;
  bool m_userRangeX = false;
  double m_ymin = 0.;
  double m_ymax = 0.;
  bool m_userRangeY = false;

  // Histograms
  std::unique_ptr<TH1D> m_hSignal;
  std::unique_ptr<TH1D> m_hSignalElectrons;
  std::unique_ptr<TH1D> m_hSignalIons;

  // Threshold crossings
  std::unique_ptr<TGraph> m_gCrossings;

  // Find an unused histogram name.
  std::string FindHistogramName(const std::string& base) const;
};
}
#endif
