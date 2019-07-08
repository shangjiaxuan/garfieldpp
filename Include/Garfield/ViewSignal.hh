#ifndef G_VIEW_SIGNAL
#define G_VIEW_SIGNAL

#include <memory>
#include <string>

#include <TGraph.h>
#include <TH1D.h>

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
    */
  void PlotSignal(const std::string& label, const bool total = true,
                  const bool electron = false, const bool ion = false, 
                  const bool delayed = false);

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

 private:
  // Sensor
  Sensor* m_sensor = nullptr;

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
  std::unique_ptr<TH1D> m_hDelayedSignal;
  std::unique_ptr<TH1D> m_hDelayedSignalElectrons;
  std::unique_ptr<TH1D> m_hDelayedSignalIons;

  // Threshold crossings
  std::unique_ptr<TGraph> m_gCrossings;
};
}
#endif
