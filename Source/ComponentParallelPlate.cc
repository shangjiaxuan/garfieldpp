#include "Garfield/ComponentParallelPlate.hh"

#include <TF1.h>
#include <TF2.h>

#include <cmath>
#include <iostream>

#include "Garfield/GarfieldConstants.hh"

namespace Garfield {

ComponentParallelPlate::ComponentParallelPlate() : Component("ParallelPlate") {}

void ComponentParallelPlate::Setup(double g, double b, double eps, double V,
                                   double sigma) {
  std::cout << m_className << "::Setup: Geometry set.\n";
  m_g = g;
  m_b = b;
  m_eps = eps;
  m_V = V;
  m_sigma = sigma;
  // For large times the resistive layer will act as a perfect conductor.
  if (sigma == 0) {
    m_ezg = -m_eps * m_V / (m_b + m_eps * m_g);
    m_ezb = -m_V / (m_b + m_eps * m_g);
    return;
  }

  m_ezg = -m_V / m_g;
  m_ezb = 0.;
}

double ComponentParallelPlate::IntegrateField(const Electrode& el, int comp,
                                              const double x, const double y,
                                              const double z) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      if (comp == fieldcomponent::zcomp)
        return m_eps * m_Vw / (m_b + m_eps * m_g);
      return 0.;
      break;
    }
    case structureelectrode::Pixel: {
      auto WFieldPixel = [=](double* k, double* p) {
        double kx = k[0];
        double ky = k[1];

        delete p;

        double K = std::sqrt(kx * kx + ky * ky);

        double intsol = 1;

        switch (comp) {
          case fieldcomponent::xcomp: {
            intsol *= 1 / (ky * cosh(m_g * K) * sinh(m_b * K) +
                           m_eps * ky * cosh(m_b * K) * sinh(m_g * K));

            intsol *= cos(ky * (y - el.ypos)) * sin((kx * el.lx) / 2) *
                      sin((ky * el.ly) / 2) * sin(kx * (x - el.xpos)) *
                      sinh(K * (m_g - z));
            break;
          }
          case fieldcomponent::ycomp: {
            intsol *= 1 / (kx * cosh(m_g * K) * sinh(m_b * K) +
                           m_eps * kx * cosh(m_b * K) * sinh(m_g * K));
            intsol *= sin(ky * (y - el.ypos)) * sin((kx * el.lx) / 2) *
                      sin((ky * el.ly) / 2) * cos(kx * (x - el.xpos)) *
                      sinh(K * (m_g - z));
            break;
          }
          case fieldcomponent::zcomp: {
            intsol *= 1 / (ky * kx * cosh(m_g * K) * sinh(m_b * K) +
                           m_eps * ky * kx * cosh(m_b * K) * sinh(m_g * K));
            intsol *= K * cos(ky * (y - el.ypos)) * sin((kx * el.lx) / 2) *
                      sin((ky * el.ly) / 2) * cos(kx * (x - el.xpos)) *
                      cosh(K * (m_g - z));
            break;
          }
        }
        return intsol;
      };
      TF2* fw =
          new TF2("WFieldPixel", WFieldPixel, 0, 10 * m_g, 0, 10 * m_g, 0);

      double sol = fw->Integral(0, 10 * m_g, 0, 10 * m_g, 1.e-6);

      delete fw;

      return (4 * m_eps * m_Vw / (Pi * Pi)) * sol;

      break;
    }
    case structureelectrode::Strip: {
      auto WFieldStrip = [=](double* k, double* p) {
        double kk = k[0];

        delete p;

        double intsol = 1 / ((cosh(m_g * kk) * sinh(m_b * kk) +
                              m_eps * cosh(m_b * kk) * sinh(m_g * kk)));
        switch (comp) {
          case fieldcomponent::xcomp: {
            intsol *= (sin(kk * el.lx / 2) * sin(kk * (x - el.xpos)) *
                       sinh(kk * (m_g - z)));
            break;
          }
          case fieldcomponent::zcomp: {
            intsol *= (sin(kk * el.lx / 2) * cos(kk * (x - el.xpos)) *
                       cosh(kk * (m_g - z)));
            break;
          }
        }

        return intsol;
      };

      if (comp == ycomp) {
        return 0.;
      } else {
        TF1* fw = new TF1("WFieldStrip", WFieldStrip, 0, 10 * m_g, 0);

        double sol = fw->Integral(0, 10 * m_g);

        delete fw;

        return (2 * m_eps * m_Vw / (Pi)) * sol;
      }
      break;
    }
    default: {
      std::cout
          << m_className
          << ":: IntegrateField: The structure of the electrode is not set."
          << "\n"
          << std::endl;
      return 0.;
    }
  }
}

double ComponentParallelPlate::IntegrateDelayedField(const Electrode& el,
                                                     int comp, const double x,
                                                     const double y,
                                                     const double z,
                                                     const double t) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      if (comp == fieldcomponent::zcomp)
        return m_eps * m_Vw *
               (1 - exp(-t * m_g * m_sigma / (m_eps0 * (m_b + m_eps * m_g)))) /
               (m_b + m_eps * m_g);
      return 0.;
      break;
    }
    case structureelectrode::Pixel: {
      auto WFieldPixel = [=](double* k, double* p) {
        double kx = k[0];
        double ky = k[1];

        delete p;

        double K = std::sqrt(kx * kx + ky * ky);

        double tau = m_eps0 *
                     (m_eps + cosh(m_g * K) * sinh(m_b * K) /
                                  (cosh(m_b * K) * sinh(m_g * K))) *
                     (1 / m_sigma);

        double intsol = 1 / (cosh(m_g * K) * sinh(m_b * K) +
                             m_eps * cosh(m_b * K) * sinh(m_g * K));

        switch (comp) {
          case fieldcomponent::xcomp: {
            intsol *= (1 - exp(-t / tau)) * cos(ky * (y - el.ypos)) *
                      cosh(m_g * K) * sin((kx * el.lx) / 2) *
                      sin(kx * (x - el.xpos)) * sinh(K * (m_g - z)) *
                      tanh(m_b * K) / (ky * sinh(m_g * K));
            break;
          }
          case fieldcomponent::ycomp: {
            intsol *= (1 - exp(-t / tau)) * sin(ky * (y - el.ypos)) *
                      cosh(m_g * K) * sin((kx * el.lx) / 2) *
                      cos(kx * (x - el.xpos)) * cosh(K * (m_g - z)) *
                      tanh(m_b * K) / (kx * sinh(m_g * K));
            break;
          }
          case fieldcomponent::zcomp: {
            intsol *= (1 - exp(-t / tau)) * cos(ky * (y - el.ypos)) *
                      cosh(m_g * K) * sin((K * el.lx) / 2) *
                      cos(kx * (x - el.xpos)) * cosh(K * (m_g - z)) *
                      tanh(m_b * K) / (kx * ky * sinh(m_g * K));
            break;
          }
        }
        return intsol;
      };

      TF2* fw =
          new TF2("WFieldPixel", WFieldPixel, 0, 10 * m_g, 0, 10 * m_g, 0);

      double sol = fw->Integral(0, 10 * m_g, 0, 10 * m_g, 1.e-6);

      delete fw;

      return (4 * m_eps * m_Vw / (Pi * Pi)) * sol;

      break;
    }
    case structureelectrode::Strip: {
      auto WFieldStrip = [=](double* k, double* p) {
        double kk = k[0];
        delete p;

        double tau = m_eps0 *
                     (m_eps + cosh(m_g * kk) * sinh(m_b * kk) /
                                  (cosh(m_b * kk) * sinh(m_g * kk))) *
                     (1 / m_sigma);

        double intsol = 1 / (cosh(m_g * kk) * sinh(m_b * kk) +
                             m_eps * cosh(m_b * kk) * sinh(m_g * kk));
        switch (comp) {
          case fieldcomponent::xcomp: {
            intsol *= (1 - exp(-t / tau)) * cosh(m_g * kk) *
                      sin((kk * el.lx) / 2) * sin(kk * (x - el.xpos)) *
                      sinh(kk * (m_g - z)) * tanh(m_b * kk) / sinh(m_g * kk);
            break;
          }
          case fieldcomponent::zcomp: {
            intsol *= (1 - exp(-t / tau)) * cosh(m_g * kk) *
                      sin((kk * el.lx) / 2) * cos(kk * (x - el.xpos)) *
                      cosh(kk * (m_g - z)) * tanh(m_b * kk) / sinh(m_g * kk);
            break;
          }
        }

        return intsol;
      };

      if (comp == ycomp) {
        return 0.;
      } else {
        TF1* fw = new TF1("WFieldStrip", WFieldStrip, 0, 10 * m_g, 0);

        double sol = fw->Integral(0, 10 * m_g);

        delete fw;

        return (2 * m_eps * m_Vw / Pi) * sol;
      }
      break;
    }
    default: {
      std::cout
          << m_className
          << ":: IntegrateField: The structure of the electrode is not set."
          << "\n"
          << std::endl;
      return 0.;
    }
  }
}

double ComponentParallelPlate::IntegratePromptPotential(const Electrode& el,
                                                        const double x,
                                                        const double y,
                                                        const double z) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      double sol = m_eps * m_Vw * (m_g - z) / (m_b + m_eps * m_g);
      return abs(sol) > m_precision ? sol : 0.;
      break;
    }
    case structureelectrode::Pixel: {
      auto WPotentialPixel = [=](double* k, double* p) {
        double kx = k[0];
        double ky = k[1];

        delete p;

        double K = std::sqrt(kx * kx + ky * ky);

        double intsol = 1;

        intsol *= cos(kx * (x - el.xpos)) * sin(kx * el.lx / 2) *
                  cos(ky * (y - el.ypos)) * sin(ky * el.ly / 2) *
                  sinh(K * (m_g - z)) /
                  (kx * ky *
                   (sinh(m_b * K) * cosh(m_g * K) +
                    m_eps * sinh(m_g * K) * cosh(m_b * K)));

        return intsol;
      };

      TF2* pw = new TF2("WPotentialPixel", WPotentialPixel, 0, 10 * m_g, 0,
                        10 * m_g, 0);

      double sol = pw->Integral(0, 2 * m_g, 0, 2 * m_g, 1.e-6);

      delete pw;

      return (4 * m_eps * m_Vw / (Pi * Pi)) * sol;
      break;
    }
    case structureelectrode::Strip: {
      auto WPotentialStrip = [=](double* k, double* p) {
        double kk = k[0];

        delete p;

        double intsol = 1 / (kk * (cosh(m_g * kk) * sinh(m_b * kk) +
                                   m_eps * cosh(m_b * kk) * sinh(m_g * kk)));
        intsol *= (sin(kk * el.lx / 2) * cos(kk * (x - el.xpos)) *
                   sinh(kk * (m_g - z)));

        return intsol;
      };

      TF1* pw = new TF1("WPotentialStrip", WPotentialStrip, 0, 10 * m_g, 0);

      double sol = pw->Integral(0, 10 * m_g);
      delete pw;
      return (2 * m_eps * m_Vw / (Pi)) * sol;
      break;
    }
    default: {
      std::cout << m_className
                << ":: IntegratePromptPotential: The structure of the "
                   "electrode is not set."
                << "\n"
                << std::endl;
      return 0.;
    }
  }
}

double ComponentParallelPlate::IntegrateDelayedPotential(const Electrode& el,
                                                         const double x,
                                                         const double y,
                                                         const double z,
                                                         const double t) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      double tau =
          m_eps0 * (m_eps + m_b / m_g) /
          (m_sigma);  // Note to self: You dropt the eps) here for convenience.

      double sol = m_Vw * (1 - exp(-t / tau)) * (m_b * (m_g - z)) /
                   (m_g * (m_b + m_eps * m_g));
      return abs(sol) > m_precision ? sol : 0.;
      break;
    }
    case structureelectrode::Pixel: {
      auto WPotentialPixel = [=](double* k, double* p) {
        double kx = k[0];
        double ky = k[1];

        delete p;

        double K = std::sqrt(kx * kx + ky * ky);
        double tau = m_eps0 *
                     (m_eps + cosh(m_g * K) * sinh(m_b * K) /
                                  (cosh(m_b * K) * sinh(m_g * K))) *
                     (1 / m_sigma);  // Note to self: You dropt the eps) here
                                     // for convenience.

        double intsol = 1 / (kx * ky *
                             (sinh(m_b * K) * cosh(m_g * K) +
                              m_eps * sinh(m_g * K) * cosh(m_b * K)));

        intsol *= cos(kx * (x - el.xpos)) * sin(kx * el.lx / 2) *
                  cos(ky * (y - el.ypos)) * sin(ky * el.ly / 2) *
                  sinh(K * (m_g - z)) * tanh(m_b * K) * cosh(m_g * K) *
                  (1 - exp(-t / tau)) / sinh(m_g * K);

        return intsol;
      };

      TF2* pw = new TF2("WPotentialPixel", WPotentialPixel, 0, 10 * m_g, 0,
                        10 * m_g, 0);

      double sol = pw->Integral(0, 2 * m_g, 0, 2 * m_g, 1.e-20);

      delete pw;

      return (4 * m_Vw / (Pi * Pi)) * sol;
      break;
    }
    case structureelectrode::Strip: {
      auto WPotentialStrip = [=](double* k, double* p) {
        double kk = k[0];

        delete p;

        double tau = m_eps0 *
                     (m_eps + cosh(m_g * kk) * sinh(m_b * kk) /
                                  (cosh(m_b * kk) * sinh(m_g * kk))) *
                     (1 / m_sigma);

        double intsol = 1 / (kk * (cosh(m_g * kk) * sinh(m_b * kk) +
                                   m_eps * cosh(m_b * kk) * sinh(m_g * kk)));
        intsol *= (sin(kk * el.lx / 2) * cos(kk * (x - el.xpos)) *
                   sinh(kk * (m_g - z)) * cosh(m_g * kk) * tanh(m_b * kk)) *
                  (1 - exp(-t / tau)) / sinh(m_g * kk);

        return intsol;
      };

      TF1* pw = new TF1("WPotentialStrip", WPotentialStrip, 0, 10 * m_g, 0);

      double sol = pw->Integral(0, 8 * m_g);
      delete pw;
      return (2 * m_Vw / (Pi)) * sol;
      break;
    }
    default: {
      std::cout << m_className
                << ":: IntegrateDelayedPotential: The structure of the "
                   "electrode is not set."
                << "\n"
                << std::endl;
      return 0.;
    }
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double& ex,
                                           double& ey, double& ez, Medium*& m,
                                           int& status) {
  ex = ey = 0.;

  if (z < 0) {
    ez = m_ezb;
  } else {
    ez = m_ezb;
  }

  m = m_geometry ? m_geometry->GetMedium(x, y, z) : m_medium;

  if (!m) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField: No medium at (" << x << ", "
                << y << ", " << z << ").\n";
    }
    status = -6;
    return;
  }

  if (z > 0) {
    status = 0;
  } else {
    status = -5;
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double& ex,
                                           double& ey, double& ez, double& v,
                                           Medium*& m, int& status) {
  ex = ey = 0.;

  if (z > 0.) {
    ez = m_ezg;
  } else {
    ez = m_ezb;
  }

  if (m_sigma == 0) {
    v = -m_eps * m_V * (m_g - z) / (m_b + m_eps * m_g);
  } else {
    v = -m_eps * m_V * (m_g - z) / (m_eps * m_g);
  }

  m = m_geometry ? m_geometry->GetMedium(x, y, z) : m_medium;
  if (!m) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField: No medium at (" << x << ", "
                << y << ", " << z << ").\n";
    }
    status = -6;
    return;
  }

  if (z > 0) {
    status = 0;
  } else {
    status = -5;
  }
}

bool ComponentParallelPlate::GetVoltageRange(double& vmin, double& vmax) {
  if (m_V == 0) return false;

  if (m_V < 0) {
    vmin = m_V;
    vmax = 0;
  } else {
    vmin = 0;
    vmax = m_V;
  }
  return true;
}

void ComponentParallelPlate::WeightingField(const double x, const double y,
                                            const double z, double& wx,
                                            double& wy, double& wz,
                                            const std::string& label) {
  wx = 0;
  wy = 0;
  wz = 0;

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      wx = electrode.flip *
           IntegrateField(electrode, fieldcomponent::xcomp, x, y, z);
      wy = electrode.flip *
           IntegrateField(electrode, fieldcomponent::ycomp, x, y, z);
      wz = electrode.flip *
           IntegrateField(electrode, fieldcomponent::zcomp, x, y, z);
    }
  }
}

double ComponentParallelPlate::WeightingPotential(const double x,
                                                  const double y,
                                                  const double z,
                                                  const std::string& label) {
  double ret = 0.;

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      if (!electrode.m_usegrid) {
        ret += electrode.flip * IntegratePromptPotential(electrode, x, y, z);
      } else {
        ret += FindWeightingPotentialInGrid(electrode, x, y, z);
      }
    }
  }
  return ret;
}

double ComponentParallelPlate::DelayedWeightingPotential(
    const double x, const double y, const double z, const double t,
    const std::string& label) {
  if (m_sigma == 0) {
    std::cout << m_className
              << ":: DelayedWeightingPotential: No conductivity set."
              << "\n"
              << std::endl;
    return 0.;
  }

  double ret = 0.;

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      if (!electrode.m_usegrid) {
        ret +=
            electrode.flip * IntegrateDelayedPotential(electrode, x, y, z, t);
      } else {
        ret += FindDelayedWeightingPotentialInGrid(electrode, x, y, z, t);
      }
    }
  }

  return ret;
}

void ComponentParallelPlate::DelayedWeightingField(
    const double x, const double y, const double z, const double t, double& wx,
    double& wy, double& wz, const std::string& label) {
  wx = 0.;
  wy = 0.;
  wz = 0.;

  if (m_sigma == 0) {
    std::cout << m_className << ":: DelayedWeightingField: No conductivity set."
              << "\n"
              << std::endl;
  }

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      wx = electrode.flip *
           IntegrateDelayedField(electrode, fieldcomponent::xcomp, x, y, z, t);
      wy = electrode.flip *
           IntegrateDelayedField(electrode, fieldcomponent::ycomp, x, y, z, t);
      wz = electrode.flip *
           IntegrateDelayedField(electrode, fieldcomponent::zcomp, x, y, z, t);
    }
  }
}

void ComponentParallelPlate::Reset() {
  m_readout.clear();
  m_readout_p.clear();

  m_g = 0.;
  m_b = 0.;
  m_eps = 1.;
  m_V = 0.;
}

void ComponentParallelPlate::UpdatePeriodicity() {
  if (m_debug) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Periodicities are not supported.\n";
  }
}

void ComponentParallelPlate::AddPixel(double x, double y, double lx_input,
                                      double ly_input,
                                      const std::string& label) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddPixel:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode pixel;
  pixel.label = label;
  pixel.ind = structureelectrode::Pixel;
  pixel.xpos = x;
  pixel.ypos = y;
  pixel.lx = lx_input;
  pixel.ly = ly_input;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(pixel));
  std::cerr << m_className << "::AddPixel:\n"
            << "The pixel electrode structure has been added.\n";
}

void ComponentParallelPlate::AddStrip(double x, double lx_input,
                                      const std::string& label) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddStrip:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode strip;
  strip.label = label;
  strip.ind = structureelectrode::Strip;
  strip.xpos = x;
  strip.lx = lx_input;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(strip));

  std::cerr << m_className << "::AddStrip:\n"
            << "The strip electrode structure has been added.\n";
}

void ComponentParallelPlate::AddPlane(const std::string& label, bool anode) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddPlane:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode plate;
  plate.label = label;
  plate.ind = structureelectrode::Plane;

  if (!anode) plate.flip = -1;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(plate));

  std::cerr << m_className << "::AddPlane:\n"
            << "The plane electrode structure has been added.\n";
}

Medium* ComponentParallelPlate::GetMedium(const double x, const double y,
                                          const double z) {
  if (m_geometry) {
    return m_geometry->GetMedium(x, y, z);
  } else if (m_medium) {
    return m_medium;
  }
  return nullptr;
}

void ComponentParallelPlate::SetWeightingPotentialGrid(
    const std::string& label, const double xmin, const double xmax,
    const double xsteps, const double ymin, const double ymax,
    const double ysteps, const double zmin, const double zmax,
    const double zsteps, const double tmin, const double tmax,
    const double tsteps) {
  for (auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      electrode.gridXSteps = xsteps;
      electrode.gridYSteps = ysteps;
      electrode.gridZSteps = zsteps;
      electrode.gridTSteps = tsteps;

      if (xsteps == 0) electrode.gridXSteps = 1;
      if (ysteps == 0) electrode.gridYSteps = 1;

      electrode.gridX0 = xmin;
      electrode.gridY0 = ymin;
      electrode.gridZ0 = zmin;
      electrode.gridT0 = tmin;

      electrode.gridXStepSize = (xmax - xmin) / xsteps;
      electrode.gridYStepSize = (ymax - ymin) / ysteps;
      electrode.gridZStepSize = (zmax - zmin) / zsteps;
      electrode.gridTStepSize = (tmax - tmin) / tsteps;

      std::vector<double> nhz(zsteps, 0);
      std::vector<std::vector<double>> nhy(ysteps, nhz);
      std::vector<std::vector<std::vector<double>>> nhx(xsteps, nhy);
      electrode.gridPromptV = nhx;

      std::vector<double> nht(tsteps, 0);
      std::vector<std::vector<double>> nhzd(zsteps, nht);
      std::vector<std::vector<std::vector<double>>> nhyd(ysteps, nhzd);
      std::vector<std::vector<std::vector<std::vector<double>>>> nhxd(xsteps,
                                                                      nhyd);
      electrode.gridDelayedV = nhxd;

      for (int ix = 0; ix < xsteps; ix++) {
        for (int iy = 0; iy < xsteps; iy++) {
          for (int iz = 0; iz < xsteps; iz++) {
            if (iz * zsteps + zmin >= 0)
              electrode.gridPromptV[ix][iy][iz] =
                  electrode.flip * IntegratePromptPotential(
                                       electrode, ix * xsteps + xmin,
                                       iy * ysteps + ymin, iz * zsteps + zmin);

            for (int it = 0; it < tsteps; it++) {
              if (iz * zsteps + zmin >= 0)
                electrode.gridDelayedV[ix][iy][iz][it] =
                    electrode.flip * IntegrateDelayedPotential(
                                         electrode, ix * xsteps + xmin,
                                         iy * ysteps + ymin, iz * zsteps + zmin,
                                         it * tsteps + tmin);
            }
          }
        }
      }

      electrode.m_usegrid = true;
    }
  }
}

void ComponentParallelPlate::SetWeightingPotentialGrids(
    const double xmin, const double xmax, const double xsteps,
    const double ymin, const double ymax, const double ysteps,
    const double zmin, const double zmax, const double zsteps,
    const double tmin, const double tmax, const double tsteps) {
  for (const auto& electrode : m_readout_p) {
    SetWeightingPotentialGrid(electrode.label, xmin, xmax, xsteps, ymin, ymax,
                              ysteps, zmin, zmax, zsteps, tmin, tmax, tsteps);
  }
}

double ComponentParallelPlate::FindWeightingPotentialInGrid(const Electrode& el,
                                                            const double x,
                                                            const double y,
                                                            const double z) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      return el.flip * IntegratePromptPotential(el, x, y, z);
      break;
    }
    case structureelectrode::Strip: {
      int ix = floor((x - el.gridX0) / el.gridXStepSize);
      int iz = floor((z - el.gridZ0) / el.gridZStepSize);

      if (ix < 0 || ix >= el.gridXSteps || iz < 0 || iz >= el.gridZSteps)
        return IntegratePromptPotential(el, x, y, z);

      double ret = 0;

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          ret += FindWeightFactor(
                     el, abs((ix + i) * el.gridXStepSize + el.gridX0 - x), 0,
                     abs((iz + j) * el.gridZStepSize + el.gridZ0 - z)) *
                 el.gridPromptV[ix + i][0][iz + j];
        }
      }

      return ret;
      break;
    }
    case structureelectrode::Pixel: {
      int ix = floor((x - el.gridX0) / el.gridXStepSize);
      int iy = floor((y - el.gridY0) / el.gridYStepSize);
      int iz = floor((z - el.gridZ0) / el.gridZStepSize);

      if (ix < 0 || ix >= el.gridXSteps || iz < 0 || iz >= el.gridYSteps ||
          iz < 0 || iz >= el.gridZSteps)
        return IntegratePromptPotential(el, x, y, z);

      double ret = 0;

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          for (int k = 0; k < 2; k++) {
            ret += FindWeightFactor(
                       el, abs((ix + i) * el.gridXStepSize + el.gridX0 - x),
                       abs((iy + k) * el.gridYStepSize + el.gridY0 - y),
                       abs((iz + j) * el.gridZStepSize + el.gridZ0 - z)) *
                   el.gridPromptV[ix + i][iy + k][iz + j];
          }
        }
      }
      return ret;
      break;
    }
  }
  return 0.;
}

double ComponentParallelPlate::FindDelayedWeightingPotentialInGrid(
    const Electrode& el, const double x, const double y, const double z,
    const double t) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      return el.flip * IntegrateDelayedPotential(el, x, y, z, t);
      break;
    }
    case structureelectrode::Strip: {
      int ix = floor((x - el.gridX0) / el.gridXStepSize);
      int iz = floor((z - el.gridZ0) / el.gridZStepSize);
      int it = floor((t - el.gridT0) / el.gridTStepSize);

      if (ix < 0 || ix >= el.gridXSteps || iz < 0 || iz >= el.gridZSteps ||
          it < 0 || it >= el.gridTSteps)
        return IntegrateDelayedPotential(el, x, y, z, t);

      double ret = 0;

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          for (int l = 0; l < 2; l++) {
            ret += FindWeightFactor(
                       el, abs((ix + i) * el.gridXStepSize + el.gridX0 - x), 0,
                       abs((iz + j) * el.gridZStepSize + el.gridZ0 - z),
                       abs((it + l) * el.gridTStepSize + el.gridT0 - t)) *
                   el.gridDelayedV[ix + i][0][iz + j][it + l];
          }
        }
      }

      return ret;
      break;
    }
    case structureelectrode::Pixel: {
      int ix = floor((x - el.gridX0) / el.gridXStepSize);
      int iy = floor((y - el.gridY0) / el.gridYStepSize);
      int iz = floor((z - el.gridZ0) / el.gridZStepSize);
      int it = floor((t - el.gridT0) / el.gridTStepSize);

      if (ix < 0 || ix >= el.gridXSteps || iz < 0 || iz >= el.gridYSteps ||
          iz < 0 || iz >= el.gridZSteps || it < 0 || it >= el.gridTSteps)
        return IntegrateDelayedPotential(el, x, y, z, t);

      double ret = 0;

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          for (int k = 0; k < 2; k++) {
            for (int l = 0; l < 2; l++) {
              ret += FindWeightFactor(
                         el, abs((ix + i) * el.gridXStepSize + el.gridX0 - x),
                         abs((iy + k) * el.gridYStepSize + el.gridY0 - y),
                         abs((iz + j) * el.gridZStepSize + el.gridZ0 - z),
                         abs((it + l) * el.gridTStepSize + el.gridT0 - t)) *
                     el.gridDelayedV[ix + i][iy + k][iz + j][it + l];
            }
          }
        }
      }
      return ret;
      break;
    }
  }
  return 0.;
}

double ComponentParallelPlate::FindWeightFactor(const Electrode& el,
                                                const double dx,
                                                const double dy,
                                                const double dz) {
  double fact = 0;

  switch (el.ind) {
    case structureelectrode::Strip: {
      fact = (el.gridXStepSize - dx) * (el.gridZStepSize - dz) /
             (el.gridXStepSize * el.gridZStepSize);
      break;
    }
    case structureelectrode::Pixel: {
      fact = (el.gridXStepSize - dx) * (el.gridYStepSize - dy) *
             (el.gridZStepSize - dz) /
             (el.gridXStepSize * el.gridYStepSize * el.gridZStepSize);
      break;
    }
  }

  return fact;
}

double ComponentParallelPlate::FindWeightFactor(const Electrode& el,
                                                const double dx,
                                                const double dy,
                                                const double dz,
                                                const double dt) {
  double fact = 0;

  switch (el.ind) {
    case structureelectrode::Strip: {
      fact = (el.gridXStepSize - dx) * (el.gridZStepSize - dz) *
             (el.gridXStepSize - dt) /
             (el.gridXStepSize * el.gridZStepSize * el.gridTStepSize);
      break;
    }
    case structureelectrode::Pixel: {
      fact = (el.gridXStepSize - dx) * (el.gridYStepSize - dy) *
             (el.gridZStepSize - dz) * (el.gridXStepSize - dt) /
             (el.gridXStepSize * el.gridYStepSize * el.gridZStepSize *
              el.gridTStepSize);
      break;
    }
  }

  return fact;
}

}  // namespace Garfield
