#include "Garfield/ComponentParallelPlate.hh"

#include <TF1.h>
#include <TF2.h>

#include <cmath>
#include <limits>
#include <iostream>

#include "Garfield/GarfieldConstants.hh"

namespace Garfield {

ComponentParallelPlate::ComponentParallelPlate() : Component("ParallelPlate") {}

void ComponentParallelPlate::Setup(double g, double b, double eps, double V,
                                   double sigma) {
    
  // Here I switch conventions with the z-axis the direction of drift.
    if(g<=0||b<=0){
        std::cerr << m_className << "::Setup: Parameters b and g must be larger than zero.\n";
        return;
    }
  m_g = g;
  m_b = b;
  if (eps < 1.) {
    std::cerr << m_className << "::Setup: Epsilon must be >= 1.\n";
    return;
  }
  m_eps = eps;
  m_V = V;
  // TODO: can sigma be negative?
  m_sigma = sigma;
  if (sigma == 0) {
    m_ezg = -m_eps * m_V / (m_b + m_eps * m_g);
    m_ezb = -m_V / (m_b + m_eps * m_g);
  } else if(sigma < 0){
      std::cerr << m_className << "::Setup: Parameter sigma must be larger than zero.\n";
      return;
  }else {
    // For large times the resistive layer will act as a perfect conductor.
    m_ezg = -m_V / m_g;
    m_ezb = 0.;
  }
  std::cout << m_className << "::Setup: Geometry set.\n";
}

bool ComponentParallelPlate::GetBoundingBox(double& x0, double& y0, double& z0,
                                            double& x1, double& y1,
                                            double& z1) {
  // If a geometry is present, try to get the bounding box from there.
    
    // Here I switch conventions back with the y-axis the direction of drift.
  if (m_geometry) {
    if (m_geometry->GetBoundingBox(x0, y0, z0, x1, y1, z1)) return true;
  }
  z0 = -std::numeric_limits<double>::infinity();
  x0 = -std::numeric_limits<double>::infinity();
  z1 = +std::numeric_limits<double>::infinity();
  x1 = +std::numeric_limits<double>::infinity();
  // TODO: check!
  y0 = 0.;
  y1 = m_g;
  return true;
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
      auto WFieldPixel = [=](double* k, double* /*p*/) {
        double kx = k[0];
        double ky = k[1];

        double K = std::sqrt(kx * kx + ky * ky);

        double intsol = 1.;

        switch (comp) {
          case fieldcomponent::xcomp: {
            intsol *= 1. / (ky * cosh(m_g * K) * sinh(m_b * K) +
                           m_eps * ky * cosh(m_b * K) * sinh(m_g * K));

            intsol *= cos(ky * (y - el.ypos)) * sin((kx * el.lx) / 2) *
                      sin((ky * el.ly) / 2) * sin(kx * (x - el.xpos)) *
                      sinh(K * (m_g - z));
            break;
          }
          case fieldcomponent::ycomp: {
            intsol *= 1. / (kx * cosh(m_g * K) * sinh(m_b * K) +
                           m_eps * kx * cosh(m_b * K) * sinh(m_g * K));
            intsol *= sin(ky * (y - el.ypos)) * sin((kx * el.lx) / 2) *
                      sin((ky * el.ly) / 2) * cos(kx * (x - el.xpos)) *
                      sinh(K * (m_g - z));
            break;
          }
          case fieldcomponent::zcomp: {
            intsol *= 1. / (ky * kx * cosh(m_g * K) * sinh(m_b * K) +
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
      const double sol = fw->Integral(0, 10 * m_g, 0, 10 * m_g, 1.e-6);
      delete fw;
      return (4 * m_eps * m_Vw / Pi2) * sol;

      break;
    }
    case structureelectrode::Strip: {
      auto WFieldStrip = [=](double* k, double* /*p*/) {
        double kk = k[0];

        double intsol = 1. / ((cosh(m_g * kk) * sinh(m_b * kk) +
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
      }
      TF1* fw = new TF1("WFieldStrip", WFieldStrip, 0, 10 * m_g, 0);
      double sol = fw->Integral(0, 10 * m_g);
      delete fw;
      return (2 * m_eps * m_Vw / Pi) * sol;
      break;
    }
    default: {
      std::cerr << m_className << "::IntegrateField: Unknown electrode type.\n";
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
      auto WFieldPixel = [=](double* k, double* /*p*/) {
        double kx = k[0];
        double ky = k[1];

        double K = std::sqrt(kx * kx + ky * ky);

        double tau = m_eps0 *
                     (m_eps + cosh(m_g * K) * sinh(m_b * K) /
                                  (cosh(m_b * K) * sinh(m_g * K))) *
                     (1 / m_sigma);

        double intsol = 1. / (cosh(m_g * K) * sinh(m_b * K) +
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
      const double sol = fw->Integral(0, 10 * m_g, 0, 10 * m_g, 1.e-6);
      delete fw;
      return (4 * m_eps * m_Vw / Pi2) * sol;
      break;
    }
    case structureelectrode::Strip: {
      auto WFieldStrip = [=](double* k, double* /*p*/) {
        double kk = k[0];
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
      } 
      TF1* fw = new TF1("WFieldStrip", WFieldStrip, 0, 10 * m_g, 0);
      const double sol = fw->Integral(0, 10 * m_g);
      delete fw;
      return (2 * m_eps * m_Vw / Pi) * sol;
      break;
    }
    default: {
      std::cerr << m_className << "::IntegrateDelayedField:\n"
                << "    Unknown electrode type.\n";
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
      return std::abs(sol) > m_precision ? sol : 0.;
      break;
    }
    case structureelectrode::Pixel: {
      auto WPotentialPixel = [=](double* k, double* /*p*/) {
        double kx = k[0];
        double ky = k[1];

        double K = std::sqrt(kx * kx + ky * ky);

        double intsol = 1.;

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
      const double sol = pw->Integral(0, 2 * m_g, 0, 2 * m_g, 1.e-6);
      delete pw;
      return (4 * m_eps * m_Vw / Pi2) * sol;
      break;
    }
    case structureelectrode::Strip: {
      auto WPotentialStrip = [=](double* k, double* /*p*/) {
        double kk = k[0];

        double intsol = 1. / (kk * (cosh(m_g * kk) * sinh(m_b * kk) +
                                   m_eps * cosh(m_b * kk) * sinh(m_g * kk)));
        intsol *= (sin(kk * el.lx / 2) * cos(kk * (x - el.xpos)) *
                   sinh(kk * (m_g - z)));

        return intsol;
      };

      TF1* pw = new TF1("WPotentialStrip", WPotentialStrip, 0, 10 * m_g, 0);
      const double sol = pw->Integral(0, 10 * m_g);
      delete pw;
      return (2 * m_eps * m_Vw / Pi) * sol;
      break;
    }
    default: {
      std::cerr << m_className << "::IntegratePromptPotential:\n"
                << "    Unknown electrode type.\n";
      return 0.;
    }
  }
}

double ComponentParallelPlate::IntegrateDelayedPotential(const Electrode& el,
                                                         const double x,
                                                         const double y,
                                                         const double z,
                                                         const double t) {
    
    // TODO: Find better integration methode!
    
  switch (el.ind) {
    case structureelectrode::Plane: {
      double tau =
          m_eps0 * (m_eps + m_b / m_g) /
          (m_sigma);  // Note to self: You dropt the eps) here for convenience.

      double sol = m_Vw * (1 - exp(-t / tau)) * (m_b * (m_g - z)) /
                   (m_g * (m_b + m_eps * m_g));
      return std::abs(sol) > m_precision ? sol : 0.;
      break;
    }
    case structureelectrode::Pixel: {
      auto WPotentialPixel = [=](double* k, double* /*p*/) {
        double kx = k[0];
        double ky = k[1];

        double K = std::sqrt(kx * kx + ky * ky);
        double tau = m_eps0 *
                     (m_eps + cosh(m_g * K) * sinh(m_b * K) /
                                  (cosh(m_b * K) * sinh(m_g * K))) *
                     (1 / m_sigma);  // Note to self: You dropt the eps) here
                                     // for convenience.

        double intsol = 1. / (kx * ky *
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
      const double sol = pw->Integral(0, 2 * m_g, 0, 2 * m_g, 1.e-20);
      delete pw;
      return (4 * m_Vw / Pi2) * sol;
      break;
    }
    case structureelectrode::Strip: {
      auto WPotentialStrip = [=](double* k, double* /*p*/) {
        double kk = k[0];

        double tau = m_eps0 *
                     (m_eps + cosh(m_g * kk) * sinh(m_b * kk) /
                                  (cosh(m_b * kk) * sinh(m_g * kk))) *
                     (1 / m_sigma);

        double intsol = 1. / (kk * (cosh(m_g * kk) * sinh(m_b * kk) +
                                   m_eps * cosh(m_b * kk) * sinh(m_g * kk)));
        intsol *= (sin(kk * el.lx / 2) * cos(kk * (x - el.xpos)) *
                   sinh(kk * (m_g - z)) * cosh(m_g * kk) * tanh(m_b * kk)) *
                  (1 - exp(-t / tau)) / sinh(m_g * kk);

        return intsol;
      };

      TF1* pw = new TF1("WPotentialStrip", WPotentialStrip, 0, 10 * m_g, 0);
      const double sol = pw->Integral(0, 8 * m_g);
      delete pw;
      return (2 * m_Vw / Pi) * sol;
      break;
    }
    default: {
      std::cerr << m_className << "::IntegrateDelayedPotential:\n"
                << "    Unknown electrode type.\n";
      return 0.;
    }
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double& ex,
                                           double& ey, double& ez, Medium*& m,
                                           int& status) {
    // Here I switch conventions back with the y-axis the direction of drift.
    
  ex = ez = 0.;

  if (y < 0) {
    ey = m_ezb;
  } else {
    ey = m_ezb;
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

  if (y > 0) {
    status = 0;
  } else {
    status = -5;
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double& ex,
                                           double& ey, double& ez, double& v,
                                           Medium*& m, int& status) {
    // Here I switch conventions back with the y-axis the direction of drift.
    
  ex = ez = 0.;

  if (y > 0.) {
    ey = m_ezg;
  } else {
    ey = m_ezb;
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

  if (y > 0) {
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
    
    // Here I switch conventions back with the y-axis the direction of drift.
    
  wx = 0;
  wy = 0;
  wz = 0;

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      wx = electrode.flip *
           IntegrateField(electrode, fieldcomponent::ycomp, z, x, y);
      wy = electrode.flip *
           IntegrateField(electrode, fieldcomponent::zcomp, z, x, y);
      wz = electrode.flip *
           IntegrateField(electrode, fieldcomponent::xcomp, z, x, y);
    }
  }
}

double ComponentParallelPlate::WeightingPotential(const double x,
                                                  const double y,
                                                  const double z,
                                                  const std::string& label) {
    
    // Here I switch conventions back with the y-axis the direction of drift.
    
  double ret = 0.;

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      if (!electrode.m_usegrid) {
        ret += electrode.flip * IntegratePromptPotential(electrode, z, x, y);
      } else {
        ret += FindWeightingPotentialInGrid(electrode, z, x, y);
      }
    }
  }
  return ret;
}

double ComponentParallelPlate::DelayedWeightingPotential(
    const double x, const double y, const double z, const double t,
    const std::string& label) {
    
    // Here I switch conventions back with the y-axis the direction of drift.
    
  if (m_sigma == 0) {
    if (m_debug) {
      std::cout << m_className << "::DelayedWeightingPotential:\n"
                << "    Conductivity is set to zero.\n";
    }
    return 0.;
  }

  double ret = 0.;

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      if (!electrode.m_usegrid) {
        ret +=
            electrode.flip * IntegrateDelayedPotential(electrode, z, x, y,t);
      } else {
        ret += FindDelayedWeightingPotentialInGrid(electrode, z, x, y,t);
      }
    }
  }

  return ret;
}

void ComponentParallelPlate::DelayedWeightingField(
    const double x, const double y, const double z, const double t, double& wx,
    double& wy, double& wz, const std::string& label) {
    
    // Here I switch conventions back with the y-axis the direction of drift.
    
  wx = 0.;
  wy = 0.;
  wz = 0.;

  if (m_sigma == 0) {
    if (m_debug) {
      std::cout << m_className << "::DelayedWeightingField:\n"
                << "    Conductivity is set to zero.\n";
    }
    return;
  }

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      wx = electrode.flip *
           IntegrateDelayedField(electrode, fieldcomponent::ycomp, z, x, y, t);
      wy = electrode.flip *
           IntegrateDelayedField(electrode, fieldcomponent::zcomp, z, x, y, t);
      wz = electrode.flip *
           IntegrateDelayedField(electrode, fieldcomponent::xcomp, z, x, y, t);
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

void ComponentParallelPlate::AddPixel(double x, double z, double lx_input,
                                      double lz_input,
                                      const std::string& label) {
    
    // Here I switch conventions back with the y-axis the direction of drift.
    
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddPixel:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode pixel;
  pixel.label = label;
  pixel.ind = structureelectrode::Pixel;
  pixel.xpos = z;
  pixel.ypos = x;
  pixel.lx = lz_input;
  pixel.ly = lx_input;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(pixel));
  std::cout << m_className << "::AddPixel: Added pixel electrode.\n";
}

void ComponentParallelPlate::AddStrip(double z, double lz_input,
                                      const std::string& label) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddStrip:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode strip;
  strip.label = label;
  strip.ind = structureelectrode::Strip;
  strip.xpos = z;
  strip.lx = lz_input;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(strip));

  std::cout << m_className << "::AddStrip: Added strip electrode.\n";
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

  std::cout << m_className << "::AddPlane: Added plane electrode.\n";
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
    
    // TODO: Use existing classes for a grid based field map!
 
}

void ComponentParallelPlate::SetWeightingPotentialGrids(const double xmin, const double xmax,
const double xsteps, const double ymin,
const double ymax, const double ysteps,
const double zmin, const double zmax,
const double zsteps, const double tmin,
const double tmax, const double tsteps){
    
    // TODO: Use existing classes for a grid based field map!
   
}

double ComponentParallelPlate::FindWeightingPotentialInGrid(const Electrode& el,
                                                            const double x,
                                                            const double y,
                                                            const double z) {
    
    // TODO: Use existing classes for a grid based field map!
    
    return 0.;
}

double ComponentParallelPlate::FindDelayedWeightingPotentialInGrid(
    const Electrode& el, const double x, const double y, const double z,
    const double t) {
    
    // TODO: Use existing classes for a grid based field map!
    
    return 0.;
}

}  // namespace Garfield
