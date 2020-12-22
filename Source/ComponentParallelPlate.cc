#include "Garfield/ComponentParallelPlate.hh"

#include <TF1.h>
#include <TF2.h>

#include <cmath>
#include <iostream>

#include "Garfield/GarfieldConstants.hh"

#define LOG(x) std::cout << m_className << ":: " << x << "\n" << std::endl

namespace Garfield {

ComponentParallelPlate::ComponentParallelPlate() : Component("ParallelPlate") {}
void ComponentParallelPlate::Setup(double g, double b, double eps, double V) {
  m_g = g;
  m_b = b;
  m_eps = eps;
  m_V = V;
  // For large times the resistive layer will act as a perfect conductor.
  m_ezg = -m_V / m_g;
  m_ezb = 0.;

  LOG("Geometry set.");
}

void ComponentParallelPlate::Setup(double g, double b, double eps, double V,
                                   double sigma) {
  m_g = g;
  m_b = b;
  m_eps = eps;
  m_V = V;
  m_sigma = sigma;

  m_ezg = -m_eps * m_V / (m_b + m_eps * m_g);
  m_ezb = -m_V / (m_b + m_eps * m_g);

  LOG("Geometry set.");
}

double ComponentParallelPlate::IntegrateField(const Electrode& el, int comp,
                                              const double X, const double Y,
                                              const double Z) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      return m_eps * m_Vw / (m_b + m_eps * m_g);
      break;
    }
    case structureelectrode::Pixel: {
      auto WFieldPixel = [](double* k, double* p) {
        double kx = k[0];
        double ky = k[1];

        int comp = (int)p[0];
        double X = p[1];
        double Y = p[2];
        double Z = p[3];
        double EPS = p[5];
        double G = p[6];
        double B = p[7];
        double LX = p[8];
        double LY = p[9];

        double K = std::sqrt(kx * kx + ky * ky);

        double intsol = 1;

        switch (comp) {
          case fieldcomponent::xcomp: {
            intsol *= 1 / (ky * cosh(G * K) * sinh(B * K) +
                           EPS * ky * cosh(B * K) * sinh(G * K));

            intsol *= cos(ky * Y) * sin((kx * LX) / 2) * sin((ky * LY) / 2) *
                      sin(kx * X) * sinh(K * (G - Z));
            break;
          }
          case fieldcomponent::ycomp: {
            intsol *= 1 / (kx * cosh(G * K) * sinh(B * K) +
                           EPS * kx * cosh(B * K) * sinh(G * K));
            intsol *= sin(ky * Y) * sin((kx * LX) / 2) * sin((ky * LY) / 2) *
                      cos(kx * X) * sinh(K * (G - Z));
            break;
          }
          case fieldcomponent::zcomp: {
            intsol *= 1 / (ky * kx * cosh(G * K) * sinh(B * K) +
                           EPS * ky * kx * cosh(B * K) * sinh(G * K));
            intsol *= K * cos(ky * Y) * sin((kx * LX) / 2) *
                      sin((ky * LY) / 2) * cos(kx * X) * cosh(K * (G - Z));
            break;
          }
        }
        return intsol;
      };

      TF2* fw =
          new TF2("WFieldPixel", WFieldPixel, 0, 10 * m_g, 0, 10 * m_g, 10);

      fw->SetParameters(comp, X - el.xpos, Y - el.ypos, Z, m_Vw, m_eps, m_g,
                        m_b, el.lx, el.ly);  // Set a = 0, b = 1
      double sol = fw->Integral(0, 10 * m_g, 0, 10 * m_g, 1.e-6);

      delete fw;

      return (4 * m_eps * m_Vw / (Pi * Pi)) * sol;

      break;
    }
    case structureelectrode::Strip: {
      auto WFieldStrip = [](double* k, double* p) {
        double kk = k[0];
        int comp = (int)p[0];
        double X = p[1];
        double Z = p[2];
        double EPS = p[4];
        double G = p[5];
        double B = p[6];
        double LX = p[7];

        double intsol =
            1 /
            ((cosh(G * kk) * sinh(B * kk) + EPS * cosh(B * kk) * sinh(G * kk)));
        switch (comp) {
          case fieldcomponent::xcomp: {
            intsol *= (sin(kk * LX / 2) * sin(kk * X) * sinh(kk * (G - Z)));
            break;
          }
          case fieldcomponent::zcomp: {
            intsol *= (sin(kk * LX / 2) * cos(kk * X) * cosh(kk * (G - Z)));
            break;
          }
        }

        return intsol;
      };

      if (comp == ycomp) {
        return 0.;
      } else {
        TF1* fw = new TF1("WFieldStrip", WFieldStrip, 0, 10 * m_g, 8);
        fw->SetParameters(comp, X - el.xpos, Z, m_Vw, m_eps, m_g, m_b,
                          el.lx);  // Set a = 0, b = 1

        double sol = fw->Integral(0, 10 * m_g);

        delete fw;

        return (2 * m_eps * m_Vw / (Pi)) * sol;
      }
      break;
    }
    default: {
      LOG("The structure of the electrode is not set");
      return 0.;
    }
  }
}

double ComponentParallelPlate::IntegratePromptPotential(const Electrode& el,
                                                        const double X,
                                                        const double Y,
                                                        const double Z) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      return m_eps * m_Vw * (m_g - Z) / (m_b + m_eps * m_g);
      break;
    }
    case structureelectrode::Pixel: {
      auto WPotentialPixel = [](double* k, double* p) {
        double kx = k[0];
        double ky = k[1];

        double X = p[0];
        double Y = p[1];
        double Z = p[2];
        double EPS = p[4];
        double G = p[5];
        double B = p[6];
        double LX = p[7];
        double LY = p[8];

        double K = std::sqrt(kx * kx + ky * ky);

        double intsol = 1;

        intsol *=
            cos(kx * X) * sin(kx * LX / 2) * cos(ky * Y) * sin(ky * LY / 2) *
            sinh(K * (G - Z)) /
            (kx * ky *
             (sinh(B * K) * cosh(G * K) + EPS * sinh(G * K) * cosh(B * K)));

        return intsol;
      };

      TF2* pw = new TF2("WPotentialPixel", WPotentialPixel, 0, 10 * m_g, 0,
                        10 * m_g, 9);

      pw->SetParameters(X - el.xpos, Y - el.ypos, Z, m_Vw, m_eps, m_g, m_b,
                        el.lx, el.ly);  // Set a = 0, b = 1
      double sol = pw->Integral(0, 2 * m_g, 0, 2 * m_g, 1.e-6);

      delete pw;

      return (4 * m_eps * m_Vw / (Pi * Pi)) * sol;
      break;
    }
    case structureelectrode::Strip: {
      auto WPotentialStrip = [](double* k, double* p) {
        double kk = k[0];
        double X = p[0];
        double Z = p[1];
        double EPS = p[3];
        double G = p[4];
        double B = p[5];
        double LX = p[6];

        double intsol = 1 / (kk * (cosh(G * kk) * sinh(B * kk) +
                                   EPS * cosh(B * kk) * sinh(G * kk)));
        intsol *= (sin(kk * LX / 2) * cos(kk * X) * sinh(kk * (G - Z)));

        return intsol;
      };

      TF1* pw = new TF1("WPotentialStrip", WPotentialStrip, 0, 10 * m_g, 7);

      pw->SetParameters(X - el.xpos, Z, m_Vw, m_eps, m_g, m_b,
                        el.lx);  // Set a = 0, b = 1
      double sol = pw->Integral(0, 10 * m_g);
      delete pw;
      return (2 * m_eps * m_Vw / (Pi)) * sol;
      break;
    }
    default: {
      LOG("The structure of the electrode is not set");
      return 0.;
    }
  }
}

double ComponentParallelPlate::IntegrateDelayedPotential(const Electrode& el,
                                                         const double X,
                                                         const double Y,
                                                         const double Z,
                                                         const double t) {
  switch (el.ind) {
    case structureelectrode::Plane: {
      double tau =
          (m_eps + m_b / m_g) /
          (m_sigma);  // Note to self: You dropt the eps) here for convenience.
      return m_Vw * (1 - exp(-t / tau)) * (m_b * (m_g - Z)) /
             (m_g *
              (m_b + m_eps * m_g));  // Note to self: Check epsr depencences.
      break;
    }
    case structureelectrode::Pixel: {
      auto WPotentialPixel = [](double* k, double* p) {
        double kx = k[0];
        double ky = k[1];

        double X = p[0];
        double Y = p[1];
        double Z = p[2];
        double EPS = p[4];
        double G = p[5];
        double B = p[6];
        double LX = p[7];
        double LY = p[8];
        double T = p[9];
        double SIGMA = p[10];

        double K = std::sqrt(kx * kx + ky * ky);
        double tau =
            (1 + cosh(G * K) * sinh(B * K) / (cosh(B * K) * sinh(G * K))) *
            (1 /
             SIGMA);  // Note to self: You dropt the eps) here for convenience.

        double intsol =
            1 / (kx * ky *
                 (sinh(B * K) * cosh(G * K) + EPS * sinh(G * K) * cosh(B * K)));

        intsol *= cos(kx * X) * sin(kx * LX / 2) * cos(ky * Y) *
                  sin(ky * LY / 2) * sinh(K * (G - Z)) * tanh(B * K) *
                  cosh(G * K) * (1 - exp(-T / tau)) / sinh(G * K);

        return intsol;
      };

      TF2* pw = new TF2("WPotentialPixel", WPotentialPixel, 0, 10 * m_g, 0,
                        10 * m_g, 11);

      pw->SetParameters(X - el.xpos, Y - el.ypos, Z, m_Vw, m_eps, m_g, m_b,
                        el.lx, el.ly, t, m_sigma);  // Set a = 0, b = 1
      double sol = pw->Integral(0, 2 * m_g, 0, 2 * m_g, 1.e-20);

      delete pw;

      return (4 * m_Vw / (Pi * Pi)) * sol;
      break;
    }
    case structureelectrode::Strip: {
      auto WPotentialStrip = [](double* k, double* p) {
        double kk = k[0];

        double X = p[0];
        double Z = p[1];
        double EPS = p[3];
        double G = p[4];
        double B = p[5];
        double LX = p[6];
        double T = p[7];
        double SIGMA = p[8];

        double tau =
            (EPS +
             cosh(G * kk) * sinh(B * kk) / (cosh(B * kk) * sinh(G * kk))) *
            (1 /
             SIGMA);  // Note to self: You dropt the eps) here for convenience.

        double intsol = 1 / (kk * (cosh(G * kk) * sinh(B * kk) +
                                   EPS * cosh(B * kk) * sinh(G * kk)));
        intsol *= (sin(kk * LX / 2) * cos(kk * X) * sinh(kk * (G - Z)) *
                   cosh(G * kk) * tanh(B * kk)) *
                  (1 - exp(-T / tau)) / sinh(G * kk);

        return intsol;
      };

      TF1* pw = new TF1("WPotentialStrip", WPotentialStrip, 0, 10 * m_g, 9);

      pw->SetParameters(X - el.xpos, Z, m_Vw, m_eps, m_g, m_b, el.lx, t,
                        m_sigma);  // Set a = 0, b = 1
      double sol = pw->Integral(0, 8 * m_g);
      delete pw;
      return (2 * m_Vw / (Pi)) * sol;
      break;
    }
    default: {
      LOG("The structure of the electrode is not set");
      return 0.;
    }
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double& ex,
                                           double& ey, double& ez, Medium*& m,
                                           int& status) {
  ex = ey = 0.;
  ez = -m_eps * m_V / (m_b + m_eps * m_g);

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

  v = -m_eps * m_V * (m_g - z) / (m_b + m_eps * m_g);

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
    LOG("Stopt drifting. (C2) ");
    status = -5;
  }
}

bool ComponentParallelPlate::GetVoltageRange(double& vmin, double& vmax) {
  return m_V == vmax - vmin;
}

void ComponentParallelPlate::WeightingField(const double x, const double y,
                                            const double z, double& wx,
                                            double& wy, double& wz,
                                            const std::string& label) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end()) {
    LOG("Label not found, please try again. (F)");
  } else {
    const auto index = it - m_readout.begin();

    auto element = m_readout_p[index];
    wx = IntegrateField(element, fieldcomponent::xcomp, x, y, z);
    wy = IntegrateField(element, fieldcomponent::ycomp, x, y, z);
    wz = IntegrateField(element, fieldcomponent::zcomp, x, y, z);
  }
}

double ComponentParallelPlate::WeightingPotential(const double x,
                                                  const double y,
                                                  const double z,
                                                  const std::string& label) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end()) {
    LOG("Label not found, please try again. (P)");
    return 0.;
  } else {
    const auto index = it - m_readout.begin();

    auto element = m_readout_p[index];
    return IntegratePromptPotential(element, x, y, z);
  }
}

double ComponentParallelPlate::DelayedWeightingPotential(
    const double x, const double y, const double z, const double t,
    const std::string& label) {
  if (m_sigma == 0) {
    LOG("No conductivity set!");
    return 0.;
  }
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end()) {
    LOG("Label not found, please try again. (P)");
    return 0.;
  } else {
    const auto index = it - m_readout.begin();

    auto element = m_readout_p[index];

    if (m_sigma == 0) {
      LOG("No conductivity set!");
      return 0.;
    }

    return IntegrateDelayedPotential(element, x, y, z, t);
  }
}

void ComponentParallelPlate::DelayedWeightingField(
    const double x, const double y, const double z, const double t, double& wx,
    double& wy, double& wz, const std::string& label) {
  wx = 0.;
  wy = 0.;
  wz = 0.;

  if (m_sigma == 0) {
    LOG("No conductivity set!");
  }
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end()) {
    LOG("Label not found, please try again. (P)");
  } else {
    const auto index = it - m_readout.begin();

    auto element = m_readout_p[index];

    if (m_sigma == 0) {
      LOG("No conductivity set!");
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
    LOG("Label " << label << " already excists, please use an other lable.");
  } else {
    Electrode* pixel = new Electrode();

    pixel->type = label;
    pixel->ind = structureelectrode::Pixel;
    pixel->xpos = x;
    pixel->ypos = y;
    pixel->lx = lx_input;
    pixel->ly = ly_input;

    m_readout.push_back(label);
    m_readout_p.push_back(*pixel);

    LOG("The pixel electrode structure has been added.");
  }
}

void ComponentParallelPlate::AddStrip(double x, double lx_input,
                                      const std::string& label) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end() && m_readout.size() > 0) {
    LOG("Label " << label << " already excists, please use an other lable.");
  } else {
    Electrode* strip = new Electrode();

    strip->type = label;
    strip->ind = structureelectrode::Strip;
    strip->xpos = x;
    strip->lx = lx_input;

    m_readout.push_back(label);
    m_readout_p.push_back(*strip);

    LOG("The strip electrode structure has been added.");
  }
}

void ComponentParallelPlate::AddPlane(const std::string& label) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end() && m_readout.size() > 0) {
    LOG("Label " << label << " already excists, please use an other lable.");
  } else {
    Electrode* plate = new Electrode();

    plate->type = label;
    plate->ind = structureelectrode::Plane;

    m_readout.push_back(label);
    m_readout_p.push_back(*plate);

    LOG("The plane electrode structure has been added.");
  }
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

}  // namespace Garfield
