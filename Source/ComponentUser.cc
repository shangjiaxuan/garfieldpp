#include <iostream>

#include "Garfield/ComponentUser.hh"

namespace Garfield {

ComponentUser::ComponentUser() : ComponentBase() {
  m_className = "ComponentUser";
}

void ComponentUser::ElectricField(const double x, const double y,
                                  const double z, double& ex, double& ey,
                                  double& ez, Medium*& m, int& status) {
  if (!m_efield) {
    ex = ey = ez = 0.;
    m = nullptr;
    status = -10;
    return;
  }

  m_efield(x, y, z, ex, ey, ez);
  m = GetMedium(x, y, z);
  if (!m) {
    if (m_debug) {
      std::cerr << m_className << "::ElectricField:\n    (" << x << ", " << y
                << ", " << z << ") is not inside a medium.\n";
    }
    status = -6;
    return;
  }

  status = m->IsDriftable() ? 0 : -5;
}

void ComponentUser::ElectricField(const double x, const double y,
                                  const double z, double& ex, double& ey,
                                  double& ez, double& v, Medium*& m,
                                  int& status) {
  if (!m_efield) {
    ex = ey = ez = v = 0.;
    m = nullptr;
    status = -10;
    return;
  }
  m_efield(x, y, z, ex, ey, ez);

  if (m_potential) {
    m_potential(x, y, z, v);
  } else {
    v = 0.;
  }

  m = GetMedium(x, y, z);
  if (!m) {
    if (m_debug) {
      std::cerr << m_className << "::ElectricField:\n    (" << x << ", " << y
                << ", " << z << ") is not inside a medium.\n";
    }
    status = -6;
    return;
  }

  status = m->IsDriftable() ? 0 : -5;
}

bool ComponentUser::GetVoltageRange(double& vmin, double& vmax) {
  vmin = vmax = 0.;
  return false;
}

void ComponentUser::MagneticField(const double x, const double y,
                                  const double z, double& bx, double& by,
                                  double& bz, int& status) {
  if (!m_bfield) {
    bx = by = bz = 0.;
    status = -10;
    return;
  }
  m_bfield(x, y, z, bx, by, bz);
  status = 0;
}

void ComponentUser::WeightingField(const double x, const double y,
                                   const double z, double& wx, double& wy,
                                   double& wz, const std::string& label) {
  wx = wy = wz = 0.;
  if (!m_wfield) return;
  m_wfield(x, y, z, wx, wy, wz, label);
}

double ComponentUser::WeightingPotential(const double x, const double y,
                                         const double z,
                                         const std::string& label) {
  double v = 0.;
  if (m_wpot) m_wpot(x, y, z, v, label);
  return v;
}

void ComponentUser::DelayedWeightingField(const double x, const double y, 
                                          const double z, const double t,
                                          double& wx, double& wy, double& wz,
                                          const std::string& label) {
  wx = wy = wz = 0.;
  if (m_dwfield) m_dwfield(x, y, z, t, wx, wy, wz, label);
}

void ComponentUser::SetElectricField(
    std::function<void(const double, const double, const double, 
                       double&, double&, double&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetElectricField: Function is empty.\n";
    return;
  }
  m_efield = f;
  m_ready = true;
}

void ComponentUser::SetPotential(
    std::function<void(const double, const double, const double, 
                       double&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetPotential: Function is empty.\n";
    return;
  }
  m_potential = f;
}

void ComponentUser::SetWeightingField(
    std::function<void(const double, const double, const double, 
                       double&, double&, double&, const std::string&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetWeightingField: Function is empty.\n";
    return;
  }
  m_wfield = f;
}

void ComponentUser::SetWeightingPotential(
    std::function<void(const double, const double, const double, 
                       double&, const std::string&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetWeightingPotential: Function is empty.\n";
    return;
  }
  m_wpot = f;
}

void ComponentUser::SetDelayedWeightingField(
    std::function<void(const double, const double, const double, const double,
                       double&, double&, double&, const std::string&)> f) {

  if (!f) {
    std::cerr << m_className << "::SetDelayedWeightingField: Function is empty.\n";
    return;
  }
  m_dwfield = f;
}

void ComponentUser::SetMagneticField(
    std::function<void(const double, const double, const double, 
                       double&, double&, double&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetMagneticField: Function is empty.\n";
    return;
  }
  m_bfield = f;
}

void ComponentUser::Reset() {
  m_efield = nullptr;
  m_potential = nullptr;
  m_wfield = nullptr;
  m_wpot = nullptr;
  m_dwfield = nullptr;
  m_bfield = nullptr;
  m_ready = false;
}

void ComponentUser::UpdatePeriodicity() {
  if (m_debug) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Periodicities are not supported.\n";
  }
}
}
