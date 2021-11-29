#ifndef PARTICLE_DEF_H
#define PARTICLE_DEF_H

#include <string>
#include <list>

namespace Heed {

/// Definition of particles. Only the basic information: the name, the notation,
/// the mass, the charge, and other auxiliary data.
///
/// The organization is similar to AtomDef from directory matter, with
/// the exception that the internal data are not declared as private.
/// Of course, the user should not change them.
///
/// The principle of definitions of particles is dictionary or a database:
/// the particles are not repeated,
/// each particle is presented in the total system no more than one time.
/// The system knows each particle presented in it.
/// The particle characteristics can be obtained by literal notation.
/// The system declines the secondary initialization.
/// The copying is not declined.
/// When the user program wants to refer to particle,
/// it has to use either char* (string) notation, or pointer (or reference)
/// to one of these objects.
/// The user pogram can initialize the new particles.
/// The standard particles are initialized right here, below.
///
/// 1999 - 2004,   I. Smirnov

class particle_def {
 public:
  std::string name = "none";
  /// Short name to make data summary files short.
  std::string notation = "none";
  double mass = 0.;
  double charge = 0.;
  float spin = 0.;
  /// Default constructor.
  particle_def() = default;
  /// Constructor.
  particle_def(const std::string& fname, const std::string& fnotation,
               double fmass, double fcharge, float fspin);
  /// Constructor to create an anti-particle.
  particle_def(const std::string& fname, const std::string& fnotation,
               particle_def& p);
  /// Copy constructor.
  particle_def(const particle_def& f) 
      : name(f.name), notation(f.notation), mass(f.mass), charge(f.charge), 
        spin(f.spin) {
  }
  /// Assignment operator.
  particle_def& operator=(const particle_def&) = default;

  /// Destructor.
  ~particle_def() = default;

  /// Function for making an anti-particle.
  particle_def anti_particle(const particle_def& p);
  void print(std::ostream& file, int l) const;

  void set_mass(const double m);
  void set_charge(const double z);
};
std::ostream& operator<<(std::ostream& file, const particle_def& f);

extern particle_def electron_def;
extern particle_def positron_def;
extern particle_def muon_minus_def;
extern particle_def muon_plus_def;
extern particle_def proton_def;
extern particle_def anti_proton_def;
extern particle_def neutron_def;
extern particle_def anti_neutron_def;
extern particle_def P11_def;
extern particle_def D13_def;
extern particle_def S11_def;

// light unflavored mesons
extern particle_def pi_plus_meson_def;
extern particle_def pi_0_meson_def;
extern particle_def pi_minus_meson_def;
extern particle_def eta_meson_def;
extern particle_def K_plus_meson_def;
extern particle_def K_minus_meson_def;

extern particle_def deuteron_def;
extern particle_def alpha_particle_def;

// "exotic" particles with properties specified by user
extern particle_def user_particle_def;

}

#endif
