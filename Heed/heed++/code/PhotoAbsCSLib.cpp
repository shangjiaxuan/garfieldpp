#include "heed++/code/PhotoAbsCSLib.h"
#include "heed++/code/PhysicalConstants.h"

#include <iostream>

// 2004, I. Smirnov

namespace {

std::string getDataBasePath() { 

  std::string path;
  // First try if the environment variable HEED_DATABASE is defined.
  char* heed_database = std::getenv("HEED_DATABASE");
  if (heed_database) {
    path = std::string(heed_database);
  } else {
    // If HEED_DATABASE is not defined, try GARFIELD_INSTALL.
    char* garfield_install = std::getenv("GARFIELD_INSTALL");
    if (garfield_install) {
      path = std::string(garfield_install) + "/share/Heed/database";
    } else {
      // Try GARFIELD_HOME.
      char* garfield_home = std::getenv("GARFIELD_HOME");
      if (garfield_home) {
        path = std::string(garfield_home) + "/Heed/heed++/database";
      } else {
        std::cerr << "Heed: Could not retrieve database path.\n";
      } 
    }
  }
  if (!path.empty()) {
    std::cout << "Heed:\n    Database path: " << path << "\n";
  }
  return path;
}

Heed::ExAtomPhotoAbsCS generate_Ar_PACS(const std::string& shelllist_dir,
                                        const std::string& pacs_table_dir) {

  Heed::ExAtomPhotoAbsCS Argon_PACS_mod_esc(18, shelllist_dir + "shelllist.dat",
                                            pacs_table_dir + "Ar.dat");

  // ExAtomPhotoAbsCS Argon_PACS_mod_esc(18,
  //                                     shelllist_dir + "shelllist.dat",
  //                                     shelllist_dir + "mw3.dat");

  // ExAtomPhotoAbsCS Argon_PACS_mod_esc(18, "argon",
  //                                     shelllist_dir + "ftbf18.dat", 2);

  Heed::AtomicSecondaryProducts* asp = Argon_PACS_mod_esc.get_asp(1);
  std::vector<double> electron_energy;
  std::vector<double> photon_energy;
  electron_energy.emplace_back(0.000200);
  // electron_energy.emplace_back(0.002670);
  asp->add_channel(0.65, electron_energy, photon_energy);
  electron_energy.resize(2);
  electron_energy[0] = 0.000050;
  electron_energy[1] = 0.000200;
  asp->add_channel(0.35, electron_energy, photon_energy, 1);
  // mcout<<"L1:\n";
  // asp->print(mcout, 2);

  asp = Argon_PACS_mod_esc.get_asp(2);
  electron_energy.resize(1);
  electron_energy[0] = 0.000200;
  asp->add_channel(1.0, electron_energy, photon_energy, 1);
  // mcout<<"L2:\n";
  // asp->print(mcout, 2);

  asp = Argon_PACS_mod_esc.get_asp(3);
  electron_energy.resize(1);
  electron_energy[0] = 0.000200;
  asp->add_channel(1.0, electron_energy, photon_energy, 1);
  // mcout<<"L3:\n";
  // asp->print(mcout, 2);

  return Argon_PACS_mod_esc;
}
}

namespace Heed {

using CLHEP::gram;
using CLHEP::mole;


std::map<std::string, ExAtomPhotoAbsCS> PhotoAbsCSLib::apacs;
std::map<std::string, SimpleAtomPhotoAbsCS> PhotoAbsCSLib::hpacs;

AtomPhotoAbsCS* PhotoAbsCSLib::getAPACS(const std::string& name) {
  if (apacs.empty()) initialise();
  if (name == "H" || name.find("H for") == 0) {
    return &hpacs[name];
  }
  if (apacs.count(name) > 0) return &apacs[name];
  return nullptr;
}

void PhotoAbsCSLib::initialise() {
  // Hydrogen
  hpacs.emplace("H", 
      SimpleAtomPhotoAbsCS(1, std::make_shared<HydrogenPhotoAbsCS>()));
  hpacs.emplace("H for H2", 
      SimpleAtomPhotoAbsCS(1, std::make_shared<PhenoPhotoAbsCS>("Hydrogen_for_H2", 1, 15.43e-6, 3.228)));
  hpacs.emplace("H for CH4", 
      SimpleAtomPhotoAbsCS(1, std::make_shared<PhenoPhotoAbsCS>("Hydrogen_for_CH4", 1, 12.65e-06, 3.228)));
  hpacs.emplace("H for NH4",
      SimpleAtomPhotoAbsCS(1, std::make_shared<PhenoPhotoAbsCS>("Hydrogen_for_NH4", 1, 10.0e-06, 3.228))); 
  // hpacs.emplace("H for CH4", 
  //     SimpleTablePhotoAbsCS("Hydrogen_for_CH4", 1, 12.65e-6,
  //                           shelllist_dir + "H_for_CH4.dat");

  const std::string shelllist_dir = getDataBasePath() + "/";
  const std::string pacs_table_dir = shelllist_dir + "henke/";

  std::string shells = shelllist_dir + "shelllist.dat"; 
  const std::map<std::string, int> atoms = {
    {"He",  2}, {"Li",  3}, {"Be",  4}, {"B",   5}, {"C",   6}, 
    {"O",   8}, {"F",   9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12},
    {"Al", 13}, {"Si", 14}, {"P",  15}, {"S",  16}, {"Cl", 17}, 
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Br", 35}, {"Kr", 36}, 
    {"Cd", 48}, {"Te", 52}, {"Xe", 54}, {"Cs", 55}, {"Hg", 80}, 
    {"U",  92}};
  for (const auto& atom : atoms) {
    std::string pacstable = pacs_table_dir + atom.first + ".dat";
    apacs.emplace(atom.first, ExAtomPhotoAbsCS(atom.second, shells, pacstable));
  }

  apacs.emplace("Ar", generate_Ar_PACS(shelllist_dir, pacs_table_dir));
  // "Standard" Argon:
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, shells,
  //                                      pacs_table_dir + "Ar.dat"));
  // Optional variants:
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, shells,
  //                                      shelllist_dir + "mw3.dat"));
  // Variant for debug, pointwise cross section
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, "argon",
  //                                      shelllist_dir + "ftbf18.dat", 2));
  // Variant for debug, fitted cross section
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, "argon",
  //                                      shelllist_dir + "shelltscf.dat",
  //                                      2, 0, 0.0);
  // Variant for debug, fitted cross section with replacement from Henke
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, "argon",
  //                                      shelllist_dir + "shelltscf.dat",
  //                                      pacs_table_dir + "Ar.dat",
  //                                      40.0e-6, 2, 0.0);
  // Another variant for debug, fitted cross section with replacement from
  // Marr and West, should be similar to old Fortran verion
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, "argon",
  //                                      shelllist_dir + "shelltscf.dat",
  //                                      shelllist_dir + "mw3.dat",
  //                                      40.0e-6, 2, 0.0);

  // For debug, FitBT
  // apacs.emplace("C", ExAtomPhotoAbsCS(6, "carbon",
  //                                     shelllist_dir + "shelltscf.dat",
  //                                     2, 0, 0.0);
  // apacs.emplace("C for CH4", ExAtomPhotoAbsCS(6, shells,
  //                                             pacs_table_dir + "C.dat",
  //                                             "C_for_CH4", 12.65e-06);

  apacs.emplace("C for CH4", ExAtomPhotoAbsCS(6, shells,
                                              shelllist_dir + "C_for_CH4.dat",
                                              "C_for_CH4", 12.65e-6));
  apacs.emplace("C for C2H4", ExAtomPhotoAbsCS(6, shells,
                                               pacs_table_dir + "C.dat",
                                               "C_for_C2H4", 10.51e-06));
  apacs.emplace("C for C2H6", ExAtomPhotoAbsCS(6, shells,
                                               pacs_table_dir + "C.dat",
                                               "C_for_C2H6", 11.52e-06));
  apacs.emplace("C_for_C4H10", ExAtomPhotoAbsCS(6, shells,
                                                pacs_table_dir + "C.dat",
                                                "C_for_C4H10", 10.55e-06));
  apacs.emplace("C for Methylal", ExAtomPhotoAbsCS(6, shells,
                                                   pacs_table_dir + "C.dat",
                                                   "C_for_Methylal", 10.0e-06));
  apacs.emplace("C for CF4", ExAtomPhotoAbsCS(6, shells,
                                              pacs_table_dir + "C.dat", 
                                              "C_for_CF4", 16.23e-06));
  apacs.emplace("C for CO2", ExAtomPhotoAbsCS(6, shells,
                                              pacs_table_dir + "C.dat", 
                                              "C_for_CO2", 13.79e-06));
  apacs.emplace("N", ExAtomPhotoAbsCS(7, shells,
                                      pacs_table_dir + "N.dat", 
                                      "N_for_N2", 15.581e-6));
  apacs.emplace("O for CO2", ExAtomPhotoAbsCS(8, shells,
                                              pacs_table_dir + "O.dat", 
                                              "O_for_CO2", 13.79e-6));

  std::string sshells = shelllist_dir + "shelllist_solid.dat";

  apacs.emplace("Diamond", ExAtomPhotoAbsCS(6, sshells,
                                            pacs_table_dir + "C.dat", 
                                            "Diamond"));
  apacs.emplace("Si crystal", ExAtomPhotoAbsCS(14, sshells,
                                               pacs_table_dir + "Si.dat",
                                               "Si_crystal"));
  apacs.emplace("Ge crystal", ExAtomPhotoAbsCS(32, shells,
                                               pacs_table_dir + "Ge.dat",
                                               "Ge_crystal", 0.67e-06));
  apacs.emplace("Si G4", ExAtomPhotoAbsCS(14, sshells,                                                                     shelllist_dir + "Si_G4.dat", "Si_G4"));
  apacs.emplace("Ga for GaAs", ExAtomPhotoAbsCS(31, sshells,
                                                pacs_table_dir + "Ga.dat",
                                                "Ga_for_GaAs"));
  apacs.emplace("As for GaAs", ExAtomPhotoAbsCS(33, sshells,
                                                pacs_table_dir + "As.dat",
                                                "As_for_GaAs"));
  apacs.emplace("Cd for CdTe", ExAtomPhotoAbsCS(48, sshells,
                                                pacs_table_dir + "Cd.dat",
                                                "Cd_for_CdTe"));
  apacs.emplace("Te for CdTe", ExAtomPhotoAbsCS(52, sshells,
                                                pacs_table_dir + "Te.dat",
                                                "Te_for_CdTe"));
}

MolecPhotoAbsCS H2_MPACS(PhotoAbsCSLib::getAPACS("H for H2"), 2);
// MolecPhotoAbsCS H2_MPACS(PhotoAbsCSLib::getAPACS("H"), 2);
MolecPhotoAbsCS He_MPACS(PhotoAbsCSLib::getAPACS("He"), 1, 41.3e-6);
MolecPhotoAbsCS Ne_MPACS(PhotoAbsCSLib::getAPACS("Ne"), 1, 35.4e-6);
MolecPhotoAbsCS Ar_MPACS(PhotoAbsCSLib::getAPACS("Ar"), 1, 26.4e-6);
MolecPhotoAbsCS Kr_MPACS(PhotoAbsCSLib::getAPACS("Kr"), 1, 24.4e-6);
MolecPhotoAbsCS Xe_MPACS(PhotoAbsCSLib::getAPACS("Xe"), 1, 22.1e-6);

MolecPhotoAbsCS N2_MPACS(PhotoAbsCSLib::getAPACS("N"), 2, 34.8e-6);
MolecPhotoAbsCS O2_MPACS(PhotoAbsCSLib::getAPACS("O"), 2, 30.8e-6);
MolecPhotoAbsCS NH3_MPACS(PhotoAbsCSLib::getAPACS("N"), 1, 
                          PhotoAbsCSLib::getAPACS("H for NH4"), 3, 26.6e-6);
MolecPhotoAbsCS N2O_MPACS(PhotoAbsCSLib::getAPACS("N"), 2, 
                          PhotoAbsCSLib::getAPACS("O"), 1, 34.8e-6);
MolecPhotoAbsCS CO2_MPACS(PhotoAbsCSLib::getAPACS("C for CO2"), 1, 
                          PhotoAbsCSLib::getAPACS("O for CO2"), 2, 33.0e-6);

// MolecPhotoAbsCS CH4_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 1,
//                           PhotoAbsCSLib::getAPACS("H for H2"), 4, 27.3e-6);
MolecPhotoAbsCS CH4_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 1, 
                          PhotoAbsCSLib::getAPACS("H for CH4"), 4, 27.3e-6);
MolecPhotoAbsCS CF4_MPACS(PhotoAbsCSLib::getAPACS("C for CF4"), 1, 
                          PhotoAbsCSLib::getAPACS("F"), 4);

// !!! The following line may need to be refined
// (to adjust outer shell energies).
MolecPhotoAbsCS SF4_MPACS(PhotoAbsCSLib::getAPACS("S"), 1, 
                          PhotoAbsCSLib::getAPACS("F"), 4);
// !!! The following line may need to be refined
// (to adjust outer shell energies).
MolecPhotoAbsCS SF6_MPACS(PhotoAbsCSLib::getAPACS("S"), 1, 
                          PhotoAbsCSLib::getAPACS("F"), 6);

MolecPhotoAbsCS C2H2_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 2, 
                           PhotoAbsCSLib::getAPACS("H for H2"), 2, 25.8e-6);
MolecPhotoAbsCS C2H4_MPACS(PhotoAbsCSLib::getAPACS("C for C2H4"), 2, 
                           PhotoAbsCSLib::getAPACS("H for H2"), 4, 25.8e-6);
MolecPhotoAbsCS C2H6_MPACS(PhotoAbsCSLib::getAPACS("C for C2H6"), 2, 
                           PhotoAbsCSLib::getAPACS("H for H2"), 6, 25.0e-6);
MolecPhotoAbsCS C3H8_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 3, 
                           PhotoAbsCSLib::getAPACS("H for H2"), 8, 24.0e-6);
MolecPhotoAbsCS C4H10_MPACS(PhotoAbsCSLib::getAPACS("C for C4H10"), 4, 
                            PhotoAbsCSLib::getAPACS("H for H2"), 10, 23.4e-6);

// !!! The following line may need to be refined
// (to adjust outer shell energies).
MolecPhotoAbsCS C2F4H2_MPACS(PhotoAbsCSLib::getAPACS("C for CF4"), 2, 
                             PhotoAbsCSLib::getAPACS("F"), 4,
                             PhotoAbsCSLib::getAPACS("H for H2"), 2);

// MolecPhotoAbsCS C2H2_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 2, 
//                            PhotoAbsCSLib::getAPACS("H for CH4"), 2);
// MolecPhotoAbsCS C2H4_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 2, 
//                            PhotoAbsCSLib::getAPACS("H for CH4"), 4);
// MolecPhotoAbsCS C2H6_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 2, 
//                            PhotoAbsCSLib::getAPACS("H for CH4"), 6);
// MolecPhotoAbsCS C3H8_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 3, 
//                            PhotoAbsCSLib::getAPACS("H for CH4"), 8);
// MolecPhotoAbsCS C4H10_MPACS(PhotoAbsCSLib::getAPACS("C for CH4"), 4, 
//                             PhotoAbsCSLib::getAPACS("H for CH4"), 10);
MolecPhotoAbsCS Methylal_MPACS(PhotoAbsCSLib::getAPACS("O"), 2, 
                               PhotoAbsCSLib::getAPACS("C for Methylal"), 3,
                               PhotoAbsCSLib::getAPACS("H for H2"), 8,
                               10.0e-6 * 23.4 / 10.55);  // similar to C4H10
/*
The value of W for noble gases is slightly less than
twice the ionization potential.
For organic gases it is very close to mean ionization potential
averaged with taking into account of atomic charges of carbon and hydrogen.
and assuming that the ionization potential of the hydrogen is the same
as in pure molecular hydrogen H2.
*/

// Additional molecular photoabsorption-cross sections
// for consistency with Magboltz
// Where available, the W values are taken from ICRU report 31
MolecPhotoAbsCS C5H12_MPACS(PhotoAbsCSLib::getAPACS("C for C4H10"), 5, 
                            PhotoAbsCSLib::getAPACS("H for H2"), 12, 23.2e-6);
MolecPhotoAbsCS H2O_MPACS(PhotoAbsCSLib::getAPACS("H for H2"), 2, 
                          PhotoAbsCSLib::getAPACS("O"), 1, 29.6e-6);
MolecPhotoAbsCS NO_MPACS(PhotoAbsCSLib::getAPACS("N"), 1, 
                         PhotoAbsCSLib::getAPACS("O"), 1);
MolecPhotoAbsCS CO_MPACS(PhotoAbsCSLib::getAPACS("C for CO2"), 1, 
                         PhotoAbsCSLib::getAPACS("O"), 1);
MolecPhotoAbsCS DME_MPACS(PhotoAbsCSLib::getAPACS("C for Methylal"), 2, 
                          PhotoAbsCSLib::getAPACS("H for H2"), 6,
                          PhotoAbsCSLib::getAPACS("O"), 1);
MolecPhotoAbsCS C2F6_MPACS(PhotoAbsCSLib::getAPACS("C for C2H6"), 2, 
                           PhotoAbsCSLib::getAPACS("F"), 6);
MolecPhotoAbsCS C3H6_MPACS(PhotoAbsCSLib::getAPACS("C for C2H6"), 3, 
                           PhotoAbsCSLib::getAPACS("H for H2"), 6);
MolecPhotoAbsCS CH3OH_MPACS(PhotoAbsCSLib::getAPACS("C for C2H6"), 1, 
                            PhotoAbsCSLib::getAPACS("H for H2"), 4,
                            PhotoAbsCSLib::getAPACS("O"), 1, 24.7e-6);
MolecPhotoAbsCS C2H5OH_MPACS(PhotoAbsCSLib::getAPACS("C for C2H6"), 2, 
                             PhotoAbsCSLib::getAPACS("H for H2"), 6,
                             PhotoAbsCSLib::getAPACS("O"), 1, 24.8e-6);
MolecPhotoAbsCS C3H7OH_MPACS(PhotoAbsCSLib::getAPACS("C for C2H6"), 3, 
                             PhotoAbsCSLib::getAPACS("H for H2"), 8,
                             PhotoAbsCSLib::getAPACS("O"), 1);
MolecPhotoAbsCS Cs_MPACS(PhotoAbsCSLib::getAPACS("Cs"), 1);
MolecPhotoAbsCS F2_MPACS(PhotoAbsCSLib::getAPACS("F"), 2);
MolecPhotoAbsCS CS2_MPACS(PhotoAbsCSLib::getAPACS("C for CO2"), 1, 
                          PhotoAbsCSLib::getAPACS("S"), 2);
MolecPhotoAbsCS COS_MPACS(PhotoAbsCSLib::getAPACS("C for CO2"), 1, 
                          PhotoAbsCSLib::getAPACS("O"), 1, 
                          PhotoAbsCSLib::getAPACS("S"), 1);
MolecPhotoAbsCS BF3_MPACS(PhotoAbsCSLib::getAPACS("B"), 1, 
                          PhotoAbsCSLib::getAPACS("F"), 3);
MolecPhotoAbsCS C2HF5_MPACS(PhotoAbsCSLib::getAPACS("C for C2H6"), 2, 
                            PhotoAbsCSLib::getAPACS("H for H2"), 1,
                            PhotoAbsCSLib::getAPACS("F"), 5);
MolecPhotoAbsCS C2H2F4_MPACS(PhotoAbsCSLib::getAPACS("C for C2H6"), 2, 
                             PhotoAbsCSLib::getAPACS("F"), 4,
                             PhotoAbsCSLib::getAPACS("H for H2"), 2);
MolecPhotoAbsCS CHF3_MPACS(PhotoAbsCSLib::getAPACS("C for CF4"), 1, 
                           PhotoAbsCSLib::getAPACS("H for H2"), 1,
                           PhotoAbsCSLib::getAPACS("F"), 3);
MolecPhotoAbsCS CF3Br_MPACS(PhotoAbsCSLib::getAPACS("C for CF4"), 1, 
                            PhotoAbsCSLib::getAPACS("F"), 3,
                            PhotoAbsCSLib::getAPACS("Br"), 1);
MolecPhotoAbsCS C3F8_MPACS(PhotoAbsCSLib::getAPACS("C for CF4"), 3, 
                           PhotoAbsCSLib::getAPACS("F"), 8);
MolecPhotoAbsCS O3_MPACS(PhotoAbsCSLib::getAPACS("O"), 3);
MolecPhotoAbsCS Hg_MPACS(PhotoAbsCSLib::getAPACS("Hg"), 1);
MolecPhotoAbsCS H2S_MPACS(PhotoAbsCSLib::getAPACS("H for H2"), 2, 
                          PhotoAbsCSLib::getAPACS("S"), 1);
MolecPhotoAbsCS GeH4_MPACS(PhotoAbsCSLib::getAPACS("Ge"), 1, 
                           PhotoAbsCSLib::getAPACS("H for H2"), 4);
MolecPhotoAbsCS SiH4_MPACS(PhotoAbsCSLib::getAPACS("Si"), 1, 
                           PhotoAbsCSLib::getAPACS("H for H2"), 4);
}
