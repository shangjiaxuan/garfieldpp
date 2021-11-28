#ifndef PHOTOABSCSLIB_H
#define PHOTOABSCSLIB_H

#include <map>

#include "heed++/code/PhotoAbsCS.h"

namespace Heed {

/// Library of photoabsorption cross sections for some frequently used atoms
/// and molecules.
/// 2004, I. Smirnov

class PhotoAbsCSLib {

 public:
  static AtomPhotoAbsCS* getAPACS(const std::string& name); 

 private:
  static std::map<std::string, ExAtomPhotoAbsCS> apacs;
  static std::map<std::string, SimpleAtomPhotoAbsCS> hpacs;
  static void initialise();
};

extern MolecPhotoAbsCS H2_MPACS;
extern MolecPhotoAbsCS He_MPACS;
extern MolecPhotoAbsCS N2_MPACS;
extern MolecPhotoAbsCS O2_MPACS;
extern MolecPhotoAbsCS Ne_MPACS;
extern MolecPhotoAbsCS Ar_MPACS;
extern MolecPhotoAbsCS Kr_MPACS;
extern MolecPhotoAbsCS Xe_MPACS;
extern MolecPhotoAbsCS NH3_MPACS;
extern MolecPhotoAbsCS N2O_MPACS;
extern MolecPhotoAbsCS CO2_MPACS;
extern MolecPhotoAbsCS CH4_MPACS;
extern MolecPhotoAbsCS CF4_MPACS;

// The definition of the following in PhotoAbsLib.c may be refined
// (to adjust outer shell energies).
extern MolecPhotoAbsCS SF4_MPACS;
extern MolecPhotoAbsCS SF6_MPACS;

extern MolecPhotoAbsCS C2H2_MPACS;
extern MolecPhotoAbsCS C2H4_MPACS;
extern MolecPhotoAbsCS C2H6_MPACS;
extern MolecPhotoAbsCS C3H8_MPACS;
extern MolecPhotoAbsCS C4H10_MPACS;

// The definition of the following in PhotoAbsLib.c may be refined
// (to adjust outer shell energies).
extern MolecPhotoAbsCS C2F4H2_MPACS;
extern MolecPhotoAbsCS Methylal_MPACS;

// Additional molecular photoabsorption-cross sections
// for compatibility with Magboltz
extern MolecPhotoAbsCS C5H12_MPACS;
extern MolecPhotoAbsCS H2O_MPACS;
extern MolecPhotoAbsCS NO_MPACS;
extern MolecPhotoAbsCS CO_MPACS;
extern MolecPhotoAbsCS DME_MPACS;
extern MolecPhotoAbsCS C2F6_MPACS;
extern MolecPhotoAbsCS C3H6_MPACS;
extern MolecPhotoAbsCS CH3OH_MPACS;
extern MolecPhotoAbsCS C2H5OH_MPACS;
extern MolecPhotoAbsCS C3H7OH_MPACS;
extern MolecPhotoAbsCS Cs_MPACS;
extern MolecPhotoAbsCS F2_MPACS;
extern MolecPhotoAbsCS CS2_MPACS;
extern MolecPhotoAbsCS COS_MPACS;
extern MolecPhotoAbsCS BF3_MPACS;
extern MolecPhotoAbsCS C2HF5_MPACS;
extern MolecPhotoAbsCS C2H2F4_MPACS;
extern MolecPhotoAbsCS CHF3_MPACS;
extern MolecPhotoAbsCS CF3Br_MPACS;
extern MolecPhotoAbsCS C3F8_MPACS;
extern MolecPhotoAbsCS O3_MPACS;
extern MolecPhotoAbsCS Hg_MPACS;
extern MolecPhotoAbsCS H2S_MPACS;
extern MolecPhotoAbsCS GeH4_MPACS;
extern MolecPhotoAbsCS SiH4_MPACS;
}

#endif
