// Interface to Magboltz (version 9)

#ifndef G_MAGBOLTZ_INTERFACE
#define G_MAGBOLTZ_INTERFACE

#include <cstdint>

#ifndef __CINT__

#include "exports.h"

namespace Garfield {

namespace Magboltz {

constexpr unsigned int nEnergySteps = 4000;
constexpr unsigned int nMaxIonisationTerms = 30;
constexpr unsigned int nMaxInelasticTerms = 250;
constexpr unsigned int nMaxAttachmentTerms = 8;
constexpr unsigned int nMaxNullTerms = 10;
constexpr unsigned int nMaxLevelsPerComponent = 300;
constexpr unsigned int nCharName = 25;
constexpr unsigned int nCharDescr = 50;
constexpr unsigned int nMaxLevels = 960;
constexpr unsigned int nMaxComponents = 6;

extern "C" {

// Magboltz COMMON blocks

// Magnetic field
MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double eovb;
  double wb;
  double btheta, bmag;
} bfld_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  std::int64_t nGas;
  std::int64_t nStep;
  std::int64_t nAniso;
  double efinal;
  double estep;
  double akt;
  double ary;
  double tempc;
  double torr;
  std::int64_t ipen;
} inpt_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double tmax;
  double small;
  double api;
  double estart;
  double theta, phi;
  double rstart;
  double efield;
  std::int64_t nmax;
} setp_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double amgas[6];
  double vtmb[6];
  double tcfmx;
  double tcfmxg[6];
  std::int64_t ithrm;
} thrm_;

// Physical constants
MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double echarg;
  double emass;
  double amu;
  double pir2;
} cnsts_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double eg[nEnergySteps];
  double eroot[nEnergySteps];
  double qt1[nEnergySteps];
  double qt2[nEnergySteps];
  double qt3[nEnergySteps];
  double qt4[nEnergySteps];
} mix2_;

MEGABOLTZ_EXTERNAL_SYMBOL struct { double den[nEnergySteps]; } dens_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double time[300];
  std::int64_t icoll[30];
  double spec[nEnergySteps];
  double tmax1;
  double ave;
  double den;
  double xid;
  double x;
  double y;
  double z;
  double st;
  std::int64_t nnull;
  std::int64_t icoln[nMaxLevels];
  std::int64_t icolnn[60];
} outpt_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double time[300];
  std::int64_t icoll[5][nMaxComponents];
  double spec[nEnergySteps];
  double tmax1;
  double ave;
  double den;
  double xid;
  double x;
  double y;
  double z;
  double st;
  std::int64_t nnull;
  std::int64_t icoln[290][nMaxComponents];
  std::int64_t icolnn[10][nMaxComponents];
} outptt_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  char dscrpt[nMaxLevels][nCharDescr];
  char dscrptn[60][nCharDescr];
} scrip_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  char dscrpt[nMaxLevelsPerComponent][nMaxComponents][nCharDescr];
  char dscrptn[10][nMaxComponents][nCharDescr];
} script_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double cf[nMaxLevels][nEnergySteps];
  double ein[nMaxLevels];
  double tcf[nEnergySteps];
  std::int64_t iarry[nMaxLevels];
  double rgas[nMaxLevels];
  double ipn[nMaxLevels];
  double wpl[nMaxLevels];
  std::int64_t last;
  std::int64_t isize;
  double penfra[nMaxLevels][3];
  double tcfmax[8];
} large_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double cf[290][nEnergySteps][nMaxComponents];
  double ein[290][nMaxComponents];
  double tcf[nEnergySteps][nMaxComponents];
  std::int64_t iarry[290][nMaxComponents];
  double rgas[290][nMaxComponents];
  double ipn[290][nMaxComponents];
  double wpl[290][nMaxComponents];
  std::int64_t last[nMaxComponents];
  std::int64_t isize[nMaxComponents];
  double penfra[290][3][nMaxComponents];
  double tcfmax[nMaxComponents];
} larget_;

// Definition of the gas mixture
MEGABOLTZ_EXTERNAL_SYMBOL struct { std::int64_t ngasn[6]; } gasn_;

MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double an1, an2, an3, an4, an5, an6, an;
  double frac[6];
} ratio_;

// Calculation results
// Drift velocity
MEGABOLTZ_EXTERNAL_SYMBOL struct { double wx, wy, wz; } vel_;
MEGABOLTZ_EXTERNAL_SYMBOL struct { double dwx, dwy, dwz; } velerr_;

// Diffusion
MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double difxx, difyy, difzz;
  double difyz, difxy, difxz;
} diflab_;
MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double dxxer, dyyer, dzzer;
  double dyzer, dxyer, dxzer;
} diferb_;
MEGABOLTZ_EXTERNAL_SYMBOL struct { double difln, diftr; } difvel_;
MEGABOLTZ_EXTERNAL_SYMBOL struct { double dfler, dfter; } diferl_;

// Townsend and attachment coefficient
MEGABOLTZ_EXTERNAL_SYMBOL struct { double alpha, att; } ctowns_;
MEGABOLTZ_EXTERNAL_SYMBOL struct { double alper, atter; } ctwner_;
MEGABOLTZ_EXTERNAL_SYMBOL struct {
  double ralpha, ralper;
  double tofene, tofener, tofwv, tofwver;
  double tofdl, tofdler, tofdt, tofdter;
  double tofwr, tofwrer;
  double rattof, ratofer;
} tofout_;

MEGABOLTZ_EXTERNAL_SYMBOL void gasmix_(std::int64_t* ngs, double* q, double* qin, std::int64_t* nin, double* e,
             double* ei, char* name, double* virl, double* eb, double* peqel,
             double* peqin, double* penfra, std::int64_t* kel, std::int64_t* kin,
             double* qion, double* peqion, double* eion, std::int64_t* nion,
             double* qatt, std::int64_t* natt, double* qnull, std::int64_t* nnull,
             double* scln, std::int64_t* nc0, double* ec0, double* wk, double* efl,
             std::int64_t* ng1, double* eg1, std::int64_t* ng2, double* eg2,
             char scrpt[nMaxLevelsPerComponent][nCharDescr],
             char scrptn[nMaxNullTerms][nCharDescr],
             short namelen, short scrpt_len, short scrptn_len);

MEGABOLTZ_EXTERNAL_SYMBOL void colf_(double* freq, double* freel, double* freion, double* freatt,
           double* frein, std::int64_t *ntotal);

MEGABOLTZ_EXTERNAL_SYMBOL void colft_(double* freq, double* freel, double* freion, double* freatt,
            double* frein, std::int64_t *ntotal);

MEGABOLTZ_EXTERNAL_SYMBOL void magboltz_();
}
}
}
#endif
#endif
