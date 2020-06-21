// Interface to Magboltz (version 9)

#ifndef G_MAGBOLTZ_INTERFACE
#define G_MAGBOLTZ_INTERFACE

#ifndef __CINT__

namespace Garfield {

namespace Magboltz {

constexpr unsigned int nEnergySteps = 4000;
constexpr unsigned int nMaxIonisationTerms = 30;
constexpr unsigned int nMaxInelasticTerms = 250;
constexpr unsigned int nMaxAttachmentTerms = 8;
constexpr unsigned int nMaxNullTerms = 10;
constexpr unsigned int nMaxLevelsPerComponent = 300;
constexpr unsigned int nCharDescr = 50;
constexpr unsigned int nMaxLevels = 960;
constexpr unsigned int nMaxComponents = 6;

extern "C" {

// Magboltz COMMON blocks

// Magnetic field
extern struct {
  double eovb;
  double wb;
  double btheta, bmag;
} bfld_;

extern struct {
  long long nGas;
  long long nStep;
  long long nAniso;
  double efinal;
  double estep;
  double akt;
  double ary;
  double tempc;
  double torr;
  long long ipen;
} inpt_;

extern struct {
  double tmax;
  double small;
  double api;
  double estart;
  double theta, phi;
  double rstart;
  double efield;
  long long nmax;
} setp_;

extern struct {
  double amgas[6];
  double vtmb[6];
  double tcfmx;
  double tcfmxg[6];
  long long ithrm;
} thrm_;

// Physical constants
extern struct {
  double echarg;
  double emass;
  double amu;
  double pir2;
} cnsts_;

extern struct {
  double eg[nEnergySteps];
  double eroot[nEnergySteps];
  double qt1[nEnergySteps];
  double qt2[nEnergySteps];
  double qt3[nEnergySteps];
  double qt4[nEnergySteps];
} mix2_;

extern struct { double den[nEnergySteps]; } dens_;

extern struct {
  double time[300];
  long long icoll[30];
  double spec[nEnergySteps];
  double tmax1;
  double ave;
  double den;
  double xid;
  double x;
  double y;
  double z;
  double st;
  long long nnull;
  long long icoln[nMaxLevels];
  long long icolnn[60];
} outpt_;

extern struct {
  double time[300];
  long long icoll[5][nMaxComponents];
  double spec[nEnergySteps];
  double tmax1;
  double ave;
  double den;
  double xid;
  double x;
  double y;
  double z;
  double st;
  long long nnull;
  long long icoln[290][nMaxComponents];
  long long icolnn[10][nMaxComponents];
} outptt_;

extern struct {
  char dscrpt[nMaxLevels][nCharDescr];
  char dscrptn[60][nCharDescr];
} scrip_;

extern struct {
  char dscrpt[nMaxLevelsPerComponent][nMaxComponents][nCharDescr];
  char dscrptn[10][nMaxComponents][nCharDescr];
} script_;

extern struct {
  double cf[nMaxLevels][nEnergySteps];
  double ein[nMaxLevels];
  double tcf[nEnergySteps];
  long long iarry[nMaxLevels];
  double rgas[nMaxLevels];
  double ipn[nMaxLevels];
  double wpl[nMaxLevels];
  long long last;
  long long isize;
  double penfra[nMaxLevels][3];
  double tcfmax[8];
} large_;

extern struct {
  double cf[290][nEnergySteps][nMaxComponents];
  double ein[290][nMaxComponents];
  double tcf[nEnergySteps][nMaxComponents];
  long long iarry[290][nMaxComponents];
  double rgas[290][nMaxComponents];
  double ipn[290][nMaxComponents];
  double wpl[290][nMaxComponents];
  long long last[nMaxComponents];
  long long isize[nMaxComponents];
  double penfra[290][3][nMaxComponents];
  double tcfmax[nMaxComponents];
} larget_;

// Definition of the gas mixture
extern struct { long long ngasn[6]; } gasn_;

extern struct {
  double an1, an2, an3, an4, an5, an6, an;
  double frac[6];
} ratio_;

// Calculation results
// Drift velocity
extern struct { double wx, wy, wz; } vel_;
extern struct { double dwx, dwy, dwz; } velerr_;

// Diffusion
extern struct {
  double difxx, difyy, difzz;
  double difyz, difxy, difxz;
} diflab_;
extern struct {
  double dxxer, dyyer, dzzer;
  double dyzer, dxyer, dxzer;
} diferb_;
extern struct { double difln, diftr; } difvel_;
extern struct { double dfler, dfter; } diferl_;

// Townsend and attachment coefficient
extern struct { double alpha, att; } ctowns_;
extern struct { double alper, atter; } ctwner_;
extern struct {
  double ralpha, ralper;
  double tofene, tofener, tofwv, tofwver;
  double tofdl, tofdler, tofdt, tofdter;
  double tofwr, tofwrer;
  double rattof, ratofer;
} tofout_;

void gasmix_(long long* ngs, double* q, double* qin, long long* nin, double* e,
             double* ei, char* name, double* virl, double* eb, double* peqel,
             double* peqin, double* penfra, long long* kel, long long* kin,
             double* qion, double* peqion, double* eion, long long* nion,
             double* qatt, long long* natt, double* qnull, long long* nnull,
             double* scln, long long* nc0, double* ec0, double* wk, double* efl,
             long long* ng1, double* eg1, long long* ng2, double* eg2,
             char scrpt[nMaxLevelsPerComponent][nCharDescr],
             char scrptn[nMaxNullTerms][nCharDescr]);

void colf_(double* freq, double* freel, double* freion, double* freatt, 
           double* frein, long long *ntotal);

void colft_(double* freq, double* freel, double* freion, double* freatt, 
            double* frein, long long *ntotal);

void magboltz_();
}
}
}
#endif
#endif
