/*
(c) 2009 Supratik Mukhopadhyay, Nayana Majumdar
*/
// Gets number of primitives (polygonal surfaces (2D) and wires (1D)) and
// their details
// Assumptions:
// MaxNbVertices: 4 (maximum nb of vertices a polygon can have)
// Vertices have been counted starting from zero - it may be better to count
// them starting with 1
// Elements have a count from 0 to NbElement - 1: can create confusion

#define DEFINE_INTFACEGLOBAL
#define DEFINE_neBEMGLOBAL
#define DEFINE_NRGLOBAL

#include <stdio.h>
#include <sys/stat.h>  // use of stat function
#include <unistd.h>
#include <float.h>

#ifdef __cplusplus
#include <vector>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "neBEMInterface.h"
#include "Isles.h"
#include "NR.h"
#include "neBEM.h"

#ifdef __cplusplus
namespace neBEM {
#endif

// Called from a code requesting neBEM services
int neBEMInitialize(void) {

  // Version information
  strcpy(neBEMVersion, "1.9.08");
  strcpy(ISLESVersion, "1.4.8");
  printf("Using neBEM version %s and ISLES version %s\n", neBEMVersion,
         ISLESVersion);

  // The following is not an absolute necessity, but can ease the burden of the
  // user by setting default values of some of the important variable. Please
  // take a look at the src/Interface/ExampleDev2neBEM.c to find out what values
  // are being set, and what they mean.
  int fstatus = neBEMSetDefaults();
  if (fstatus != 0) {
    neBEMMessage("neBEMInitialize - neBEMSetDefaults");
    return -1;
  }

  // Change some of the global variables according to requirements.
  // Note that the solution flags and counters are set before the neBEMSolve
  // is invoked.
  LengthScale = 1.0;
  DebugLevel = 0;

  // The following function allows input files to be read and the geometry
  // set according to these files, as was customary in the previous versions
  if (OptDeviceFile) {
    printf("Reading geometry details from %s\n", DeviceInputFile);
    fstatus = neBEMGetInputsFromFiles();
    if (fstatus != 0) {
      neBEMMessage("neBEMInitialize - neBEMGetInputFromFiles");
      return -1;
    }
  }

  // creation of the mother output directory should be necessary only once
  if (neBEMState == 0) {
    fstatus = CreateDirStr();
    if (fstatus != 0) {
      neBEMMessage("neBEMInitialize - CreateDirStr");
      return -1;
    }
  }

  // Create Isles log file for keeping track of approximations (numerical
  // quadrature) in case evaluation of algebraic expressions fails.
  char IslesFile[256];
  strcpy(IslesFile, PPOutDir);
  strcat(IslesFile, "/Isles.log");
  fIsles = fopen(IslesFile, "w");
  if (fIsles == NULL) {
    neBEMMessage("neBEMInitialize - IslesFile");
    return -1;
  }

  // following integers keep track of the success / failure of the exact
  // expressions in estimating realistic physical properties.
  IslesCntr = ExactCntr = FailureCntr = ApproxCntr = 0;

  // Set up parameters related to neBEM computations
  int RqstdThreads = 1;
  FILE *processFile = fopen("neBEMProcess.inp", "r");
  if (processFile == NULL) {
    printf("neBEMProcess.inp absent ... assuming defaults ...\n");
    PrimAfter = 0;
    if (NbThreads > 0) RqstdThreads = NbThreads;
  } else {
    fscanf(processFile, "PrimAfter: %d\n", &PrimAfter);
    fscanf(processFile, "RqstdThreads: %d\n", &RqstdThreads);
    fclose(processFile);
  }

#ifdef _OPENMP
  int MaxProcessors = omp_get_num_procs();
  if (RqstdThreads > 1) {
    if (RqstdThreads < MaxProcessors) {
      // one processor left alone
      omp_set_num_threads(RqstdThreads);
    } else {
      printf("RqstdThreads: %d\n", RqstdThreads);
      RqstdThreads = MaxProcessors - 1;
      omp_set_num_threads(RqstdThreads);
      printf("Adjusted RqstdThreads: %d\n", RqstdThreads);
    }
  } else {
    // Work with one thread
    RqstdThreads = 1;  // cannot be zero or negative!
    omp_set_num_threads(RqstdThreads);
    printf("RqstdThreads: %d => No Multi-threading ...\n", RqstdThreads);
  }

  // OpenMP related information
  printf("PrimAfter: %d\n", PrimAfter);
  printf("RqstdThreads: %d, MaxProcessors: %d\n", RqstdThreads, MaxProcessors);
  printf("Maximum number of threads to be used for parallelization: %d\n",
         omp_get_max_threads()); 
  printf("Number of threads used for neBEMInitialize: %d\n", omp_get_num_threads());
#endif
  // Set up parameters related to voxelized data export for Garfield++
  FILE *voxelInpFile = fopen("neBEMVoxel.inp", "r");
  if (voxelInpFile == NULL) {
    printf("neBEMVoxel.inp absent ... assuming OptVoxel = 0 ...\n");
    OptVoxel = 0;
    OptStaggerVoxel = 0;
  } else {
    fscanf(voxelInpFile, "OptVoxel: %d\n", &OptVoxel);
    fscanf(voxelInpFile, "OptStaggerVoxel: %d\n", &OptStaggerVoxel);
    fscanf(voxelInpFile, "Xmin: %le\n", &Voxel.Xmin);
    fscanf(voxelInpFile, "Xmax: %le\n", &Voxel.Xmax);
    fscanf(voxelInpFile, "Ymin: %le\n", &Voxel.Ymin);
    fscanf(voxelInpFile, "Ymax: %le\n", &Voxel.Ymax);
    fscanf(voxelInpFile, "Zmin: %le\n", &Voxel.Zmin);
    fscanf(voxelInpFile, "Zmax: %le\n", &Voxel.Zmax);
    fscanf(voxelInpFile, "XStagger: %le\n", &Voxel.XStagger);
    fscanf(voxelInpFile, "YStagger: %le\n", &Voxel.YStagger);
    fscanf(voxelInpFile, "ZStagger: %le\n", &Voxel.ZStagger);
    fscanf(voxelInpFile, "NbOfXCells: %d\n", &Voxel.NbXCells);
    fscanf(voxelInpFile, "NbOfYCells: %d\n", &Voxel.NbYCells);
    fscanf(voxelInpFile, "NbOfZCells: %d\n", &Voxel.NbZCells);
    fclose(voxelInpFile);
  }  // inputs for Voxel

  // Set up parameters related to 3dMap data export for Garfield++
  FILE *mapInpFile = fopen("neBEMMap.inp", "r");
  if (mapInpFile == NULL) {
    printf("neBEMMap.inp absent ... assuming OptMap = 0 ...\n");
    OptMap = 0;
    OptStaggerMap = 0;
  } else {
    // While reading the input, OptMap and OptStaggerMap have to be read
    // first since that will decide whether there is a map and its version.
    fscanf(mapInpFile, "OptMap: %d\n", &OptMap);
    fscanf(mapInpFile, "OptStaggerMap: %d\n", &OptStaggerMap);
    fscanf(mapInpFile, "MapVersion: %9s\n", MapVersion);
    fscanf(mapInpFile, "Xmin: %le\n", &Map.Xmin);
    fscanf(mapInpFile, "Xmax: %le\n", &Map.Xmax);
    fscanf(mapInpFile, "Ymin: %le\n", &Map.Ymin);
    fscanf(mapInpFile, "Ymax: %le\n", &Map.Ymax);
    fscanf(mapInpFile, "Zmin: %le\n", &Map.Zmin);
    fscanf(mapInpFile, "Zmax: %le\n", &Map.Zmax);
    fscanf(mapInpFile, "XStagger: %le\n", &Map.XStagger);
    fscanf(mapInpFile, "YStagger: %le\n", &Map.YStagger);
    fscanf(mapInpFile, "ZStagger: %le\n", &Map.ZStagger);
    fscanf(mapInpFile, "NbOfXCells: %d\n", &Map.NbXCells);
    fscanf(mapInpFile, "NbOfYCells: %d\n", &Map.NbYCells);
    fscanf(mapInpFile, "NbOfZCells: %d\n", &Map.NbZCells);
    fclose(mapInpFile);
  }  // inputs for 3dMap

  // Set up parameters related to fast volume
  FILE *fastInpFile = fopen("neBEMFastVol.inp", "r");
  if (fastInpFile == NULL) {
    printf("neBEMFastVol.inp absent ... assuming OptFastVol = 0 ...\n");
    OptFastVol = 0;
    OptStaggerFastVol = 0;
    OptCreateFastPF = 0;
    OptReadFastPF = 0;
    FastVol.NbBlocks = 0;
    FastVol.NbOmitVols = 0;
    FastVol.NbIgnoreVols = 0;
  } else {
    fscanf(fastInpFile, "OptFastVol: %d\n", &OptFastVol);
    fscanf(fastInpFile, "OptStaggerFastVol: %d\n", &OptStaggerFastVol);
    fscanf(fastInpFile, "OptCreateFastPF: %d\n", &OptCreateFastPF);
    fscanf(fastInpFile, "OptReadFastPF: %d\n", &OptReadFastPF);
    fscanf(fastInpFile, "NbPtSkip: %d\n", &NbPtSkip);
    fscanf(fastInpFile, "NbStgPtSkip: %d\n", &NbStgPtSkip);
    fscanf(fastInpFile, "LX: %le\n", &FastVol.LX);
    fscanf(fastInpFile, "LY: %le\n", &FastVol.LY);
    fscanf(fastInpFile, "LZ: %le\n", &FastVol.LZ);
    fscanf(fastInpFile, "CornerX: %le\n", &FastVol.CrnrX);
    fscanf(fastInpFile, "CornerY: %le\n", &FastVol.CrnrY);
    fscanf(fastInpFile, "CornerZ: %le\n", &FastVol.CrnrZ);
    fscanf(fastInpFile, "YStagger: %le\n", &FastVol.YStagger);
    if (!OptStaggerFastVol)
      FastVol.YStagger = 0.0;  // ignore any non-zero value
    fscanf(fastInpFile, "NbOfBlocks: %d\n", &FastVol.NbBlocks);
    BlkNbXCells = ivector(1, FastVol.NbBlocks);
    BlkNbYCells = ivector(1, FastVol.NbBlocks);
    BlkNbZCells = ivector(1, FastVol.NbBlocks);
    BlkLZ = dvector(1, FastVol.NbBlocks);
    BlkCrnrZ = dvector(1, FastVol.NbBlocks);
    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      fscanf(fastInpFile, "NbOfXCells: %d\n", &BlkNbXCells[block]);
      fscanf(fastInpFile, "NbOfYCells: %d\n", &BlkNbYCells[block]);
      fscanf(fastInpFile, "NbOfZCells: %d\n", &BlkNbZCells[block]);
      fscanf(fastInpFile, "LZ: %le\n", &BlkLZ[block]);
      fscanf(fastInpFile, "CornerZ: %le\n", &BlkCrnrZ[block]);
    }  // inputs for blocks
    fscanf(fastInpFile, "NbOfOmitVols: %d\n", &FastVol.NbOmitVols);
    if (FastVol.NbOmitVols) {
      OmitVolLX = dvector(1, FastVol.NbOmitVols);
      OmitVolLY = dvector(1, FastVol.NbOmitVols);
      OmitVolLZ = dvector(1, FastVol.NbOmitVols);
      OmitVolCrnrX = dvector(1, FastVol.NbOmitVols);
      OmitVolCrnrY = dvector(1, FastVol.NbOmitVols);
      OmitVolCrnrZ = dvector(1, FastVol.NbOmitVols);
      for (int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
        fscanf(fastInpFile, "OmitVolLX: %le\n", &OmitVolLX[omit]);
        fscanf(fastInpFile, "OmitVolLY: %le\n", &OmitVolLY[omit]);
        fscanf(fastInpFile, "OmitVolLZ: %le\n", &OmitVolLZ[omit]);
        fscanf(fastInpFile, "OmitVolCornerX: %le\n", &OmitVolCrnrX[omit]);
        fscanf(fastInpFile, "OmitVolCornerY: %le\n", &OmitVolCrnrY[omit]);
        fscanf(fastInpFile, "OmitVolCornerZ: %le\n", &OmitVolCrnrZ[omit]);
      }  // inputs for OmitVols
    }    // inputs for OmitVols
    fscanf(fastInpFile, "NbOfIgnoreVols: %d\n", &FastVol.NbIgnoreVols);
    if (FastVol.NbIgnoreVols) {
      IgnoreVolLX = dvector(1, FastVol.NbIgnoreVols);
      IgnoreVolLY = dvector(1, FastVol.NbIgnoreVols);
      IgnoreVolLZ = dvector(1, FastVol.NbIgnoreVols);
      IgnoreVolCrnrX = dvector(1, FastVol.NbIgnoreVols);
      IgnoreVolCrnrY = dvector(1, FastVol.NbIgnoreVols);
      IgnoreVolCrnrZ = dvector(1, FastVol.NbIgnoreVols);
      for (int ignore = 1; ignore <= FastVol.NbIgnoreVols; ++ignore) {
        fscanf(fastInpFile, "IgnoreVolLX: %le\n", &IgnoreVolLX[ignore]);
        fscanf(fastInpFile, "IgnoreVolLY: %le\n", &IgnoreVolLY[ignore]);
        fscanf(fastInpFile, "IgnoreVolLZ: %le\n", &IgnoreVolLZ[ignore]);
        fscanf(fastInpFile, "IgnoreVolCornerX: %le\n", &IgnoreVolCrnrX[ignore]);
        fscanf(fastInpFile, "IgnoreVolCornerY: %le\n", &IgnoreVolCrnrY[ignore]);
        fscanf(fastInpFile, "IgnoreVolCornerZ: %le\n", &IgnoreVolCrnrZ[ignore]);
      }  // inputs for IgnoreVols
    }    // inputs for IgnoreVols
    for (int ignore = 1; ignore <= FastVol.NbIgnoreVols; ++ignore) {
      printf("IgnoreVolLX: %le\n", IgnoreVolLX[ignore]);
      printf("IgnoreVolLY: %le\n", IgnoreVolLY[ignore]);
      printf("IgnoreVolLZ: %le\n", IgnoreVolLZ[ignore]);
      printf("IgnoreVolCornerX: %le\n", IgnoreVolCrnrX[ignore]);
      printf("IgnoreVolCornerY: %le\n", IgnoreVolCrnrY[ignore]);
      printf("IgnoreVolCornerZ: %le\n", IgnoreVolCrnrZ[ignore]);
    }  // inputs for IgnoreVols
    fclose(fastInpFile);
  }  // else fastInpFile

  // Set up parameters related to fixed specification of weighting field
  FILE *fixedWtInpFile = fopen("neBEMFixedWtField.inp", "r");
  if (fixedWtInpFile == NULL) {
    printf(
        "neBEMFixedWtField.inp absent ... assuming OptFixedWtField = 0 ...\n");
    OptFixedWtField = 0;
    FixedWtPotential = 0.0;
    FixedWtFieldX = 0.0;
    FixedWtFieldY = 0.0;
    FixedWtFieldZ = 0.0;
  } else {
    fscanf(fixedWtInpFile, "OptFixedWtField: %d\n", &OptFixedWtField);
    fscanf(fixedWtInpFile, "FixedWtPotential: %lg\n", &FixedWtPotential);
    fscanf(fixedWtInpFile, "FixedWtFieldX: %lg\n", &FixedWtFieldX);
    fscanf(fixedWtInpFile, "FixedWtFieldY: %lg\n", &FixedWtFieldY);
    fscanf(fixedWtInpFile, "FixedWtFieldZ: %lg\n", &FixedWtFieldZ);
    fclose(fixedWtInpFile);
  }  // else fixedWtInpFile

  // Set up parameters related to weighting field fast volume
  FILE *fastWtFldInpFile = fopen("neBEMWtFldFastVol.inp", "r");
  if (fastWtFldInpFile == NULL) {
    printf(
        "neBEMWtFldFastVol.inp absent ... assuming OptWtFldFastVol = 0 ...\n");
    OptWtFldFastVol = 0;
    OptWtFldStaggerFastVol = 0;
    OptWtFldCreateFastPF = 0;
    OptWtFldReadFastPF = 0;
    WtFldFastVol.NbBlocks = 0;
    WtFldFastVol.NbOmitVols = 0;
    WtFldFastVol.NbIgnoreVols = 0;
  } else {
    fscanf(fastWtFldInpFile, "OptFastVol: %d\n", &OptWtFldFastVol);
    fscanf(fastWtFldInpFile, "OptStaggerFastVol: %d\n",
           &OptWtFldStaggerFastVol);
    fscanf(fastWtFldInpFile, "OptCreateFastPF: %d\n", &OptWtFldCreateFastPF);
    fscanf(fastWtFldInpFile, "OptReadFastPF: %d\n", &OptWtFldReadFastPF);
    fscanf(fastWtFldInpFile, "NbPtSkip: %d\n", &WtFldNbPtSkip);
    fscanf(fastWtFldInpFile, "NbStgPtSkip: %d\n", &WtFldNbStgPtSkip);
    fscanf(fastWtFldInpFile, "LX: %le\n", &WtFldFastVol.LX);
    fscanf(fastWtFldInpFile, "LY: %le\n", &WtFldFastVol.LY);
    fscanf(fastWtFldInpFile, "LZ: %le\n", &WtFldFastVol.LZ);
    fscanf(fastWtFldInpFile, "CornerX: %le\n", &WtFldFastVol.CrnrX);
    fscanf(fastWtFldInpFile, "CornerY: %le\n", &WtFldFastVol.CrnrY);
    fscanf(fastWtFldInpFile, "CornerZ: %le\n", &WtFldFastVol.CrnrZ);
    fscanf(fastWtFldInpFile, "YStagger: %le\n", &WtFldFastVol.YStagger);
    if (!OptWtFldStaggerFastVol)
      WtFldFastVol.YStagger = 0.0;  // ignore any non-zero value
    fscanf(fastWtFldInpFile, "NbOfBlocks: %d\n", &WtFldFastVol.NbBlocks);
    WtFldBlkNbXCells = ivector(1, WtFldFastVol.NbBlocks);
    WtFldBlkNbYCells = ivector(1, WtFldFastVol.NbBlocks);
    WtFldBlkNbZCells = ivector(1, WtFldFastVol.NbBlocks);
    WtFldBlkLZ = dvector(1, WtFldFastVol.NbBlocks);
    WtFldBlkCrnrZ = dvector(1, WtFldFastVol.NbBlocks);
    for (int block = 1; block <= WtFldFastVol.NbBlocks; ++block) {
      fscanf(fastWtFldInpFile, "NbOfXCells: %d\n", &WtFldBlkNbXCells[block]);
      fscanf(fastWtFldInpFile, "NbOfYCells: %d\n", &WtFldBlkNbYCells[block]);
      fscanf(fastWtFldInpFile, "NbOfZCells: %d\n", &WtFldBlkNbZCells[block]);
      fscanf(fastWtFldInpFile, "LZ: %le\n", &WtFldBlkLZ[block]);
      fscanf(fastWtFldInpFile, "CornerZ: %le\n", &WtFldBlkCrnrZ[block]);
    }  // inputs for blocks
    fscanf(fastWtFldInpFile, "NbOfOmitVols: %d\n", &WtFldFastVol.NbOmitVols);
    if (WtFldFastVol.NbOmitVols) {
      WtFldOmitVolLX = dvector(1, WtFldFastVol.NbOmitVols);
      WtFldOmitVolLY = dvector(1, WtFldFastVol.NbOmitVols);
      WtFldOmitVolLZ = dvector(1, WtFldFastVol.NbOmitVols);
      WtFldOmitVolCrnrX = dvector(1, WtFldFastVol.NbOmitVols);
      WtFldOmitVolCrnrY = dvector(1, WtFldFastVol.NbOmitVols);
      WtFldOmitVolCrnrZ = dvector(1, WtFldFastVol.NbOmitVols);
      for (int omit = 1; omit <= WtFldFastVol.NbOmitVols; ++omit) {
        fscanf(fastWtFldInpFile, "OmitVolLX: %le\n", &WtFldOmitVolLX[omit]);
        fscanf(fastWtFldInpFile, "OmitVolLY: %le\n", &WtFldOmitVolLY[omit]);
        fscanf(fastWtFldInpFile, "OmitVolLZ: %le\n", &WtFldOmitVolLZ[omit]);
        fscanf(fastWtFldInpFile, "OmitVolCornerX: %le\n",
               &WtFldOmitVolCrnrX[omit]);
        fscanf(fastWtFldInpFile, "OmitVolCornerY: %le\n",
               &WtFldOmitVolCrnrY[omit]);
        fscanf(fastWtFldInpFile, "OmitVolCornerZ: %le\n",
               &WtFldOmitVolCrnrZ[omit]);
      }  // inputs for OmitVols
    }    // inputs for OmitVols
    fscanf(fastWtFldInpFile, "NbOfIgnoreVols: %d\n",
           &WtFldFastVol.NbIgnoreVols);
    if (WtFldFastVol.NbIgnoreVols) {
      WtFldIgnoreVolLX = dvector(1, WtFldFastVol.NbIgnoreVols);
      WtFldIgnoreVolLY = dvector(1, WtFldFastVol.NbIgnoreVols);
      WtFldIgnoreVolLZ = dvector(1, WtFldFastVol.NbIgnoreVols);
      WtFldIgnoreVolCrnrX = dvector(1, WtFldFastVol.NbIgnoreVols);
      WtFldIgnoreVolCrnrY = dvector(1, WtFldFastVol.NbIgnoreVols);
      WtFldIgnoreVolCrnrZ = dvector(1, WtFldFastVol.NbIgnoreVols);
      for (int ignore = 1; ignore <= WtFldFastVol.NbIgnoreVols; ++ignore) {
        fscanf(fastWtFldInpFile, "IgnoreVolLX: %le\n",
               &WtFldIgnoreVolLX[ignore]);
        fscanf(fastWtFldInpFile, "IgnoreVolLY: %le\n",
               &WtFldIgnoreVolLY[ignore]);
        fscanf(fastWtFldInpFile, "IgnoreVolLZ: %le\n",
               &WtFldIgnoreVolLZ[ignore]);
        fscanf(fastWtFldInpFile, "IgnoreVolCornerX: %le\n",
               &WtFldIgnoreVolCrnrX[ignore]);
        fscanf(fastWtFldInpFile, "IgnoreVolCornerY: %le\n",
               &WtFldIgnoreVolCrnrY[ignore]);
        fscanf(fastWtFldInpFile, "IgnoreVolCornerZ: %le\n",
               &WtFldIgnoreVolCrnrZ[ignore]);
      }  // inputs for IgnoreVols
    }    // inputs for IgnoreVols
    for (int ignore = 1; ignore <= WtFldFastVol.NbIgnoreVols; ++ignore) {
      printf("WtFldIgnoreVolLX: %le\n", WtFldIgnoreVolLX[ignore]);
      printf("WtFldIgnoreVolLY: %le\n", WtFldIgnoreVolLY[ignore]);
      printf("WtFldIgnoreVolLZ: %le\n", WtFldIgnoreVolLZ[ignore]);
      printf("WtFldIgnoreVolCornerX: %le\n", WtFldIgnoreVolCrnrX[ignore]);
      printf("WtFldIgnoreVolCornerY: %le\n", WtFldIgnoreVolCrnrY[ignore]);
      printf("WtFldIgnoreVolCornerZ: %le\n", WtFldIgnoreVolCrnrZ[ignore]);
    }  // inputs for IgnoreVols
    fclose(fastWtFldInpFile);
  }  // else fastWtFldInpFile

  printf("neBEM initialized ...\n");
  fflush(stdout);
  sleep(3);  // wait for three seconds so that the user gets time to react

  neBEMState = 1;  // state 1 implied initialization of neBEM completed

  // announce success - later, add the name of the calling code
  return (0);
}  // neBEMInitialize ends

// Reads geometry details
// Note that reflection and periodicity (reflection) information is gathered
// using the neBEMGetPeriodicities() function, called from here.
// Repetition variables were introduced to facilitate GarfieldInterface.
// Now Garfield can pass parameters directly to neBEM and, hence, these
// variables have become redundant.
// They are likely to be removed in near future.
int neBEMReadGeometry(void) {
  int dbgFn = 0;
  int fstatus;

  startClock = clock();

  // For a model that was defined before and for which data was stored in a file
  if ((!NewModel) && (!NewBC) && (OptStorePrimitives)) {
    fstatus = ReadPrimitives();
    if (fstatus) {
      neBEMMessage("neBEMReadGeometry - problem reading stored Primitives.\n");
      return -1;
    }
    neBEMState = 3;  // primitives read in after initialization and Nbs
    return 0;
  }

  printf("geometry inputs ...\n");
  if (neBEMState != 1) {
    printf("reading geometry possible only after initialization ...\n");
    return -1;
  }

  NbPrimitives = neBEMGetNbPrimitives();
  OrgnlNbPrimitives = NbPrimitives;
  if (NbPrimitives == 0) {
    // nothing to do - return control to calling routine
    neBEMMessage("neBEMReadGeometry - no primitive.\n");
    return (-1);  // for the time being
  }

  NbSurfs = 0;
  NbWires = 0;
  neBEMState = 2;

  // Allocate memory for storing the geometry primitives till the elements are
  // created.
  // Explicit storage of these variables related to primitives may not be
  // necessary if the elements are created immediately after a primitive is read
  // in.
  // MaxNbVertices = 4; // value specified through SetDefaults or init files

  // neBEM has been initialized, NbPrimitives set
  PrimType = ivector(1, NbPrimitives);
  NbVertices = ivector(1, NbPrimitives);
  OrgnlToEffPrim = imatrix(1, NbPrimitives, 0, 2);  // 0 init, 1 intrfc, 2 rmv
  XVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  YVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  ZVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  XNorm = dvector(1, NbPrimitives);
  YNorm = dvector(1, NbPrimitives);
  ZNorm = dvector(1, NbPrimitives);
  PrimLX = dvector(1, NbPrimitives);
  PrimLZ = dvector(1, NbPrimitives);
  Radius = dvector(1, NbPrimitives);  // can lead to a little memory misuse
  PrimOriginX = dvector(1, NbPrimitives);
  PrimOriginY = dvector(1, NbPrimitives);
  PrimOriginZ = dvector(1, NbPrimitives);
  PrimDC = (DirnCosn3D *)malloc(NbPrimitives * sizeof(DirnCosn3D));
  VolRef1 = ivector(1, NbPrimitives);
  VolRef2 = ivector(1, NbPrimitives);
  NbSurfSegX = ivector(1, NbPrimitives);
  NbSurfSegZ = ivector(1, NbPrimitives);
  NbWireSeg = ivector(1, NbPrimitives);  // little memory misuse
  InterfaceType = ivector(1, NbPrimitives);
  Epsilon1 = dvector(1, NbPrimitives);
  Epsilon2 = dvector(1, NbPrimitives);
  Lambda = dvector(1, NbPrimitives);
  ApplPot = dvector(1, NbPrimitives);
  ApplCh = dvector(1, NbPrimitives);
  PeriodicTypeX = ivector(1, NbPrimitives);
  PeriodicTypeY = ivector(1, NbPrimitives);
  PeriodicTypeZ = ivector(1, NbPrimitives);
  PeriodicInX = ivector(1, NbPrimitives);
  PeriodicInY = ivector(1, NbPrimitives);
  PeriodicInZ = ivector(1, NbPrimitives);
  XPeriod = dvector(1, NbPrimitives);
  YPeriod = dvector(1, NbPrimitives);
  ZPeriod = dvector(1, NbPrimitives);
  MirrorTypeX = ivector(1, NbPrimitives);
  MirrorTypeY = ivector(1, NbPrimitives);
  MirrorTypeZ = ivector(1, NbPrimitives);
  MirrorDistXFromOrigin = dvector(1, NbPrimitives);
  MirrorDistYFromOrigin = dvector(1, NbPrimitives);
  MirrorDistZFromOrigin = dvector(1, NbPrimitives);
  BndPlaneInXMin = ivector(1, NbPrimitives);
  BndPlaneInXMax = ivector(1, NbPrimitives);
  BndPlaneInYMin = ivector(1, NbPrimitives);
  BndPlaneInYMax = ivector(1, NbPrimitives);
  BndPlaneInZMin = ivector(1, NbPrimitives);
  BndPlaneInZMax = ivector(1, NbPrimitives);
  XBndPlaneInXMin = dvector(1, NbPrimitives);
  XBndPlaneInXMax = dvector(1, NbPrimitives);
  YBndPlaneInYMin = dvector(1, NbPrimitives);
  YBndPlaneInYMax = dvector(1, NbPrimitives);
  ZBndPlaneInZMin = dvector(1, NbPrimitives);
  ZBndPlaneInZMax = dvector(1, NbPrimitives);
  VBndPlaneInXMin = dvector(1, NbPrimitives);
  VBndPlaneInXMax = dvector(1, NbPrimitives);
  VBndPlaneInYMin = dvector(1, NbPrimitives);
  VBndPlaneInYMax = dvector(1, NbPrimitives);
  VBndPlaneInZMin = dvector(1, NbPrimitives);
  VBndPlaneInZMax = dvector(1, NbPrimitives);
  NbElmntsOnPrim = ivector(1, NbPrimitives);
  ElementBgn = ivector(1, NbPrimitives);
  ElementEnd = ivector(1, NbPrimitives);
  AvChDen = dvector(1, NbPrimitives);
  AvAsgndChDen = dvector(1, NbPrimitives);

  // Loop over the primitives - major loop
  int nvertex, volref1, volref2, volmax = 0;
#ifdef __cplusplus
  std::vector<double> xvert(MaxNbVertices, 0.);
  std::vector<double> yvert(MaxNbVertices, 0.);
  std::vector<double> zvert(MaxNbVertices, 0.);
#else
  double xvert[MaxNbVertices], yvert[MaxNbVertices], zvert[MaxNbVertices];
#endif
  double xnorm, ynorm, znorm;  // in case of wire , radius is read as xnorm
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
#ifdef __cplusplus
    fstatus = neBEMGetPrimitive(prim, &nvertex, 
                                xvert.data(), yvert.data(), zvert.data(),
                                &xnorm, &ynorm, &znorm, &volref1, &volref2);
#else
    fstatus = neBEMGetPrimitive(prim, &nvertex, xvert, yvert, zvert,  // arrays
                                &xnorm, &ynorm, &znorm, &volref1, &volref2);
#endif
    if (fstatus != 0) {
      neBEMMessage("neBEMReadGeometry - neBEMGetPrimitve");
      return -1;
    }
    if (volmax < volref1) {
      volmax = volref1;
    }  // maxm nb of volumes
    if (volmax < volref2) {
      volmax = volref2;
    }  // maxm nb of volumes

    if (nvertex > MaxNbVertices) {
      printf("Number of vertices for primitive %d exceeds %d!\n", prim,
             MaxNbVertices);
      printf("Returning to garfield ...\n");
      return (-1);
    }

    PrimType[prim] = nvertex;  // wire:2, triangle:3, rectangle:4
    NbVertices[prim] = nvertex;
    for (int vert = 0; vert < NbVertices[prim]; ++vert) {
      XVertex[prim][vert] = xvert[vert];
      YVertex[prim][vert] = yvert[vert];
      ZVertex[prim][vert] = zvert[vert];
    }
    if (PrimType[prim] == 2) {
      // wire
      XNorm[prim] = 0.0;  // modulus not 1 - an absurd trio!
      YNorm[prim] = 0.0;
      ZNorm[prim] = 0.0;
      Radius[prim] = xnorm;
    }
    if ((PrimType[prim] == 3) || (PrimType[prim] == 4)) {
      XNorm[prim] = xnorm;
      YNorm[prim] = ynorm;
      ZNorm[prim] = znorm;
      Radius[prim] = 0.0;  // absurd radius!
    }
    VolRef1[prim] = volref1;
    VolRef2[prim] = volref2;

    // feedback for user begins - suppress later
    if (OptPrintPrimaryDetails) {
      printf("neBEM:\tprimitive %d between volumes %d, %d has %d vertices\n",
             prim, volref1, volref2, nvertex);
      for (int ivertex = 0; ivertex < nvertex; ivertex++) {
        printf("\tnode %d (%g,%g,%g)\n", ivertex, xvert[ivertex],
               yvert[ivertex], zvert[ivertex]);
      }
      printf("\tnormal vector: (%g,%g,%g)\n", xnorm, ynorm, znorm);
    }
    // feedback for user ends - suppress later

    // Now look for the volume related information for this primitve
    // This is obtained from the specified volume references
    // volref1 refers to the volume itself
    // volref2 describes the volume in the direction of the +ve normal
    // Note that materials from 1 to 10 are conductors and
    // 										 from 11 to 20 are
    // dielectrics
    int shape1, material1, boundarytype1;
    double epsilon1, potential1, charge1;
    if (volref1 == -1) {
      // Must be an error, since no device is made of vacuum
      neBEMMessage("neBEMReadGeometry - volref1 = -1!");
      return -1;
    } else {
      neBEMVolumeDescription(volref1, &shape1, &material1, &epsilon1,
                             &potential1, &charge1, &boundarytype1);
    }
    if (OptPrintVolumeDetails) {
      printf("\tvolref1: %d\n", volref1);
      printf("\t\tboundarytype1: %d, shape1: %d, material1: %d\n",
             boundarytype1, shape1, material1);
      printf("\t\tepsilon1: %lg, potential1: %lg, charge1: %lg\n", epsilon1,
             potential1, charge1);
    }
    // in the -ve normal direction - properties of the external volume
    int shape2, material2, boundarytype2;
    double epsilon2, potential2, charge2;
    if (volref2 == -1) {
      shape2 = 0;
      material2 = 11;
      epsilon2 = 1.0;
      potential2 = 0.0;
      charge2 = 0.0;
      boundarytype2 = 0;
    } else {
      neBEMVolumeDescription(volref2, &shape2, &material2, &epsilon2,
                             &potential2, &charge2, &boundarytype2);
    }
    if (OptPrintVolumeDetails) {
      printf("\tvolref2: %d\n", volref2);
      printf("\t\tboundarytype2: %d, shape2: %d, material2: %d\n",
             boundarytype2, shape2, material2);
      printf("\t\tepsilon2: %lg, potential2: %lg, charge2: %lg\n", epsilon2,
             potential2, charge2);
    }

    // Put default values to variables that depend on the interface type
    // Is there any risk involved in putting these defaults?
    // At present, they even seem necessary. For example, for floating
    // conductors or dielectric-dielectric interface, the formulation requires
    // that the RHS is zero (may be modified by the effect of known charges).
    Epsilon1[prim] = epsilon1;
    Epsilon2[prim] = epsilon2;  // 1: self, 2: external
    ApplPot[prim] = 0.0;
    Lambda[prim] = 0.0;
    ApplCh[prim] = 0.0;

    // BoundaryTypes:
    // --------------
    // Vacuum: 0
    // Conductor at specified potential: 1
    // Conductor with a specified charge: 2
    // Floating conductor (zero charge, perpendicular E): 3
    // Dielectric intrface (plastic-plastic) without "manual" charge: 4
    // Dielectric with surface charge (plastic-gas, typically): 5
    // Symmetry boundary, E parallel: 6 (may not be necessary)
    // Symmetry boundary, E perpendicular: 7 (may not be necessary)

    // InterfaceType:
    // --------------
    // To be skipped: 0
    // Conductor-dielectric: 1
    // Conductor with known charge: 2
    // Conductor at floating potential: 3
    // Dielectric-dielectric: 4
    // Dielectric with given surface charge: 5
    // Check dielectric-dielectric formulation in
    // (NumSolnOfBIEforMolES_JPBardhan.pdf):
    // Numerical solution of boundary-integral equations for molecular
    // electrostatics,
    // by Jaydeep P. Bardhan,
    // THE JOURNAL OF CHEMICAL PHYSICS 130, 094102 (2009)

    switch (boundarytype1) {  // the volume itself is volref1
      case 1:                 // conductor at specified potential
        if (boundarytype2 == 0 || boundarytype2 == 4) {
          // dielectric-conductor
          InterfaceType[prim] = 1;
          ApplPot[prim] = potential1;
        } else if (boundarytype2 == 1) {
          // conductor-conductor
          if (fabs(potential1 - potential2)  // same potential
              < 1e-6 * (1 + fabs(potential1) + fabs(potential2))) {
            printf("neBEMReadGeometry: identical potentials; skipped.\n");
            printf("Primitive skipped: #%d\n", prim);
            InterfaceType[prim] = 0;
          } else {
            // different potentials
            printf("neBEMReadGeometry: different potentials; rejected.\n");
            return -1;
          }
        } else {
          // conductor-unknown
          printf(
              "neBEMReadGeometry: conductor at given potential; rejected.\n");
          return -1;
        }
        break;

      case 2:  // conductor with a specified charge
        if ((boundarytype2 == 0) || (boundarytype2 == 4)) {
          // conductor-dielectric
          InterfaceType[prim] = 2;
          ApplCh[prim] = charge1;
        } else {
          printf("neBEMReadGeometry: charged conductor; rejected.\n");
          return -1;
        }
        break;

      case 3:  // floating conductor (zero charge, perpendicular E)
        if ((boundarytype2 == 0) || (boundarytype2 == 4)) {
          // conductor-dielectric
          InterfaceType[prim] = 3;
          if (!NbFloatingConductors)   // assuming only one floating conductor
            NbFloatingConductors = 1;  // in the system
        } else {
          printf("neBEMReadGeometry: floating conductor; rejected.\n");
          return -1;
        }
        break;

      case 4:  // dielectric interface (plastic-plastic) without "manual" charge
        if (boundarytype2 == 0) {
          // dielectric-vacuum
          // epsilon1 is self dielectric-constant
          // epsilon2 is towards positive normal
          InterfaceType[prim] = 4;  
          Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
          // consistent with Bardhan's eqn 16 where (1 / (2*Lambda)) is used
        } else if (boundarytype2 == 1) {
          // dielectric-conductor
          InterfaceType[prim] = 1;  // conductor at known potential
          ApplPot[prim] = potential2;
        } else if (boundarytype2 == 2) {
          // dielectric-conductor
          InterfaceType[prim] = 2;  // conductor with known charge
          ApplCh[prim] = charge2;
        } else if (boundarytype2 == 3) {
          // dielectric-conductor
          InterfaceType[prim] = 3;  // conductor at floating potential
        } else if (boundarytype2 == 4) {
          // dielectric-dielectric
          if (fabs(epsilon1 - epsilon2) < 1e-6 * (1 + fabs(epsilon1) + fabs(epsilon2))) {
            // identical dielectrica
            printf(
                "neBEMReadGeometry: between identical dielectrica; skipd.\n");
            printf("Primitive skipped: #%d\n", prim);
            InterfaceType[prim] = 0;
          } else {
            // distinctly different dielectrica
            // epsilon1 is self dielectric-constant
            // epsilon2 towards positive normal
            InterfaceType[prim] = 4;  
            Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
            // consistent with Bardhan's paper (1 / Lambda)
          }  
        } else if (boundarytype2 == 5) {
          // dielectric-dielectric with charge
          if (fabs(epsilon1 - epsilon2)  // identical dielectrica
              < 1e-6 * (1 + fabs(epsilon1) + fabs(epsilon2))) {
            printf(
                "neBEMReadGeometry: between identical dielectrica; skipped.\n");
            printf("Primitive skipped: #%d\n", prim);
            InterfaceType[prim] = 0;
          } else {
            // distinctly different dielectrica
            InterfaceType[prim] = 5;  // epsilon2 towards positive normal
            ApplCh[prim] = charge2;
            Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
          }
        }       // if-else if boundarytypes 0 and 4
        else {  // dielectric-unknown
          printf("neBEMReadGeometry: unknown dielectric; rejected.\n");
          return -1;
        }
        break;

      case 5:  // dielectric with surface charge (plastic-gas, typically)
        if (boundarytype2 == 0) { // dielectric-vacuum
          InterfaceType[prim] = 5;  // epsilon2 is towards +ve normal
          ApplCh[prim] = charge1;
          Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
          // consistent with Bardhan's paper (1 / Lambda)
        } else if (boundarytype2 == 4) { // dielectric-dielectric
          if (fabs(epsilon1 - epsilon2) < 1e-6 * (1 + fabs(epsilon1) + fabs(epsilon2))) {
            // identical dielectrica
            printf(
                "neBEMReadGeometry: between identical dielectrica; skipd.\n");
            printf("Primitive skipped: #%d\n", prim);
            InterfaceType[prim] = 0;
          } else {
            // distinctly different dielectrica
            InterfaceType[prim] = 5;  // epsilon2 towards positive normal
            ApplCh[prim] = charge1;
            Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
            // consistent with Bardhan's paper (1 / Lambda)
          }  
        }    // if-else if boundarytypes 0 and 4
        else {
          printf(
              "neBEMReadGeometry: charged dielectric adjacent to a conductor; "
              "rejected.\n");
          return -1;
        }
        break;

      case 6:  // symmetry boundary, E parallel
        if (boundarytype2 == 0) {
          InterfaceType[prim] = 6;
        } else {
          printf("neBEMReadGeometry: E-parallel symmetry; rejected.\n");
          return -1;
        }
        break;

      case 7:  // symmetry boundary, E perpendicular
        if (boundarytype2 == 0) {
          InterfaceType[prim] = 7;
        } else {
          printf("neBEMReadGeometry: E-perpendicular symmetry; rejected.\n");
          return -1;
        }
        break;

      default:
        printf("neBEMReadGeometry: Boundary type 1: %d\n", boundarytype1);
        printf("neBEMReadGeometry: Boundary type 2: %d\n", boundarytype2);
        printf("neBEMReadGeometry:          out of range ... exiting.\n");
        return -1;
    }  // switch boundarytype1 ends

    if (OptPrintVolumeDetails) {
      printf(
          "\tType: %d, ApplPot: %lg, Epsilon1: %lg, Epsilon2: %lg, Lambda: "
          "%lg, ApplCh: %lg\n",
          InterfaceType[prim], ApplPot[prim], Epsilon1[prim], Epsilon2[prim],
          Lambda[prim], ApplCh[prim]);
    }

    // Read the periodicities
    // Note that mirror has been taken care of separately below
    // ix: PeriodicTypeX (1 - simple, 2 - mirror, 3 - axial, 4 - rotation)
    // jx: PeriodicInX (Number of copies internal to neBEM)
    // sx: XPeriod
    // NOTE: A change in this part of the algorithm is likely to affect
    // src/Solve/neBEM.c (LHMatrix) and
    // src/Solve/ComputeProperties.c (PFAtPoint and WtPFAtPoint)
    {
      int ix, iy, iz;
      int jx, jy, jz;
      double sx, sy, sz;
      fstatus = neBEMGetPeriodicities(prim, &ix, &jx, &sx, &iy, &jy, &sy, &iz,
                                      &jz, &sz);
      if (fstatus != 0) {
        neBEMMessage("neBEMReadGeometry - neBEMGetPeriodicities");
        return -1;
      }
      if (jx < 0) jx = 0;
      if (jy < 0) jy = 0;
      if (jz < 0) jz = 0;

      PeriodicTypeX[prim] = ix;
      PeriodicTypeY[prim] = iy;
      PeriodicTypeZ[prim] = iz;
      if (0) {
        printf("For primitive: %d\n", prim);
        printf("\tPeriodicTypeX: %d, PeriodicTypeY: %d, PeriodicTypeZ: %d\n",
               ix, iy, iz);
        printf("\tPeriodicInX: %d, PeriodicInY: %d, PeriodicInZ: %d\n", jx, jy,
               jz);
        printf("\tXPeriod: %lg, YPeriod: %lg, ZPeriod: %lg\n", sx, sy, sz);
      }
      if (ix > 0) {
        // These checks need to be done separately. Otherwise, there is
        // a possibility that non-zero values of PeriodicIn* and *Period
        // are used throughout the code despite PeriodicType* is 0
        PeriodicInX[prim] = jx;  
        XPeriod[prim] = sx;
      } else {
        PeriodicInX[prim] = 0;
        XPeriod[prim] = 0.0;
      }
      if (iy > 0) {
        PeriodicInY[prim] = jy;
        YPeriod[prim] = sy;
      } else {
        PeriodicInY[prim] = 0;
        YPeriod[prim] = 0.0;
      }
      if (iz > 0) {
        PeriodicInZ[prim] = jz;
        ZPeriod[prim] = sz;
      } else {
        PeriodicInZ[prim] = 0;
        ZPeriod[prim] = 0.0;
      }
    }  // read periodicity information

    // Read mirror information
    // Mirror can be in perpendicular to any of the three cartesian axes
    // There can be more than one mirror at the same time
    // MirrorType (1 - charge density -ve of original, equivalent to method of
    // images, 2 - charge density same as original)
    {
      int ix, iy, iz;
      int jx, jy, jz;  // not used at present
      double sx, sy, sz;
      fstatus =
          neBEMGetMirror(prim, &ix, &jx, &sx, &iy, &jy, &sy, &iz, &jz, &sz);
      if (fstatus != 0) {
        neBEMMessage("neBEMReadGeometry - neBEMGetMirror");
        return -1;
      }
      if (jx < 0) jx = 0;
      if (jy < 0) jy = 0;
      if (jz < 0) jz = 0;

      MirrorTypeX[prim] = ix;
      MirrorTypeY[prim] = iy;
      MirrorTypeZ[prim] = iz;
      if (0) {
        printf("For primitive: %d\n", prim);
        printf("\tMirrorTypeX: %d, MirrorTypeY: %d, MirrorTypeZ: %d\n", ix, iy,
               iz);
        printf("\tNOT USED ==> MirrorInX: %d, MirrorInY: %d, MirrorInZ: %d\n",
               jx, jy, jz);
        printf("\tMirrorDistX: %lg, MirrorDistY: %lg, MirrorDistZ: %lg\n", sx,
               sy, sz);
        getchar();
      }
      if (ix > 0) {
        // printf("neBEMReadGeometry: Mirror have been requested.\n");
        MirrorDistXFromOrigin[prim] = sx;  // assumed to pass through the origin
      } else {
        MirrorDistXFromOrigin[prim] = 0.0;  // pass through the origin
      }
      if (iy > 0) {
        // printf("neBEMReadGeometry: Mirror have been requested.\n");
        MirrorDistYFromOrigin[prim] = sy;
      } else {
        MirrorDistYFromOrigin[prim] = 0.0;
      }
      if (iz > 0) {
        // printf("neBEMReadGeometry: Mirror have been requested.\n");
        MirrorDistZFromOrigin[prim] = sz;
      } else {
        MirrorDistZFromOrigin[prim] = 0.0;
      }
    }  // read mirror information

    // Information on bounding planes
    // ixmin=0: lower x-plane does not exist
    // ixmin=1: lower x-plane does exist
    // cxmin: coordinate of lower x-plane
    // vxmin: potential of lower x-plane
    // Similar for ixmax, iymin, iymax, izmin, izmax
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    double cxmin, cxmax, cymin, cymax, czmin, czmax;
    double vxmin, vxmax, vymin, vymax, vzmin, vzmax;
    fstatus = neBEMGetBoundingPlanes(
        &ixmin, &cxmin, &vxmin, &ixmax, &cxmax, &vxmax, &iymin, &cymin, &vymin,
        &iymax, &cymax, &vymax, &izmin, &czmin, &vzmin, &izmax, &czmax, &vzmax);
    if (fstatus != 0) {
      neBEMMessage("neBEMReadGeometry - neBEMGetBoundingPlanes");
      return -1;
    }
    BndPlaneInXMin[prim] = ixmin;
    BndPlaneInXMax[prim] = ixmax;
    BndPlaneInYMin[prim] = iymin;
    BndPlaneInYMax[prim] = iymax;
    BndPlaneInZMin[prim] = izmin;
    BndPlaneInZMax[prim] = izmax;
    if (ixmin) {
      XBndPlaneInXMin[prim] = cxmin;
      VBndPlaneInXMin[prim] = vxmin;
    } else {
      XBndPlaneInXMin[prim] = 0.0;
      VBndPlaneInXMin[prim] = 0.0;
    }
    if (ixmax > 0) {
      XBndPlaneInXMax[prim] = cxmax;
      VBndPlaneInXMax[prim] = vxmax;
    } else {
      XBndPlaneInXMax[prim] = 0.0;
      VBndPlaneInXMax[prim] = 0.0;
    }
    if (iymin > 0) {
      YBndPlaneInYMin[prim] = cymin;
      VBndPlaneInYMin[prim] = vymin;
    } else {
      YBndPlaneInYMin[prim] = 0.0;
      VBndPlaneInYMin[prim] = 0.0;
    }
    if (iymax > 0) {
      YBndPlaneInYMax[prim] = cymax;
      VBndPlaneInYMax[prim] = vymax;
    } else {
      YBndPlaneInYMax[prim] = 0.0;
      VBndPlaneInYMax[prim] = 0.0;
    }
    if (izmin > 0) {
      ZBndPlaneInZMin[prim] = czmin;
      VBndPlaneInZMin[prim] = vzmin;
    } else {
      ZBndPlaneInZMin[prim] = 0.0;
      VBndPlaneInZMin[prim] = 0.0;
    }
    if (izmax > 0) {
      ZBndPlaneInZMax[prim] = czmax;
      VBndPlaneInZMax[prim] = vzmax;
    } else {
      ZBndPlaneInZMax[prim] = 0.0;
      VBndPlaneInZMax[prim] = 0.0;
    }
  }  // loop over the primitives - major loop

  VolMax = volmax;              // Maximum nb of volumes in the problem
  volRef = ivector(0, VolMax);  // variables to store volume related infomration
  volShape = ivector(0, VolMax);
  volMaterial = ivector(0, VolMax);
  volEpsilon = dvector(0, VolMax);
  volPotential = dvector(0, VolMax);
  volCharge = dvector(0, VolMax);
  volBoundaryType = ivector(0, VolMax);
  for (int volref = 0; volref <= VolMax; ++volref) {
    neBEMVolumeDescription(volref, &volShape[volref], &volMaterial[volref],
                           &volEpsilon[volref], &volPotential[volref],
                           &volCharge[volref], &volBoundaryType[volref]);
    if (dbgFn) {
      printf("volref: %d\n", volref);
      printf("shape: %d,  material: %d\n", volShape[volref],
             volMaterial[volref]);
      printf("eps: %lg,  pot: %lg\n", volEpsilon[volref], volPotential[volref]);
      printf("q: %lg,  type: %d\n", volCharge[volref], volBoundaryType[volref]);
    }
  }

  // Ignore unnecessary primitives from the final count
  // Ideally, all the removal conditions for a primitive should be checked in
  // one loop and the list should be updated in one single go.
  {
    for (int prim = 1; prim <= NbPrimitives; ++prim) {
      OrgnlToEffPrim[prim][0] = prim;
      OrgnlToEffPrim[prim][1] = prim;
      OrgnlToEffPrim[prim][2] = prim;
    }

    {  // Skip primitive
      // Remove skipped primitives having InterfaceType == 0.
      // Also remove primitives having too small dimensions.
      int NbSkipped = 0, effprim;
      double DVertex[4], minDVertex = 0.0;  // maximum number of vertices is 4
      for (int prim = 1; prim <= NbPrimitives; ++prim) {
        effprim = prim - NbSkipped;

        // Check dimensions of the primitive
        for (int vert = 0; vert < NbVertices[prim] - 1; ++vert) {
          DVertex[vert] =
              sqrt(((XVertex[prim][vert + 1] - XVertex[prim][vert]) *
                    (XVertex[prim][vert + 1] - XVertex[prim][vert])) +
                   ((YVertex[prim][vert + 1] - YVertex[prim][vert]) *
                    (YVertex[prim][vert + 1] - YVertex[prim][vert])) +
                   ((ZVertex[prim][vert + 1] - ZVertex[prim][vert]) *
                    (ZVertex[prim][vert + 1] - ZVertex[prim][vert])));
          if (vert == 0)
            minDVertex = DVertex[vert];
          else {
            if (DVertex[vert] < minDVertex) minDVertex = DVertex[vert];
          }
        }

        if ((InterfaceType[prim]) && (minDVertex > MINDIST)) {
          OrgnlToEffPrim[prim][1] = effprim;
          OrgnlToEffPrim[prim][2] = effprim;
          PrimType[effprim] = PrimType[prim];
          NbVertices[effprim] = NbVertices[prim];
          for (int vert = 0; vert < NbVertices[effprim]; ++vert) {
            XVertex[effprim][vert] = XVertex[prim][vert];
            YVertex[effprim][vert] = YVertex[prim][vert];
            ZVertex[effprim][vert] = ZVertex[prim][vert];
          }
          if (PrimType[effprim] == 2)  // wire
          {
            XNorm[effprim] = 0.0;  // modulus not 1 - an absurd trio!
            YNorm[effprim] = 0.0;
            ZNorm[effprim] = 0.0;
            Radius[effprim] = Radius[prim];
          }
          if ((PrimType[effprim] == 3) || (PrimType[effprim] == 4)) {
            XNorm[effprim] = XNorm[prim];
            YNorm[effprim] = YNorm[prim];
            ZNorm[effprim] = ZNorm[prim];
            Radius[effprim] = 0.0;  // absurd radius!
          }
          VolRef1[effprim] = VolRef1[prim];
          VolRef2[effprim] = VolRef2[prim];

          InterfaceType[effprim] = InterfaceType[prim];
          Epsilon1[effprim] = Epsilon1[prim];
          Epsilon2[effprim] = Epsilon2[prim];
          Lambda[effprim] = Lambda[prim];
          ApplPot[effprim] = ApplPot[prim];
          ApplCh[effprim] = ApplCh[prim];
          PeriodicTypeX[effprim] = PeriodicTypeX[prim];
          PeriodicTypeY[effprim] = PeriodicTypeY[prim];
          PeriodicTypeZ[effprim] = PeriodicTypeZ[prim];
          PeriodicInX[effprim] = PeriodicInX[prim];
          PeriodicInY[effprim] = PeriodicInY[prim];
          PeriodicInZ[effprim] = PeriodicInZ[prim];
          XPeriod[effprim] = XPeriod[prim];
          YPeriod[effprim] = YPeriod[prim];
          ZPeriod[effprim] = ZPeriod[prim];
          MirrorTypeX[effprim] = MirrorTypeX[prim];
          MirrorTypeY[effprim] = MirrorTypeY[prim];
          MirrorTypeZ[effprim] = MirrorTypeZ[prim];
          MirrorDistXFromOrigin[effprim] = MirrorDistXFromOrigin[prim];
          MirrorDistYFromOrigin[effprim] = MirrorDistYFromOrigin[prim];
          MirrorDistZFromOrigin[effprim] = MirrorDistZFromOrigin[prim];
          BndPlaneInXMin[effprim] = BndPlaneInXMin[prim];
          BndPlaneInXMax[effprim] = BndPlaneInXMax[prim];
          BndPlaneInYMin[effprim] = BndPlaneInYMin[prim];
          BndPlaneInYMax[effprim] = BndPlaneInYMax[prim];
          BndPlaneInZMin[effprim] = BndPlaneInZMin[prim];
          BndPlaneInZMax[effprim] = BndPlaneInZMax[prim];
          XBndPlaneInXMin[effprim] = XBndPlaneInXMin[prim];
          XBndPlaneInXMax[effprim] = XBndPlaneInXMax[prim];
          YBndPlaneInYMin[effprim] = YBndPlaneInYMin[prim];
          YBndPlaneInYMax[effprim] = YBndPlaneInYMax[prim];
          ZBndPlaneInZMin[effprim] = ZBndPlaneInZMin[prim];
          ZBndPlaneInZMax[effprim] = ZBndPlaneInZMax[prim];
          VBndPlaneInXMin[effprim] = VBndPlaneInXMin[prim];
          VBndPlaneInXMax[effprim] = VBndPlaneInXMax[prim];
          VBndPlaneInYMin[effprim] = VBndPlaneInYMin[prim];
          VBndPlaneInYMax[effprim] = VBndPlaneInYMax[prim];
          VBndPlaneInZMin[effprim] = VBndPlaneInZMin[prim];
          VBndPlaneInZMax[effprim] = VBndPlaneInZMax[prim];
        }  // InterfaceType
        else {
          OrgnlToEffPrim[prim][1] = 0;  // removed from the list
          OrgnlToEffPrim[prim][2] = 0;
          ++NbSkipped;
          if (DebugLevel == 101) {
            printf("Skipped primitive %d, InterfaceType: %d, minDVertex: %lg\n",
                   prim, InterfaceType[prim], minDVertex);
          }
        }
      }  // loop over primitives to remove the skipped primitives
      NbPrimitives -= NbSkipped;
      printf("Number of primitives skipped: %d, Effective NbPrimitives: %d\n",
             NbSkipped, NbPrimitives);
    }  // Skip primitives

    /*
    {	// look for overlaps among primitives
    // Also remove primitives so that periodicities do not lead to overlapped
    // (full or partial?) primitives.
    // Check overlap of primitives due to repetition
    // Algorithm:
    // 1. Consider those primitives which have normal along X.
    // 2. Consider their initial positions in X.
    // 3. Find out their positions after one positive repetition in both X and
    Y.
    // and check which of them fall on already existing primitives having normal
    // along X.
    // A. Similarly for one negative repetition.
    // II. Same check for all the primitives having normal along Y will have to
    be
    // carried out.
    // Similar approach could be useful for an arbitrary primitive if we work in
    // the local coordinate system of the primitive and check all other
    primitives
    // in that LCS.
    // Complications are likely to arise due to corner and edge overlaps which
    will
    // happen to all neighbouring primitives on the same plane.
    int Overlap[NbPrimitives+1][NbPrimitives+1];	// this C++ syntax seems to
    work!
    // int **Overlap;
    // Overlap = imatrix(1, NbPrimitives, 1, NbPrimitives);
    for(unsigned int prim = 1; prim <= NbPrimitives; ++prim)
            for(unsigned int chkprim = 1;
                                                                                                                    chkprim <= NbPrimitives; ++chkprim)
                    Overlap[prim][chkprim] = 0;

    for(unsigned int prim = 1; prim <= NbPrimitives; ++prim)
            {
            if(dbgFn)
                    printf("\nNew XNorm[%d]: %lg: ", prim, XNorm[prim]);

            if(fabs(fabs(XNorm[prim]) - 1.0) >= 1.0e-12)	// primitive || to
    YZ plane continue;

            for(unsigned int chkprim = prim+1;
                                                                                                                    chkprim <= NbPrimitives; ++chkprim)
                    {
                    if(dbgFn)
                            printf("XNorm[%d]: %lg, ", chkprim, XNorm[chkprim]);

                    if(fabs(fabs(XNorm[chkprim]) - 1.0) >= 1.0e-12)	//
    primitive || to YZ plane continue;

                    if(fabs(XVertex[prim][0] - XVertex[chkprim][0]) <= 1.0e-12)
                            {	// same plane; check YZ vertices for full /
    partial overlap double smallY = YVertex[prim][0]; double bigY =
    YVertex[prim][0]; double smallZ = ZVertex[prim][0]; double bigZ =
    ZVertex[prim][0]; for(unsigned int vert = 1; vert<=
    NbVertices[prim]; ++vert)
                                    {
                                    if(smallY > YVertex[prim][vert])
                                            smallY = YVertex[prim][vert];
                                    if(bigY < YVertex[prim][vert])
                                            bigY = YVertex[prim][vert];
                                    if(smallZ > ZVertex[prim][vert])
                                            smallZ = ZVertex[prim][vert];
                                    if(bigZ < ZVertex[prim][vert])
                                            bigZ = ZVertex[prim][vert];
                                    }

                            double smallcY = YVertex[chkprim][0];
                            double bigcY = YVertex[chkprim][0];
                            double smallcZ = ZVertex[chkprim][0];
                            double bigcZ = ZVertex[chkprim][0];
                            for(unsigned int vert = 1;
                                                                                                                                    vert<= NbVertices[chkprim]; ++vert)
                                    {
                                    if(smallcY > YVertex[chkprim][vert])
                                            smallcY = YVertex[chkprim][vert];
                                    if(bigcY < YVertex[chkprim][vert])
                                            bigcY = YVertex[chkprim][vert];
                                    if(smallcZ > ZVertex[chkprim][vert])
                                            smallcZ = ZVertex[chkprim][vert];
                                    if(bigcZ < ZVertex[chkprim][vert])
                                            bigcZ = ZVertex[chkprim][vert];
                                    }

                            // Because of the rearrangement, the big and small
    vertices are now
                            // distinct and unambiguous. Comparisons can now be
    made a lot more
                            // easily.
                            // One important assumption is implicit,however.
                            // The edges of the primitive are assumed to
    parallel to the
                            // Y and Z axes. We are assuming rectangular
    primitives as a
                            // consequence.
                            // One to one match of all the vertices
                            // Subsequent checks should be carried out only if
    Overlap is still 0

                            // No overlap
                            if(((smallcY > smallY) && (smallcY > bigY))
                                            && ((bigcY > smallY) && (bigcY >
    bigY))) Overlap[prim][chkprim] = -1;	// no Y overlap
                            if(Overlap[prim][chkprim] != 0) continue;
                            if(((smallcZ > smallZ) && (smallcZ > bigZ))
                                            && ((bigcZ > smallZ) && (bigcZ >
    bigZ))) Overlap[prim][chkprim] = -2;	// no Z overlap
                            if(Overlap[prim][chkprim] != 0) continue;

                            // Check for shared corners
                            if((fabs(bigcY-smallY) <= MINDIST)
                                            && (fabs(bigcZ-smallZ) <= MINDIST))
                                    Overlap[prim][chkprim] = -3;
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((fabs(bigcY-smallY) <= MINDIST)
                                            && (fabs(smallcZ-bigZ) <= MINDIST))
                                    Overlap[prim][chkprim] = -4;
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((fabs(smallcY-bigY) <= MINDIST)
                                            && (fabs(smallcZ-bigZ) <= MINDIST))
                                    Overlap[prim][chkprim] = -5;
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((fabs(smallcY-bigY) <= MINDIST)
                                            && (fabs(bigcZ-bigZ) <= MINDIST))
                                    Overlap[prim][chkprim] = -6;
                            if(Overlap[prim][chkprim] != 0) continue;

                            // Check for completely / partially shared edges
                            if((fabs(bigcY-smallY) <= MINDIST)
                                            && (((smallcZ > smallZ) && (smallcZ
    < bigZ))	// = implies corner
                                            || ((bigcZ > smallZ) && (bigcZ <
    bigZ))))	// = implies corner Overlap[prim][chkprim] = -7;
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((fabs(smallcZ-bigZ) <= MINDIST)
                                            && (((smallcY > smallY) && (smallcY
    < bigY))	// = implies corner
                                            || ((bigcY > smallY) && (bigcY <
    bigY))))	// = implies corner Overlap[prim][chkprim] = -8;
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((fabs(smallcY-bigY) <= MINDIST)
                                            && (((smallcZ > smallZ) && (smallcZ
    < bigZ))	// = implies corner
                                            || ((bigcZ > smallZ) && (bigcZ <
    bigZ))))	// = implies corner Overlap[prim][chkprim] = -9;
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((fabs(bigcZ-smallZ) <= MINDIST)
                                            && (((smallcY > smallY) && (smallcY
    < bigY))	// = implies corner
                                            || ((bigcY > smallY) && (bigcY <
    bigY))))	// = implies corner Overlap[prim][chkprim] = -10;
                            if(Overlap[prim][chkprim] != 0) continue;

                            // One-to-one overlap - case 0
                            if((fabs(smallcY-smallY) <= MINDIST)
                                            && (fabs(bigcY-bigY) <= MINDIST)
                                            && (fabs(smallcZ-smallZ) <= MINDIST)
                                            && (fabs(bigcZ-bigZ) <= MINDIST))
                                    {	// chkprim has a one-to-one overlap with
    prim Overlap[prim][chkprim] = 1;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;

                            // checkprim is completely contained within prim -
    case 1 if((smallcY > smallY) && (smallcY < bigY)
                                             && (smallcZ > smallZ) && (smallcZ <
    bigZ)
                                            && (bigcY > smallY) && (bigcY <
    bigY)
                                             && (bigcZ > smallZ) && (bigcZ <
    bigZ)) {	// chkprim completely contained within prim
                                    Overlap[prim][chkprim] = 2;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;

                            // checkprim completely contains prim - case 2
                            if((smallY > smallcY) && (smallY < bigcY)
                                             && (smallZ > smallcZ) && (smallZ <
    bigcZ)
                                            && (bigY > smallcY) && (bigY <
    bigcY)
                                             && (bigZ > smallcZ) && (bigZ <
    bigcZ)) {	// chkprim completely contains prim Overlap[prim][chkprim] = 3;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;

                            // checkprim completely contained within prim in Z
    direction - case 3 if((smallcZ > smallZ) && (bigcZ < bigZ))
                                    {
                                    if((smallcY > smallY) && (smallcY < bigY)
                                                    && (bigcY > smallY) &&
    (bigcY > bigY)) {	// but not Y: bigY outside prim Overlap[prim][chkprim] =
    4;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    if((smallcY < smallY) && (smallcY < bigY)
                                                    && (bigcY > smallY) &&
    (bigcY < bigY)) {	// but not Y: smallY outside prim Overlap[prim][chkprim]
    = 5;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    if((smallcY < smallY) && (smallcY < bigY)
                                                    && (bigcY > smallY) &&
    (bigcY > bigY)) {	// but not Y: chkprim wider than prim
                                            Overlap[prim][chkprim] = 6;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    }

                            // checkprim completely contained within prim in Y
    direction - case 4 if((smallcY > smallY) && (bigcY < bigY))
                                    {
                                    if((smallcZ > smallZ) && (smallcZ < bigZ)
    // bigcZ outside
                                                    && (bigcZ > smallZ) &&
    (bigcZ > bigZ)) {	// but not Z: bigZ outside prim Overlap[prim][chkprim] =
    7;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    if((smallcZ < smallZ) && (smallcZ < bigZ)
    // smallcZ outside
                                                    && (bigcZ > smallZ) &&
    (bigcZ < bigZ)) {	// but not Z: smallZ outside prim Overlap[prim][chkprim]
    = 8;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    if((smallcZ < smallZ) && (smallcZ < bigZ)
    // smallcZ outside
                                                    && (bigcZ > smallZ) &&
    (bigcZ > bigZ)) // bigcZ outside {	// but not Z: chkprim taller than prim
                                            Overlap[prim][chkprim] = 9;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    }

                            // checkprim completely contains prim in Z direction
    - case 5 if((smallcZ < smallZ) && (bigcZ > bigZ))
                                    {
                                    if((smallY > smallcY) && (smallY < bigcY)
    // bigY outside
                                                    && (bigY > smallcY) && (bigY
    > bigcY)) {	// but not Y: bigY outside chkprim Overlap[prim][chkprim] = 10;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    if((smallY < smallcY) && (smallY < bigcY)
    // smallY outside
                                                    && (bigY > smallcY) && (bigY
    < bigcY)) {	// but not Y: smallY outside chkprim Overlap[prim][chkprim] =
    11;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    if((smallY < smallcY) && (bigY > bigcY))
    // outside chkprim extent {	// but not Y: prim wider than chkprim
                                            Overlap[prim][chkprim] = 12;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    }

                            // checkprim completely contains prim in Y direction
    - case 6 if((smallcY < smallY) && (bigcY > bigY))
                                    {
                                    if((smallZ > smallcZ) && (smallZ < bigcZ)
    // bigZ outside
                                                    && (bigZ > smallcZ) && (bigZ
    > bigcZ)) {	// but not Z: bigZ outside chkprim Overlap[prim][chkprim] = 13;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    if((smallZ < smallcZ) && (smallZ < bigcZ)
    // smallZ outside
                                                    && (bigZ > smallcZ) && (bigZ
    < bigcZ)) {	// but not Z: smallZ outside chkprim Overlap[prim][chkprim] =
    14;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    if((smallZ < smallcZ) && (smallZ < bigcZ)
    // smallZ and bigZ outside
                                                    && (bigZ > smallcZ) && (bigZ
    < bigcZ)) {	// but not Z: prim taller than chkprim Overlap[prim][chkprim] =
    14;
                                            }
                                    if(Overlap[prim][chkprim] != 0) continue;
                                    }

                            // check cases where only one corner of chkprim is
    within prim
                            // all the other three being beyond.
                            if((smallcY > smallY) && (smallcY < bigY)
                                            && (smallcZ > smallZ) && (smallcZ <
    bigZ)
                                            && (bigcY > bigY) && (bigcZ > bigZ))
                                    {	// only smallcY,smallcZ corner within
    prim Overlap[prim][chkprim] = 15;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((bigcY > smallY) && (bigcY < bigY)
                                            && (smallcZ > smallZ) && (smallcZ <
    bigZ)
                                            && (smallcY < smallY) && (bigcZ >
    bigZ)) {	// only bigcY,smallcZ corner within prim Overlap[prim][chkprim]
    = 16;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((smallcY > smallY) && (smallcY < bigY)
                                            && (bigcZ > smallZ) && (bigcZ <
    bigZ)
                                            && (bigcY > bigY) && (smallcZ <
    smallZ)) {	// only smallcY,bigcZ corner within prim Overlap[prim][chkprim]
    = 17;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((bigcY > smallY) && (bigcY < bigY)
                                            && (bigcZ > smallZ) && (bigcZ <
    bigZ)
                                            && (smallcY < smallY) && (smallcZ <
    smallZ)) {	// only bigcY,bigcZ corner within prim Overlap[prim][chkprim] =
    18;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;

                            // check cases where only one corner of prim is
    within chkprim
                            // all the other three being beyond.
                            if((smallY > smallcY) && (smallY < bigcY)
                                            && (smallZ > smallcZ) && (smallZ <
    bigcZ)
                                            && (bigcY < bigY) && (bigcZ < bigZ))
                                    {	// only smallY,smallZ corner within prim
                                    Overlap[prim][chkprim] = 19;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((bigY > smallcY) && (bigY < bigcY)
                                            && (smallZ > smallcZ) && (smallZ <
    bigcZ)
                                            && (smallcY > smallY) && (bigcZ <
    bigZ)) {	// only bigY,smallZ corner within prim Overlap[prim][chkprim] =
    20;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((smallY > smallcY) && (smallY < bigcY)
                                            && (bigZ > smallcZ) && (bigZ <
    bigcZ)
                                            && (bigcY < bigY) && (smallcZ >
    smallZ)) {	// only smallY,bigZ corner within prim Overlap[prim][chkprim] =
    21;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;
                            if((bigY > smallcY) && (bigY < bigcY)
                                            && (bigZ > smallcZ) && (bigZ <
    bigcZ)
                                            && (smallcY > smallY) && (smallcZ >
    smallZ)) {	// only bigY,bigZ corner within prim Overlap[prim][chkprim] =
    22;
                                    }
                            if(Overlap[prim][chkprim] != 0) continue;
                            }	// same plane; check overlap

                    }	// loop over check primitives
            }	// loop over primitives

    if(dbgFn) printf("\n");
    for(unsigned int prim = 1; prim <= NbPrimitives; ++prim)
            for(unsigned int chkprim = prim+1;
                                                                                                                    chkprim <= NbPrimitives; ++chkprim)
                    {
                    printf("prim: %d, XNorm: %lg, chkprim: %d, XNorm: %lg,
    Overlap: %d\n", prim, XNorm[prim], chkprim, XNorm[chkprim],
    Overlap[prim][chkprim]);
                    }
    }	// look for overlapped primitives
    */

    {  // remove primitives, as specified in a user supplied input file
      int NbRmPrims;
      FILE *rmprimFile = fopen("neBEMRmPrim.inp", "r");
      if (rmprimFile == NULL) {
        printf("neBEMRmPrim.inp absent ... assuming defaults ...\n");
        NbRmPrims = 0;
      } else {
        fscanf(rmprimFile, "NbRmPrims: %d\n", &NbRmPrims);
        if (NbRmPrims) {
          int tint;
#ifdef __cplusplus
          std::vector<double> rmXNorm(NbRmPrims + 1, 0.);
          std::vector<double> rmYNorm(NbRmPrims + 1, 0.);
          std::vector<double> rmZNorm(NbRmPrims + 1, 0.);
          std::vector<double> rmXVert(NbRmPrims + 1, 0.);
          std::vector<double> rmYVert(NbRmPrims + 1, 0.);
          std::vector<double> rmZVert(NbRmPrims + 1, 0.);
#else
          double rmXNorm[NbRmPrims + 1], rmYNorm[NbRmPrims + 1];
          double rmZNorm[NbRmPrims + 1];
          double rmXVert[NbRmPrims + 1], rmYVert[NbRmPrims + 1];
          double rmZVert[NbRmPrims + 1];
#endif
          for (int rmprim = 1; rmprim <= NbRmPrims; ++rmprim) {
            fscanf(rmprimFile, "Prim: %d\n", &tint);
            fscanf(rmprimFile, "rmXNorm: %le\n", &rmXNorm[rmprim]);
            fscanf(rmprimFile, "rmYNorm: %le\n", &rmYNorm[rmprim]);
            fscanf(rmprimFile, "rmZNorm: %le\n", &rmZNorm[rmprim]);
            fscanf(rmprimFile, "rmXVert: %le\n", &rmXVert[rmprim]);
            fscanf(rmprimFile, "rmYVert: %le\n", &rmYVert[rmprim]);
            fscanf(rmprimFile, "rmZVert: %le\n", &rmZVert[rmprim]);
            printf(
                "rmprim: %d, rmXNorm: %lg, rmYNorm: %lg, rmZNorm: %lg, "
                "rmXVert: %lg, rmYVert: %lg, rmZVert: %lg\n",
                rmprim, rmXNorm[rmprim], rmYNorm[rmprim], rmZNorm[rmprim],
                rmXVert[rmprim], rmYVert[rmprim], rmZVert[rmprim]);
          }
#ifdef __cplusplus
          std::vector<int> remove(NbPrimitives + 1, 0);
#else
          int remove[NbPrimitives + 1];
#endif
          // Check updated prim list
          for (int prim = 1; prim <= NbPrimitives; ++prim) {
            remove[prim] = 0;
            if (dbgFn) {
              printf("\n\nprim: %d, XVertex: %lg, YVertex: %lg, ZVertex: %lg\n",
                     prim, XVertex[prim][0], YVertex[prim][0],
                     ZVertex[prim][0]);
              printf("XNorm: %lg, YNorm: %lg, ZNorm: %lg\n", XNorm[prim],
                     YNorm[prim], ZNorm[prim]);
            }

            for (int rmprim = 1; rmprim <= NbRmPrims; ++rmprim) {
              if (dbgFn) {
                printf(
                    "rmprim: %d, rmXVertex: %lg, rmYVertex: %lg, rmZVertex: "
                    "%lg\n",
                    rmprim, rmXVert[rmprim], rmYVert[rmprim], rmZVert[rmprim]);
                printf("rmXNorm: %lg, rmYNorm: %lg, rmZNorm: %lg\n",
                       rmXNorm[rmprim], rmYNorm[rmprim], rmZNorm[rmprim]);
              }

              // check the normal
              if ((fabs(fabs(XNorm[prim]) - fabs(rmXNorm[rmprim])) <=
                   MINDIST) &&
                  (fabs(fabs(YNorm[prim]) - fabs(rmYNorm[rmprim])) <=
                   MINDIST) &&
                  (fabs(fabs(ZNorm[prim]) - fabs(rmZNorm[rmprim])) <=
                   MINDIST)) {  // prim and rmprim are parallel
                // coplanarity check to be implemented later.
                // For the time-being, we will assume that the planes to be
                // removed have their normals parallel to a given axis. So, we
                // only check that and remove the primitive if the distace along
                // that axis match. Possible pitfall => the primitives may be
                // coplanar but non-overlapping!
                if (fabs(fabs(XNorm[prim]) - 1.0) <= 1.0e-12) {
                  // primitive || to YZ
                  if (fabs(XVertex[prim][0] - rmXVert[rmprim]) <= MINDIST) {
                    remove[prim] = 1;
                  }
                }
                if (fabs(fabs(YNorm[prim]) - 1.0) <= 1.0e-12) {
                  // primitive || to XZ
                  if (fabs(YVertex[prim][0] - rmYVert[rmprim]) <= MINDIST) {
                    remove[prim] = 1;
                  }
                }
                if (fabs(fabs(ZNorm[prim]) - 1.0) <= 1.0e-12) {
                  // primitive || to XY
                  if (fabs(ZVertex[prim][0] - rmZVert[rmprim]) <= MINDIST) {
                    remove[prim] = 1;
                  }
                }
              }  // case where prim and rmprim are parallel
              if (dbgFn) {
                printf("prim: %d, rmprim: %d, remove: %d\n", prim, rmprim,
                       remove[prim]);
              }
              if (remove[prim] == 1)
                break;  // once removed, no point checking others
            }           // for rmprim - loop over all removal specification

          }  // for prim loop over all primitives

          int NbRemoved = 0;
          FILE *fprrm = fopen("RmPrims.info", "w");
          if (fprrm == NULL) {
            printf(
                "error opening RmPrims.info file in write mode ... "
                "returning\n");
            return (-1);
          }
          // Note that some of the original primitives have already been removed
          // based on interface and dimension considerations
          int orgnlNb = 0;
          for (int prim = 1; prim <= NbPrimitives; ++prim) {
            // identify primitive number in the original list
            for (int orgnlprim = 1; orgnlprim <= OrgnlNbPrimitives;
                 ++orgnlprim) {
              if (OrgnlToEffPrim[orgnlprim][1] ==
                  prim)  // number updated for intrfc
              {
                orgnlNb = orgnlprim;
                break;
              }
            }  // loop for finding out its position in the original list

            if (remove[prim] == 1) {
              ++NbRemoved;
              OrgnlToEffPrim[orgnlNb][2] = 0;
              fprintf(fprrm, "NbRemoved: %d, Removed primitive: %d\n",
                      NbRemoved, prim);
              fprintf(fprrm, "PrimType: %d\n", PrimType[prim]);
              fprintf(fprrm, "NbVertices: %d\n", NbVertices[prim]);
              for (int vert = 0; vert < NbVertices[prim]; ++vert) {
                fprintf(fprrm, "Vertx %d: %lg, %lg, %lg\n", vert,
                        XVertex[prim][vert], YVertex[prim][vert],
                        ZVertex[prim][vert]);
              }
              fprintf(fprrm, "Normals: %lg, %lg, %lg\n", XNorm[prim],
                      YNorm[prim], ZNorm[prim]);
              continue;
            }       // if remove
            else {  // keep this one in the updated list of primitives
              int effprim = prim - NbRemoved;

              OrgnlToEffPrim[orgnlNb][2] = effprim;
              PrimType[effprim] = PrimType[prim];
              NbVertices[effprim] = NbVertices[prim];
              for (int vert = 0; vert < NbVertices[effprim]; ++vert) {
                XVertex[effprim][vert] = XVertex[prim][vert];
                YVertex[effprim][vert] = YVertex[prim][vert];
                ZVertex[effprim][vert] = ZVertex[prim][vert];
              }
              if (PrimType[effprim] == 2) {
                // wire
                XNorm[effprim] = 0.0;  // modulus not 1 - an absurd trio!
                YNorm[effprim] = 0.0;
                ZNorm[effprim] = 0.0;
                Radius[effprim] = Radius[prim];
              }
              if ((PrimType[effprim] == 3) || (PrimType[effprim] == 4)) {
                XNorm[effprim] = XNorm[prim];
                YNorm[effprim] = YNorm[prim];
                ZNorm[effprim] = ZNorm[prim];
                Radius[effprim] = 0.0;  // absurd radius!
              }
              VolRef1[effprim] = VolRef1[prim];
              VolRef2[effprim] = VolRef2[prim];

              InterfaceType[effprim] = InterfaceType[prim];
              Epsilon1[effprim] = Epsilon1[prim];
              Epsilon2[effprim] = Epsilon2[prim];
              Lambda[effprim] = Lambda[prim];
              ApplPot[effprim] = ApplPot[prim];
              ApplCh[effprim] = ApplCh[prim];
              PeriodicTypeX[effprim] = PeriodicTypeX[prim];
              PeriodicTypeY[effprim] = PeriodicTypeY[prim];
              PeriodicTypeZ[effprim] = PeriodicTypeZ[prim];
              PeriodicInX[effprim] = PeriodicInX[prim];
              PeriodicInY[effprim] = PeriodicInY[prim];
              PeriodicInZ[effprim] = PeriodicInZ[prim];
              XPeriod[effprim] = XPeriod[prim];
              YPeriod[effprim] = YPeriod[prim];
              ZPeriod[effprim] = ZPeriod[prim];
              MirrorTypeX[effprim] = MirrorTypeX[prim];
              MirrorTypeY[effprim] = MirrorTypeY[prim];
              MirrorTypeZ[effprim] = MirrorTypeZ[prim];
              MirrorDistXFromOrigin[effprim] = MirrorDistXFromOrigin[prim];
              MirrorDistYFromOrigin[effprim] = MirrorDistYFromOrigin[prim];
              MirrorDistZFromOrigin[effprim] = MirrorDistZFromOrigin[prim];
              BndPlaneInXMin[effprim] = BndPlaneInXMin[prim];
              BndPlaneInXMax[effprim] = BndPlaneInXMax[prim];
              BndPlaneInYMin[effprim] = BndPlaneInYMin[prim];
              BndPlaneInYMax[effprim] = BndPlaneInYMax[prim];
              BndPlaneInZMin[effprim] = BndPlaneInZMin[prim];
              BndPlaneInZMax[effprim] = BndPlaneInZMax[prim];
              XBndPlaneInXMin[effprim] = XBndPlaneInXMin[prim];
              XBndPlaneInXMax[effprim] = XBndPlaneInXMax[prim];
              YBndPlaneInYMin[effprim] = YBndPlaneInYMin[prim];
              YBndPlaneInYMax[effprim] = YBndPlaneInYMax[prim];
              ZBndPlaneInZMin[effprim] = ZBndPlaneInZMin[prim];
              ZBndPlaneInZMax[effprim] = ZBndPlaneInZMax[prim];
              VBndPlaneInXMin[effprim] = VBndPlaneInXMin[prim];
              VBndPlaneInXMax[effprim] = VBndPlaneInXMax[prim];
              VBndPlaneInYMin[effprim] = VBndPlaneInYMin[prim];
              VBndPlaneInYMax[effprim] = VBndPlaneInYMax[prim];
              VBndPlaneInZMin[effprim] = VBndPlaneInZMin[prim];
              VBndPlaneInZMax[effprim] = VBndPlaneInZMax[prim];
            }  // else remove == 0
          }    // loop over primitives to remove the primitives tagged to be
               // removed
          fclose(fprrm);

          NbPrimitives -= NbRemoved;
          printf(
              "Number of primitives removed: %d, Effective NbPrimitives: %d\n",
              NbRemoved, NbPrimitives);
          fflush(stdout);
        }  // if NbRmPrimes true, implying primitives need to be removed
        fclose(rmprimFile);
      }  // if the rmprimFile is not NULL, prepare to remove primitives
    }    // remove primitives as desired by the user

    FILE *fignore = fopen("IgnorePrims.info", "w");
    if (fignore == NULL) {
      printf(
          "error opening IgnorePrims.info file in write mode ... returning\n");
      return (-1);
    }

    for (int prim = 1; prim <= OrgnlNbPrimitives; ++prim) {
      fprintf(fignore, "%d %d %d\n", OrgnlToEffPrim[prim][0],
              OrgnlToEffPrim[prim][1], OrgnlToEffPrim[prim][2]);
    }

    fclose(fignore);
  }  // Ignore unnecessary primitives from the final count

  // Reduced-Order Modelling information
  printf("ROM: switch to primitive representation after %d repetitions.\n",
         PrimAfter);

  // Store model data in native neBEM format
  char NativeInFile[256];

  strcpy(NativeInFile, NativeOutDir);
  strcat(NativeInFile, "/neBEMNative.inp");
  FILE *fNativeInFile = fopen(NativeInFile, "w");
  fprintf(fNativeInFile, "#====>Input directory\n");
  fprintf(fNativeInFile, "%s\n", NativePrimDir);
  fprintf(fNativeInFile, "#====>No. of primitives:\n");
  fprintf(fNativeInFile, "%d\n", NbPrimitives);
  fprintf(fNativeInFile, "#====>No. of volumes:\n");
  fprintf(fNativeInFile, "%d\n", VolMax);
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    char NativePrimFile[256];
    char strPrimNb[11];
    sprintf(strPrimNb, "%d", prim);
    strcpy(NativePrimFile, "Primitive");
    strcat(NativePrimFile, strPrimNb);
    strcat(NativePrimFile, ".inp");

    fprintf(fNativeInFile, "#Input filename for primitive:\n");
    fprintf(fNativeInFile, "%s\n", NativePrimFile);

    strcpy(NativePrimFile, NativePrimDir);
    strcat(NativePrimFile, "/Primitive");
    strcat(NativePrimFile, strPrimNb);
    strcat(NativePrimFile, ".inp");

    FILE *fNativePrim = fopen(NativePrimFile, "w");

    fprintf(fNativePrim, "#Number of vertices:\n");
    fprintf(fNativePrim, "%d\n", NbVertices[prim]);
    for (int vertex = 0; vertex < NbVertices[prim]; ++vertex) {
      fprintf(fNativePrim, "#Vertex %d:\n", vertex);
      fprintf(fNativePrim, "%lg %lg %lg\n", XVertex[prim][vertex],
              YVertex[prim][vertex], ZVertex[prim][vertex]);
    }  // for vertex
    fprintf(fNativePrim, "#Outward normal:\n");
    fprintf(fNativePrim, "%lg %lg %lg\n", XNorm[prim], YNorm[prim],
            ZNorm[prim]);
    fprintf(fNativePrim, "#Volume references:\n");
    fprintf(fNativePrim, "%d %d\n", VolRef1[prim], VolRef2[prim]);
    fprintf(fNativePrim, "#Maximum number of segments:\n");
    fprintf(fNativePrim, "%d %d\n", 0, 0);  // use the trio target, min, max
    fprintf(fNativePrim, "#Information on X periodicity:\n");
    fprintf(fNativePrim, "%d %d %lg\n", PeriodicTypeX[prim], PeriodicInX[prim],
            XPeriod[prim]);
    fprintf(fNativePrim, "#Information on Y periodicity:\n");
    fprintf(fNativePrim, "%d %d %lg\n", PeriodicTypeY[prim], PeriodicInY[prim],
            YPeriod[prim]);
    fprintf(fNativePrim, "#Information on Z periodicity:\n");
    fprintf(fNativePrim, "%d %d %lg\n", PeriodicTypeZ[prim], PeriodicInZ[prim],
            ZPeriod[prim]);

    fclose(fNativePrim);
  }  // for prim

  fprintf(fNativeInFile, "#====>No. of boundary conditions per element:\n");
  fprintf(fNativeInFile, "%d\n", 1);  // CHECK!!!
  fprintf(fNativeInFile, "#====>Device output directory name:\n");
  fprintf(fNativeInFile, "NativeResults\n");
  fprintf(fNativeInFile, "#====>Map input file:\n");
  fprintf(fNativeInFile, "MapFile.inp\n");
  fclose(fNativeInFile);

  for (int volume = 0; volume <= VolMax; ++volume) {
    // Note that materials from 1 to 10 are conductors and
    //                     from 11 to 20 are dielectrics
    int shape, material, boundarytype;
    double epsilon, potential, charge;
    neBEMVolumeDescription(volume, &shape, &material, &epsilon, &potential,
                           &charge, &boundarytype);

    char NativeVolFile[256];
    char strVolNb[11];
    sprintf(strVolNb, "%d", volume);

    strcpy(NativeVolFile, NativePrimDir);
    strcat(NativeVolFile, "/Volume");
    strcat(NativeVolFile, strVolNb);
    strcat(NativeVolFile, ".inp");

    FILE *fNativeVol = fopen(NativeVolFile, "w");

    fprintf(fNativeVol, "#Shape of the volume:\n");
    fprintf(fNativeVol, "%d\n", shape);
    fprintf(fNativeVol, "#Material of the volume:\n");
    fprintf(fNativeVol, "%d\n", material);
    fprintf(fNativeVol, "#Relative permittivity:\n");
    fprintf(fNativeVol, "%lg\n", epsilon);
    fprintf(fNativeVol, "#Potential:\n");
    fprintf(fNativeVol, "%lg\n", potential);
    fprintf(fNativeVol, "#Charge:\n");
    fprintf(fNativeVol, "%lg\n", charge);
    fprintf(fNativeVol, "#Boundary type:\n");
    fprintf(fNativeVol, "%d\n", boundarytype);

    fclose(fNativeVol);
  }

  // create a dummy map file
  {
    char NativeMapFile[256];

    strcpy(NativeMapFile, NativePrimDir);
    strcat(NativeMapFile, "/MapFile.inp");

    FILE *fNativeMap = fopen(NativeMapFile, "w");

    fprintf(fNativeMap, "#Number of maps\n");
    fprintf(fNativeMap, "%d\n", 0);
    fprintf(fNativeMap, "#Number of lines\n");
    fprintf(fNativeMap, "%d\n", 0);
    fprintf(fNativeMap, "#Number of points\n");
    fprintf(fNativeMap, "%d\n", 0);

    fclose(fNativeMap);
  }

  // Store primitive related data in a file for a new model, if opted for
  if (NewModel && OptStorePrimitives) {
    if (OptFormattedFile) {
      fstatus = WritePrimitives();
      if (fstatus) {
        neBEMMessage("neBEMReadGeometry - problem writing Primtives.\n");
        return -1;
      }
    }  // formatted file

    if (OptUnformattedFile) {
      neBEMMessage(
          "neBEMReadGeometry - unformatted write not inplemented yet.\n");
      return -1;
    }  // unformatted file
  }    // store primitives

  printf("neBEMReadGeometry: Geometry read!\n");

  neBEMState = 3;  // info about primitives read in after initialization and Nbs

  stopClock = clock();
  neBEMTimeElapsed(startClock, stopClock);
  printf("to read geometry\n");

  return (0);  // Success!
}  // neBEMReadGeometry ends

// Discretize the primitives using discretization information supplied by the
// user.
int neBEMDiscretize(int **NbElemsOnPrimitives) {
  int fstatus;

  // For a model and a mesh that were defined before and for which data were
  // stored in two files - one for primitives and one for elements
  // The following operation does not assume that the primitives have been read
  // in from a stored file. In essence, these two read-in-s are maintained
  // independent of each other. However, for a stored element file to be useful,
  // both the model and mesh have to be old (not new).
  if ((!NewModel) && (!NewMesh) && (!NewBC) && (OptStoreElements)) {
    fstatus = ReadElements();
    if (fstatus) {
      neBEMMessage("neBEMDiscretize - problem reading stored Elements.\n");
      return -1;
    }
    neBEMState = 4;
    return 0;
  }

  // Otherwise, continue with fresh discretization
  if (neBEMState != 3) {
    printf("discretization can continue only in State 3 ...\n");
    return -1;
  }

  // Only the number of primitives has been ascertained.
  // All the rest globally important numbers will be determined in this function
  NbSurfs = 0;
  NbWires = 0;
  NbElements = 0;

  // Here, the primitive can be analyzed and the elements necessary to
  // discretize it, may be determined. A user hint may help that can be supplied
  // during the setting up of the device. The hint may be as naive as gross,
  // coarse, medium, regular, fine depending on which the element size (in %
  // of the device / primitive) may be decided upon.
  // Earlier method of specifying the number of primitives on a priori has been
  // over-ridden.
  // Now we prescribe the length / area of each element and discretize each
  // primitive such that each element has an length / area close to the
  // suggested value.
  // Since the number of elements are being decided here, all the shapes will
  // have to be one among wire, right triangle and restangle. Any decomposition
  // of arbitrary polygons into rectangles (as many as possible) and right
  // triangles will have to be done earlier. All the counts have been
  // incremented by one to be on the safe side! The additional memory can be
  // minimized through a more careful computation of the required allocation for
  // each type of primitive.
  char MeshLogFile[256];

  strcpy(MeshLogFile, MeshOutDir);
  strcat(MeshLogFile, "/MeshLog.out");
  fMeshLog = fopen(MeshLogFile, "w");
  fprintf(fMeshLog, "Details of primitive discretization\n");

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    if (NbVertices[prim] == 4) {
      NbSurfSegX[prim] = NbElemsOnPrimitives[prim][1];
      NbSurfSegZ[prim] = NbElemsOnPrimitives[prim][2];
      fstatus = AnalyzePrimitive(prim, &NbSurfSegX[prim], &NbSurfSegZ[prim]);
      if (fstatus == 0) {
        neBEMMessage("neBEMDiscretize - AnalyzePrimitve");
        return -1;
      }
      NbElements += (NbSurfSegX[prim] + 1) * (NbSurfSegZ[prim] + 1);
    }
    if (NbVertices[prim] == 3) {
      NbSurfSegX[prim] = NbElemsOnPrimitives[prim][1];
      NbSurfSegZ[prim] = NbElemsOnPrimitives[prim][2];
      fstatus = AnalyzePrimitive(prim, &NbSurfSegX[prim], &NbSurfSegZ[prim]);
      if (fstatus == 0) {
        neBEMMessage("neBEMDiscretize - AnalyzePrimitive");
        return -1;
      }
      NbElements += (NbSurfSegX[prim] + 1) * (NbSurfSegZ[prim] + 1);
    }
    if (NbVertices[prim] == 2) {
      int itmp;
      NbWireSeg[prim] = NbElemsOnPrimitives[prim][1];
      fstatus = AnalyzePrimitive(prim, &NbWireSeg[prim], &itmp);
      if (fstatus == 0) {
        neBEMMessage("neBEMDiscretize - AnalyzePrimitive");
        return -1;
      }
      NbElements += (NbWireSeg[prim] + 1);
    }
    if (DebugLevel == 101) {
      if (NbVertices[prim] == 2) {
        printf("Primitive %d to be discretized into %d elements.\n", prim,
               NbWireSeg[prim]);
      } else {
        printf("Primitive %d to be discretized into %d X %d elements.\n", prim,
               NbSurfSegX[prim], NbSurfSegZ[prim]);
      }
    }
  }
  printf("Memory allocated for maximum %d elements.\n", NbElements);
  fclose(fMeshLog);

  // Allocate enough space to store all the elements
  if (neBEMState == 3) {
    printf("neBEMDiscretize: NbElements = %d, sizeof(Element) = %zu\n",
           NbElements, sizeof(Element));
    if (EleArr) {
      Element *tmp = (Element *)realloc(EleArr, NbElements * sizeof(Element));
      if (tmp != NULL) {
        EleArr = tmp;
        EleCntr = 0;
      } else {
        free(EleArr);
        printf("neBEMDiscretize: Re-allocating EleArr failed.\n");
        return (1);
      }
      printf("neBEMDiscretize: Re-allocated EleArr.\n");
    }  // if EleArr => re-allocation
    else {
      EleArr = (Element *)malloc(NbElements * sizeof(Element));
      if (EleArr == NULL) {
        neBEMMessage("neBEMDiscretize - EleArr malloc");
        return -1;
      }
    }  // else EleArr => fresh allocation
  }    // neBEMState == 3

  // Prepare a data file that will contain the plotting information of the
  // the primitives and the elements
  if (OptGnuplot) {
    char GnuFile[256];
    strcpy(GnuFile, MeshOutDir);
    strcat(GnuFile, "/GViewDir/gPrimView.gp");
    fgnuPrim = fopen(GnuFile, "w");
    fprintf(fgnuPrim, "set title \"neBEM primitives in gnuplot VIEWER\"\n");
    // fprintf(fgnu, "#set label 1 \'LengthScale = %d\', LengthScale, right\n");
    fprintf(fgnuPrim, "#set pm3d\n");
    fprintf(fgnuPrim, "#set style data pm3d\n");
    fprintf(fgnuPrim, "#set palette model CMY\n");
    fprintf(fgnuPrim, "set hidden3d\n");
    fprintf(fgnuPrim, "set nokey\n");
    fprintf(fgnuPrim, "set xlabel \"X\"\n");
    fprintf(fgnuPrim, "set ylabel \"Y\"\n");
    fprintf(fgnuPrim, "set zlabel \"Z\"\n");
    fprintf(fgnuPrim, "set view 70, 335, 1, 1\n");
    fprintf(fgnuPrim, "\nsplot \\\n");

    strcpy(GnuFile, MeshOutDir);
    strcat(GnuFile, "/GViewDir/gElemView.gp");
    fgnuElem = fopen(GnuFile, "w");
    fprintf(fgnuElem, "set title \"neBEM elements in gnuplot VIEWER\"\n");
    // fprintf(fgnu, "#set label 1 \'LengthScale = %d\', LengthScale, right\n");
    fprintf(fgnuElem, "#set pm3d\n");
    fprintf(fgnuElem, "#set style data pm3d\n");
    fprintf(fgnuElem, "#set palette model CMY\n");
    fprintf(fgnuElem, "set hidden3d\n");
    fprintf(fgnuElem, "set nokey\n");
    fprintf(fgnuElem, "set xlabel \"X\"\n");
    fprintf(fgnuElem, "set ylabel \"Y\"\n");
    fprintf(fgnuElem, "set zlabel \"Z\"\n");
    fprintf(fgnuElem, "set view 70, 335, 1, 1\n");
    fprintf(fgnuElem, "\nsplot \\\n");

    strcpy(GnuFile, MeshOutDir);
    strcat(GnuFile, "/GViewDir/gMeshView.gp");
    fgnuMesh = fopen(GnuFile, "w");
    fprintf(fgnuMesh, "set title \"neBEM mesh in gnuplot VIEWER\"\n");
    // fprintf(fgnu, "#set label 1 \'LengthScale = %d\', LengthScale, right\n");
    fprintf(fgnuMesh, "#set pm3d\n");
    fprintf(fgnuMesh, "#set style data pm3d\n");
    fprintf(fgnuMesh, "#set palette model CMY\n");
    fprintf(fgnuMesh, "set hidden3d\n");
    fprintf(fgnuMesh, "set nokey\n");
    fprintf(fgnuMesh, "set xlabel \"X\"\n");
    fprintf(fgnuMesh, "set ylabel \"Y\"\n");
    fprintf(fgnuMesh, "set zlabel \"Z\"\n");
    fprintf(fgnuMesh, "set view 70, 335, 1, 1\n");
    fprintf(fgnuMesh, "\nsplot \\\n");
  }

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    switch (PrimType[prim]) {
      int fstatus;
      case 3:  // triangular surface
      case 4:  // rectangular surface
        ++NbSurfs;
        fstatus = SurfaceElements(
            prim, NbVertices[prim], XVertex[prim], YVertex[prim], ZVertex[prim],
            XNorm[prim], YNorm[prim], ZNorm[prim], VolRef1[prim], VolRef2[prim],
            InterfaceType[prim], ApplPot[prim], ApplCh[prim], Lambda[prim],
            NbSurfSegX[prim], NbSurfSegZ[prim]);
        if (fstatus != 0) {
          neBEMMessage("neBEMDiscretize - SurfaceElements");
          return -1;
        }
        break;
      case 2:       // wire - a wire presumably has only 2 vertices
        ++NbWires;  // it has one radius and one segmentation information
        fstatus = WireElements(
            prim, NbVertices[prim], XVertex[prim], YVertex[prim], ZVertex[prim],
            Radius[prim], VolRef1[prim], VolRef2[prim], InterfaceType[prim],
            ApplPot[prim], ApplCh[prim], Lambda[prim], NbWireSeg[prim]);
        if (fstatus != 0) {
          neBEMMessage("neBEMDiscretize - WireElements");
          return -1;
        }
        break;

      default:
        printf("PrimType out of range in CreateElements ... exiting ...\n");
        exit(-1);
    }  // switch PrimType ends
  }    // loop on prim number ends

  if (OptGnuplot) {
    fprintf(fgnuPrim, "\n\npause-1");
    fclose(fgnuPrim);
    fprintf(fgnuElem, "\n\npause-1");
    fclose(fgnuElem);
    fprintf(fgnuMesh, "\n\npause-1");
    fclose(fgnuMesh);
  }

  // If the required memory exceeds the maximum allowed number of elements
  if (EleCntr > NbElements) {
    neBEMMessage("neBEMDiscretize - EleCntr more than NbElements!");
    return -1;
  }

  // Check whether collocation points overlap
  {
    int startcntr = 1, cntr1, cntr2, NbCollPtOverlaps = 0;
    Point3D pt1, pt2;
    double dist;
    for (cntr1 = startcntr; cntr1 <= EleCntr; ++cntr1) {
      pt1.X = (EleArr + cntr1 - 1)->BC.CollPt.X;
      pt1.Y = (EleArr + cntr1 - 1)->BC.CollPt.Y;
      pt1.Z = (EleArr + cntr1 - 1)->BC.CollPt.Z;
      for (cntr2 = cntr1 + 1; cntr2 <= EleCntr; ++cntr2) {
        pt2.X = (EleArr + cntr2 - 1)->BC.CollPt.X;
        pt2.Y = (EleArr + cntr2 - 1)->BC.CollPt.Y;
        pt2.Z = (EleArr + cntr2 - 1)->BC.CollPt.Z;

        dist = GetDistancePoint3D(&pt1, &pt2);
        if (dist <=
            MINDIST)  // we need a linked-list here so that the overlapped
        {  // element is easily deleted and the rest upgraded immediately
          ++NbCollPtOverlaps;

          // Upgrade the element array manually, starting from cntr2 and restart
          // the overlap check. At present it is only a warning to the user with
          // some relevant information.
          // Find the primitives and volumes for the overlapping elements
          // The element structure should also maintain information on the
          // volumes that an element belongs to.
          int prim1 = (EleArr + cntr1 - 1)->PrimitiveNb;
          int volele1 = VolRef1[prim1];
          int prim2 = (EleArr + cntr2 - 1)->PrimitiveNb;
          int volele2 = VolRef1[prim2];

          neBEMMessage("neBEMDiscretize - Overlapping collocation points!");
          printf("Element %d, primitive %d, volume %d overlaps with\n", cntr1,
                 prim1, volele1);
          printf("\telement %d, primitive %d, volume %d.\n", cntr2, prim2,
                 volele2);
          printf("\tposition 1: (%g , %g , %g) micron,\n", 1e6 * pt1.X,
                 1e6 * pt1.Y, 1e6 * pt1.Z);
          printf("\tposition 2: (%g , %g , %g) micron.\n", 1e6 * pt2.X,
                 1e6 * pt2.Y, 1e6 * pt2.Z);
          printf("Please redo the geometry.\n");
          return -1;
        }  // if dist <= MINDIST
      }    // for cntr2
    }      // for cntr1
  }        // check collocation point overlap

  NbElements = EleCntr;  // the final number of elements
  printf("Total final number of elements: %d\n", NbElements);

  // Store element related data in a file for a new mesh created, if opted for
  if (NewMesh && OptStoreElements) {
    if (OptFormattedFile) {
      fstatus = WriteElements();
      if (fstatus) {
        neBEMMessage("neBEMDiscretize - problem writing Elements.\n");
        return -1;
      }
    }  // formatted file

    if (OptUnformattedFile) {
      neBEMMessage(
          "neBEMDiscretize - unformatted write not inplemented yet.\n");
      return -1;
    }  // unformatted file
  }    // store elements

  neBEMState = 4;
  stopClock = clock();
  neBEMTimeElapsed(startClock, stopClock);
  printf("to complete discretization\n");

  return (0);
}  // neBEMDiscretize ends

int neBEMBoundaryConditions(void) {
  startClock = clock();

  // The boundary conditions can be set only with neBEMState == 4 or 7
  // This state is assigned either after element discretization has been
  // completed or in a condition when we are looking for modifying only the
  // boundary condition for a device having same geometry (hence, the same
  // inverted influence coefficient matrix)
  if ((neBEMState == 4) || (neBEMState == 7)) {
    int fstatus = BoundaryConditions();
    if (fstatus != 0) {
      neBEMMessage("neBEMBondaryConditions - BoundaryConditions");
      return -1;
    }
    if (neBEMState == 4) neBEMState = 5;  // create LHMatrix, invert etc
    if (neBEMState == 7) neBEMState = 8;  // assume LHS and inversion to be over
  } else {
    printf("Boundary conditions can be set only in state 4 / 7 ...\n");
    printf("returning ...\n");
    return (-1);
  }

  stopClock = clock();
  neBEMTimeElapsed(startClock, stopClock);
  printf("to setup boundary conditions.\n");

  return (0);
}  // neBEMBoundaryConditions ends

// int neBEMKnownCharges(int (*Pt2UserFunction)(void))
int neBEMKnownCharges(void) {
  int debugFn = 0;

  /*
  // How to check that the pointer points to a valid function?
  if(Pt2UserFunction == NULL)
          {
          printf("Not a valid function ... returning ...\n");
          return(-1);
          }
  else
          {
          // printf("Pt2UserFunction points to %p\n", Pt2UserFunction);
          }

  // Status of the known charges conditions is meaningful only after the
  // element discretization has been completed, i.e., beyond the 5th state.
  if(neBEMState >= 5)
          {	// the following function is declared in the Interface.h.
          int fstatus = (*Pt2UserFunction)();	// user to supply function
          if(fstatus != 0)
                  {
                  neBEMMessage("neBEMKnownCharges - Pt2UserFunction");
                  return -1;
                  }
          if(neBEMState > 5)	// assume LHS and inversion to be over?
                  neBEMState = 8;
          }
  else
          {
          printf("Known charges are meaningful only beyond state 4 ...\n");
          printf("returning ...\n");
          return(-1);
          }
  */

  // Set up parameters related to known charge calculations
  // Electrons and ions can be distributed as points, lines, areas or volumes.
  {
    FILE *KnChInpFile = fopen("neBEMKnCh.inp", "r");
    if (KnChInpFile == NULL) {
      OptKnCh = 0;
      printf(
          "neBEMKnCh.inp absent ... assuming absence of known charges ...\n");
    } else {
      fscanf(KnChInpFile, "OptKnCh: %d\n", &OptKnCh);
      if (1) printf("OptKnCh: %d\n", OptKnCh);

      if (!OptKnCh) printf("OptKnCh = 0 ... assuming no known charges ...\n");

      if (OptKnCh) {
        char PointKnChFile[256];
        char LineKnChFile[256];
        char AreaKnChFile[256];
        char VolumeKnChFile[256];
        double KnChFactor;

        fscanf(KnChInpFile, "NbPointsKnCh: %d\n", &NbPointsKnCh);
        fscanf(KnChInpFile, "PointKnChFile: %255s\n", PointKnChFile);
        fscanf(KnChInpFile, "NbLinesKnCh: %d\n", &NbLinesKnCh);
        fscanf(KnChInpFile, "LineKnChFile: %255s\n", LineKnChFile);
        fscanf(KnChInpFile, "NbAreasKnCh: %d\n", &NbAreasKnCh);
        fscanf(KnChInpFile, "AreaKnChFile: %255s\n", AreaKnChFile);
        fscanf(KnChInpFile, "NbVolumesKnCh: %d\n", &NbVolumesKnCh);
        fscanf(KnChInpFile, "VolumeKnChFile: %255s\n", VolumeKnChFile);
        fscanf(KnChInpFile, "KnChFactor: %lg\n", &KnChFactor);
        if (1) {
          printf("NbPointsKnCh: %d\n", NbPointsKnCh);
          printf("PointKnChFile: %s\n", PointKnChFile);
          printf("NbLinesKnCh: %d\n", NbLinesKnCh);
          printf("LineKnChFile: %s\n", LineKnChFile);
          printf("NbAreasKnCh: %d\n", NbAreasKnCh);
          printf("AreaKnChFile: %s\n", AreaKnChFile);
          printf("NbVolumesKnCh: %d\n", NbVolumesKnCh);
          printf("VolumeKnChFile: %s\n", VolumeKnChFile);
          printf("KnChFactor: %lg\n", KnChFactor);
        }

        if (NbPointsKnCh) {  // Points
          if (debugFn)
            printf("No. of points with known charges: %d\n", NbPointsKnCh);

          FILE *fptrPointKnChFile = fopen(PointKnChFile, "r");
          if (fptrPointKnChFile == NULL) {
            neBEMMessage("PointKnCh file absent ... returning\n");
            return -10;
          }

          PointKnChArr =
              (PointKnCh *)malloc((NbPointsKnCh + 1) * sizeof(PointKnCh));
          if (PointKnChArr == NULL) {
            neBEMMessage("Memory allocation failed ... returning\n");
            fclose(fptrPointKnChFile);
            return -10;
          }
          for (int point = 0; point <= NbPointsKnCh; ++point) {
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            PointKnChArr[point].Nb = 0;
            PointKnChArr[point].P.X = 0.0;
            PointKnChArr[point].P.Y = 0.0;
            PointKnChArr[point].P.Z = 0.0;
            PointKnChArr[point].Assigned = 0.0;
          }

          char header[256];
          fgets(header, 256, fptrPointKnChFile);  // header

          for (int point = 1; point <= NbPointsKnCh; ++point) {
            fscanf(fptrPointKnChFile, "%d %lg %lg %lg %lg\n",
                   &PointKnChArr[point].Nb, &PointKnChArr[point].P.X,
                   &PointKnChArr[point].P.Y, &PointKnChArr[point].P.Z,
                   &PointKnChArr[point].Assigned);

            PointKnChArr[point].Assigned *= KnChFactor;
          }

          fclose(fptrPointKnChFile);
          if (debugFn) printf("done for all points\n");
        }  // if NbPointsKnCh: Inputs and calculations for points ends

        if (NbLinesKnCh) {  // Lines
          if (debugFn)
            printf("No. of lines with known charges: %d\n", NbLinesKnCh);

          FILE *fptrLineKnChFile = fopen(LineKnChFile, "r");
          if (fptrLineKnChFile == NULL) {
            neBEMMessage("LineKnCh file absent ... returning\n");
            return -10;
          }

          LineKnChArr =
              (LineKnCh *)malloc((NbLinesKnCh + 1) * sizeof(LineKnCh));
          if (LineKnChArr == NULL) {
            neBEMMessage("Memory allocation failed ... returning\n");
            fclose(fptrLineKnChFile);
            return -10;
          }
          for (int line = 0; line <= NbLinesKnCh; ++line) {  
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            LineKnChArr[line].Nb = 0;
            LineKnChArr[line].Start.X = 0.0;
            LineKnChArr[line].Start.Y = 0.0;
            LineKnChArr[line].Start.Z = 0.0;
            LineKnChArr[line].Stop.X = 0.0;
            LineKnChArr[line].Stop.Y = 0.0;
            LineKnChArr[line].Stop.Z = 0.0;
            LineKnChArr[line].Radius = 0.0;
            LineKnChArr[line].Assigned = 0.0;
          }

          char header[256];
          fgets(header, 256, fptrLineKnChFile);  // header

          for (int line = 1; line <= NbLinesKnCh; ++line) {
            fscanf(fptrLineKnChFile, "%d %lg %lg %lg %lg %lg %lg %lg %lg\n",
                   &LineKnChArr[line].Nb, &LineKnChArr[line].Start.X,
                   &LineKnChArr[line].Start.Y, &LineKnChArr[line].Start.Z,
                   &LineKnChArr[line].Stop.X, &LineKnChArr[line].Stop.Y,
                   &LineKnChArr[line].Stop.Z, &LineKnChArr[line].Radius,
                   &LineKnChArr[line].Assigned);
            LineKnChArr[line].Assigned *= KnChFactor;
          }

          fclose(fptrLineKnChFile);
          if (debugFn) printf("done for all lines\n");
        }  // if NbLinesKnCh: Inputs and calculations for lines ends

        if (NbAreasKnCh) {  // Areas
          if (debugFn)
            printf("No. of areas with known charges: %d\n", NbAreasKnCh);

          FILE *fptrAreaKnChFile = fopen(AreaKnChFile, "r");
          if (fptrAreaKnChFile == NULL) {
            neBEMMessage("AreaKnCh file absent ... returning\n");
            return -10;
          }

          AreaKnChArr =
              (AreaKnCh *)malloc((NbAreasKnCh + 1) * sizeof(AreaKnCh));
          if (AreaKnChArr == NULL) {
            neBEMMessage("Memory allocation failed ... returning\n");
            fclose(fptrAreaKnChFile);
            return -10;
          }
          for (int area = 0; area <= NbAreasKnCh; ++area) {  
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            AreaKnChArr[area].Nb = 0;
            AreaKnChArr[area].NbVertices = 0;
            for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices;
                 ++vert) {
              (AreaKnChArr + area - 1)->Vertex[vert].X = 0.0;
              (AreaKnChArr + area - 1)->Vertex[vert].Y = 0.0;
              (AreaKnChArr + area - 1)->Vertex[vert].Z = 0.0;
            }
            AreaKnChArr[area].Assigned = 0.0;
          }

          char header[256];
          fgets(header, 256, fptrAreaKnChFile);  // header

          for (int area = 1; area <= NbAreasKnCh; ++area) {
            fscanf(fptrAreaKnChFile, "%d %d %le\n",
                   &(AreaKnChArr + area - 1)->Nb,
                   &(AreaKnChArr + area - 1)->NbVertices,
                   &(AreaKnChArr + area - 1)->Assigned);
            for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices;
                 ++vert) {
              fscanf(fptrAreaKnChFile, "%le %le %le\n",
                     &(AreaKnChArr + area - 1)->Vertex[vert].X,
                     &(AreaKnChArr + area - 1)->Vertex[vert].Y,
                     &(AreaKnChArr + area - 1)->Vertex[vert].Z);
            }
            AreaKnChArr[area].Assigned *= KnChFactor;
          }

          fclose(fptrAreaKnChFile);
          if (debugFn) printf("done for all areas\n");
        }  // if AreasKnCh: Inputs and calculations for areas ends

        if (NbVolumesKnCh) {  // Volumes
          if (debugFn)
            printf("No. of volumes with known charges: %d\n", NbVolumesKnCh);

          FILE *fptrVolumeKnChFile = fopen(VolumeKnChFile, "r");
          if (fptrVolumeKnChFile == NULL) {
            neBEMMessage("VolumeKnCh file absent ... returning\n");
            return -10;
          }

          VolumeKnChArr =
              (VolumeKnCh *)malloc((NbVolumesKnCh + 1) * sizeof(VolumeKnCh));
          if (VolumeKnChArr == NULL) {
            neBEMMessage("Memory allocation failed ... returning\n");
            fclose(fptrVolumeKnChFile);
            return -10;
          }
          for (int volume = 0; volume <= NbVolumesKnCh; ++volume) { 
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            VolumeKnChArr[volume].Nb = 0;
            VolumeKnChArr[volume].NbVertices = 0;
            for (int vert = 1; vert <= (VolumeKnChArr + volume - 1)->NbVertices;
                 ++vert) {
              (VolumeKnChArr + volume - 1)->Vertex[vert].X = 0.0;
              (VolumeKnChArr + volume - 1)->Vertex[vert].Y = 0.0;
              (VolumeKnChArr + volume - 1)->Vertex[vert].Z = 0.0;
            }
            VolumeKnChArr[volume].Assigned = 0.0;
          }

          char header[256];
          fgets(header, 256, fptrVolumeKnChFile);  // header

          for (int volume = 1; volume <= NbVolumesKnCh; ++volume) {
            fscanf(fptrVolumeKnChFile, "%d %d %le\n",
                   &(VolumeKnChArr + volume - 1)->Nb,
                   &(VolumeKnChArr + volume - 1)->NbVertices,
                   &(VolumeKnChArr + volume - 1)->Assigned);
            for (int vert = 1; vert <= (VolumeKnChArr + volume - 1)->NbVertices;
                 ++vert) {
              fscanf(fptrVolumeKnChFile, "%le %le %le\n",
                     &(VolumeKnChArr + volume - 1)->Vertex[vert].X,
                     &(VolumeKnChArr + volume - 1)->Vertex[vert].Y,
                     &(VolumeKnChArr + volume - 1)->Vertex[vert].Z);
            }
            VolumeKnChArr[volume].Assigned *= KnChFactor;
          }

          fclose(fptrVolumeKnChFile);
          if (debugFn) printf("done for all volumes\n");
        }  // if NbVolumesKnCh: Inputs and calculations for volumes ends

      }  // if OptKnCh

      fclose(KnChInpFile);
    }  // else KnChInpFile
  }    // parameters related to known charge calculations

  return (0);
}  // neBEMKnownCharges ends

// It is quite possible that the elements had a charge assgined to them
// even before they are being considered for getting charged up.
// Moreover, following the sequence in which the algorithm has been developed,
// the same element accumulates the available electrons and then the ions.
// The resultant of all three (prior charge, electrons and ions) turns out
// to be the assigned charge on the elements, after the execution of this
// function.
int neBEMChargingUp(int /*InfluenceMatrixFlag*/) {
  int debugFn = 0;

  // status of elements before being charged up
  if (debugFn) {
    for (int ele = 1; ele <= NbElements; ++ele) {
      printf("ele, Assigned charge: %d, %lg\n", ele,
             (EleArr + ele - 1)->Assigned);
    }
  }

  // Set up parameters related to charging-up calculations
  // The plan is to distribute the electrons and ions ending in dielectric
  // volumes to the elements on the volumes
  // Only the electrons, to begin with
  {
    FILE *ChargingUpInpFile = fopen("neBEMChargingUp.inp", "r");
    if (ChargingUpInpFile == NULL) {
      printf(
          "neBEMChargingUp.inp absent ... assuming no charging up effect "
          "...\n");
      // assign NbChUpEEle and NbChUpIEle and Prims to be zeros?
    } else {
      fscanf(ChargingUpInpFile, "OptChargingUp: %d\n", &OptChargingUp);
      if (!OptChargingUp)
        printf("OptChargingUp = 0 ... assuming no charging up effect ...\n");
      if (1) printf("OptChargingUp: %d\n", OptChargingUp);

      if (OptChargingUp) {
        char ChargingUpEFile[256];
        char ChargingUpIFile[256];
        double ChUpFactor;
        int *NbChUpEonEle, *NbChUpIonEle;

        fscanf(ChargingUpInpFile, "ChargingUpEFile: %255s\n", ChargingUpEFile);
        fscanf(ChargingUpInpFile, "ChargingUpIFile: %255s\n", ChargingUpIFile);
        fscanf(ChargingUpInpFile, "ChUpFactor: %lg\n", &ChUpFactor);
        if (1) {
          printf("ChargingUpEFile: %s\n", ChargingUpEFile);
          printf("ChargingUpIFile: %s\n", ChargingUpIFile);
          printf("ChUpFactor: %lg\n", ChUpFactor);
        }

        {  // Calculation for electrons
          FILE *fptrChargingUpEFile = fopen(ChargingUpEFile, "r");
          if (fptrChargingUpEFile == NULL) {
            neBEMMessage("ChargingUpE file absent ... returning\n");
            return -10;
          }
          int NbOfE = neBEMGetNbOfLines(ChargingUpEFile);
          if (NbOfE <= 1) {
            neBEMMessage("Too few lines in ChargingUpE ... returning\n");
            fclose(fptrChargingUpEFile);
            return -11;
          }
          NbChUpEonEle = (int *)malloc((NbElements + 1) * sizeof(int));
          if (NbChUpEonEle == NULL) {
            neBEMMessage("Memory allocation failed ... returning\n");
            fclose(fptrChargingUpEFile);
            return -11;
          }            
          for (int ele = 0; ele <= NbElements; ++ele) {
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            NbChUpEonEle[ele] = 0;
          }

          // read the header line
          char header[256];
          fgets(header, 256, fptrChargingUpEFile);

          --NbOfE;  // one line was for the header
          if (debugFn) printf("No. of electrons: %d\n", NbOfE);
          const char *tmpEFile = "/tmp/ElectronTempFile.out";
          FILE *ftmpEF = fopen(tmpEFile, "w");
          if (ftmpEF == NULL) {
            printf("cannot open temporary output file ... returning ...\n");
            free(NbChUpEonEle);
            return -100;
          }
          FILE *fPtEChUpMap = fopen("PtEChUpMap.out", "w");
          if (fPtEChUpMap == NULL) {
            printf("cannot open PtEChUpMap.out file for writing ...\n");
            free(NbChUpEonEle);
            fclose(ftmpEF);
            return 110;
          }

          char label;
          int vol, enb;  // label, volume and electron number
          double xlbend, ylbend, zlbend; // lbend == Last But END
          double xend, yend, zend;
          Point3D
              ptintsct;  // each electron likely to have an intersection point
          for (int electron = 1; electron <= NbOfE; ++electron) {
            fscanf(fptrChargingUpEFile, "%c %d %d %lg %lg %lg %lg %lg %lg\n",
                   &label, &vol, &enb, &xlbend, &xend, &ylbend, &yend, &zlbend,
                   &zend);
            xlbend *= 0.01;
            xend *= 0.01;
            ylbend *= 0.01;
            yend *= 0.01;
            zlbend *= 0.01;
            zend *= 0.01;
            ptintsct.X = 0.0;
            ptintsct.Y = 0.0;
            ptintsct.Z = 0.0;  // initialize

            // find the parametric equation of this last segment
            // if xend < xlbend, swap the directions
            // This has not been mentioned as mandatory in Vince's book
            // "Geometry for Computer Graphics", but implied in the book "A
            // Programmer's Geometry"
            if (xend < xlbend)  // not implemented now
            {
            }
            double lseg = (xend - xlbend) * (xend - xlbend) +
                          (yend - ylbend) * (yend - ylbend) +
                          (zend - zlbend) * (zend - zlbend);
            lseg = sqrt(lseg);
            double xgrd =
                (xend - xlbend) / lseg;  // normalized direction vector
            double ygrd = (yend - ylbend) / lseg;
            double zgrd = (zend - zlbend) / lseg;
            if (debugFn) {
              printf("\nelectron: %d\n", electron);
              printf("xlbend: %lg, ylbend: %lg, zlbend: %lg\n", xlbend, ylbend,
                     zlbend);
              printf("xend: %lg, yend: %lg, zend: %lg\n", xend, yend, zend);
              printf("xgrd: %lg, ygrd: %lg, zgrd: %lg\n", xgrd, ygrd, zgrd);
              fprintf(ftmpEF, "#e: %d, label: %c, vol: %d\n", electron, label,
                      vol);
              fprintf(ftmpEF, "%lg %lg %lg\n", xlbend, ylbend, zlbend);
              fprintf(ftmpEF, "%lg %lg %lg\n", xend, yend, zend);
              fprintf(ftmpEF, "#xgrd: %lg, ygrd: %lg, zgrd: %lg\n", xgrd, ygrd,
                      zgrd);
              fprintf(ftmpEF, "\n");
            }

            // determine which element gets this electron
            // At present, the logic is as follow:
            // Using the information on the last segment, find out which
            // primitive is pierced by it and at which point From the
            // intersection point, find out the element number on the primitive
            // in a local sense Using the information (start and end elements of
            // a given primitive) identify the element in a global sense. This
            // approach should be lot more efficient than checking intersection
            // element by element.
            // The intersection point is computed following algorithm
            // implemented in a matlab code (plane_imp_line_par_int_3d.m) Also
            // check which primitive in the list is the closet to the end point

            double SumOfAngles;
            int PrimIntsctd = -1,
                EleIntsctd = -1;   // intersected primitive & element
            int nearestprim = -1;  // absurd value
            double dist = 1.0e6, mindist = 1.0e6;  // absurdly high numbers
            // check all primitives
            for (int prim = 1; prim <= NbPrimitives; ++prim) {
              if (InterfaceType[prim] != 4)
                continue;  // primitive not a dielectric

              int intersect = 0, extrasect = 1;  // worst of conditions
              int InPrim = 0, InEle = 0;
              if (debugFn)
                printf("prim: %d, mindist: %lg, nearestprim: %d\n", prim,
                       mindist, nearestprim);

              // Use two nodes at a time to get two independent vectors on
              // primitive Get cross-product of these two vector - normal to the
              // plane Note that the normal is already associated with the
              // primitive of type 3 and 4; 2 is wire and does not have any
              // associated normal
              if ((PrimType[prim] == 3) || (PrimType[prim] == 4)) {
                if (debugFn) {
                  printf("prim: %d\n", prim);
                  printf("vertex0: %lg, %lg, %lg\n", XVertex[prim][0],
                         YVertex[prim][0], ZVertex[prim][0]);
                  printf("vertex1: %lg, %lg, %lg\n", XVertex[prim][1],
                         YVertex[prim][1], ZVertex[prim][1]);
                  printf("vertex2: %lg, %lg, %lg\n", XVertex[prim][2],
                         YVertex[prim][2], ZVertex[prim][2]);
                  if (PrimType[prim] == 4) {
                    printf("vertex3: %lg, %lg, %lg\n", XVertex[prim][3],
                           YVertex[prim][3], ZVertex[prim][3]);
                  }
                  fprintf(ftmpEF, "#prim: %d\n", prim);
                  fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][0],
                          YVertex[prim][0], ZVertex[prim][0]);
                  fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][1],
                          YVertex[prim][1], ZVertex[prim][1]);
                  fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][2],
                          YVertex[prim][2], ZVertex[prim][2]);
                  if (PrimType[prim] == 4) {
                    fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][3],
                            YVertex[prim][3], ZVertex[prim][3]);
                  }
                  fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][0],
                          YVertex[prim][0], ZVertex[prim][0]);
                  fprintf(ftmpEF, "\n");
                  fflush(stdout);
                }  // debugFn

                // use a, b, c (normal is ai + bj + ck) at one of the nodes to
                // get d, ax + by + cz + d = 0 is the equation of the plane
                double a = XNorm[prim];
                double b = YNorm[prim];
                double c = ZNorm[prim];
                double d = -a * XVertex[prim][0] - b * YVertex[prim][0] -
                           c * ZVertex[prim][0];

                // distance of the end point to this primitve is
                dist = (xend * a + yend * b + zend * c + d) /
                       sqrt(a * a + b * b + c * c);
                dist = fabs(dist);  // if only magnitude is required
                if (prim == 1) {
                  mindist = dist;
                  nearestprim = prim;
                }
                if ((prim == 1) && debugFn)
                  printf(
                      "after prim == 1 check mindist: %lg, nearestprim: %d\n",
                      mindist, nearestprim);

                // Point of intersection
                // Algo as on p62 (pdf 81) of Vince - Geometry for Computer
                // Graphics 1.5.13 Intersection of a line and a plane Algorithm
                // as implemented in plne_imp_line_par_int_3d.m a (nx), b (ny),
                // c (nz), d are a, b, c, d vx, vy, vz are xgrd, ygrd and zgrd
                // tx, ty, tz are xlbend, ylbend, zlbend
                // In the present case, n and v are unit vectors
                double norm1 = sqrt(a * a + b * b + c * c);
                double norm2 = sqrt(xgrd * xgrd + ygrd * ygrd + zgrd * zgrd);
                double denom =
                    a * xgrd + b * ygrd + c * zgrd;  // (vec)n.(vec)v; if 0, ||
                double tol =
                    1.0e-16;  // CHECK: -8 in original code; sizes small here
                intersect = extrasect = 0;

                if (debugFn) {
                  printf("a, b, c, d, dist: %lg, %lg, %lg, %lg, %lg\n", a, b, c,
                         d, dist);
                  printf("vector n: ai + bj + ck\n");
                  printf("vector v: xgrd, ygrd, zgrd: %lg, %lg, %lg\n", xgrd,
                         ygrd, zgrd);
                  printf("norm1, norm2, (vec n . vec v) denom: %lg, %lg, %lg\n",
                         norm1, norm2, denom);
                  printf("if vec n . vec v == 0, line and plane parallel\n");
                  fflush(stdout);
                }

                if (fabs(denom) < tol * norm1 * norm2) {
                  // line parallel to the plane
                  if (fabs(a * xlbend + b * ylbend + c * zlbend + d) <=
                      1.0e-16) {  // CHECK: was == 0.0 in original code
                    intersect = 1;
                    extrasect = 0;  // line ends on the plane
                    ptintsct.X = xlbend;
                    ptintsct.Y = ylbend;
                    ptintsct.Z = zlbend;
                  } else {
                    intersect = 0;
                    extrasect = 1;     // both wrong
                    ptintsct.X = 0.0;  // Wrong to assign 0 values
                    ptintsct.Y =
                        0.0;  // However, they are never going to be used
                    ptintsct.Z = 0.0;  // since intersect is 0
                  }
                  if (debugFn) {
                    printf("line and plane parallel ...\n");
                    printf("intersect: %d, extrasect: %d\n", intersect,
                           extrasect);
                    printf("intersection point: %lg, %lg, %lg\n", ptintsct.X,
                           ptintsct.Y, ptintsct.Z);
                  }       // if line and plane are parallel
                } else {  // if they are not parallel, they must intersect
                  intersect = 1;
                  double t =
                      -(a * xlbend + b * ylbend + c * zlbend + d) / denom;

                  // check whether t is less than the length of the segment
                  // and in the correct direction
                  // If not, then an extrapolated intersection is not of
                  // interest
                  if ((t < 0.0) ||
                      (fabs(t) > fabs(lseg)))  // wrong dirn or beyond end
                  {
                    extrasect = 1;
                    ptintsct.X = xlbend + t * xgrd;
                    ptintsct.Y = ylbend + t * ygrd;
                    ptintsct.Z = zlbend + t * zgrd;
                  } else {
                    extrasect = 0;
                    ptintsct.X = xlbend + t * xgrd;
                    ptintsct.Y = ylbend + t * ygrd;
                    ptintsct.Z = zlbend + t * zgrd;
                  }
                  if (debugFn) {
                    printf("line and plane NOT parallel ...\n");
                    printf("intersect: %d, extrasect: %d\n", intersect,
                           extrasect);
                    printf("intersection point: %lg, %lg, %lg\n", ptintsct.X,
                           ptintsct.Y, ptintsct.Z);
                    printf("t, lseg: %lg, %lg\n", t, lseg);
                    printf(
                        "for an interesting intersection, lseg > t > 0.0 "
                        "...\n\n");
                    fflush(stdout);
                  }   // must intersect
                }     // if not parallel
              }       // if PrimType is 3 or 4
              else {  // this is a wire primitive - assume no charging up issues
                dist = -1.0;  // an absurd negative distance
                intersect = 0;
                extrasect = 0;
                continue;
              }  // else PrimType 3 or 4

              if (dist < mindist) {
                mindist = dist;
                nearestprim = prim;
              }
              if (debugFn)
                printf("nearestprim: %d, mindist: %lg\n\n", nearestprim,
                       mindist);

              // implicit assumption: the first primitive that gets pierced by
              // the ray is the one that we want. There can be other primitives
              // that are pierced by the same ray. So, this logic should be
              // refined further
              if ((intersect == 1) && (extrasect == 0)) {
                // check whether the intersection point is within primitive
                // polygon
                int nvert = PrimType[prim];
                Point3D polynode[4];
                polynode[0].X = XVertex[prim][0];
                polynode[0].Y = YVertex[prim][0];
                polynode[0].Z = ZVertex[prim][0];
                polynode[1].X = XVertex[prim][1];
                polynode[1].Y = YVertex[prim][1];
                polynode[1].Z = ZVertex[prim][1];
                polynode[2].X = XVertex[prim][2];
                polynode[2].Y = YVertex[prim][2];
                polynode[2].Z = ZVertex[prim][2];
                if (PrimType[prim] == 4) {
                  polynode[3].X = XVertex[prim][3];
                  polynode[3].Y = YVertex[prim][3];
                  polynode[3].Z = ZVertex[prim][3];
                }
                // printf("neBEMChkInPoly for primitive %d\n", prim);
                SumOfAngles = neBEMChkInPoly(nvert, polynode, ptintsct);
                if (fabs(fabs(SumOfAngles) - neBEMtwopi) <= 1.0e-8) {
                  InPrim = 1;
                  PrimIntsctd = prim;
                }
                if (debugFn) {
                  // print polynode and InPrim
                  printf("Prim: %d\n", prim);
                  printf("ptintsct: %lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                         ptintsct.Z);
                  printf("nvert: %d\n", nvert);
                  printf("polynode0: %lg, %lg, %lg\n", polynode[0].X,
                         polynode[0].Y, polynode[0].Z);
                  printf("polynode1: %lg, %lg, %lg\n", polynode[1].X,
                         polynode[1].Y, polynode[1].Z);
                  printf("polynode2: %lg, %lg, %lg\n", polynode[2].X,
                         polynode[2].Y, polynode[2].Z);
                  if (nvert == 4) {
                    printf("polynode3: %lg, %lg, %lg\n", polynode[3].X,
                           polynode[3].Y, polynode[3].Z);
                  }
                  printf("SumOfAngles: %lg, InPrim: %d\n", SumOfAngles, InPrim);
                  fflush(stdout);
                }

                if (!InPrim && (prim != NbPrimitives)) {
                  continue;  // check next primitive
                }

                // Once identified, check in which element belonging to this
                // primitive contains the point of intersection
                if (InPrim) {
                  InEle = 0;
                  for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim];
                       ++ele) {
                    nvert = 0;
                    if ((EleArr + ele - 1)->G.Type == 3) nvert = 3;
                    if ((EleArr + ele - 1)->G.Type == 4) nvert = 4;
                    if (!nvert) {
                      neBEMMessage(
                          "no vertex in element! ... neBEMKnownCharges ...\n");
                      fclose(fPtEChUpMap);
                      return -20;
                    }

                    polynode[0].X = (EleArr + ele - 1)->G.Vertex[0].X;
                    polynode[0].Y = (EleArr + ele - 1)->G.Vertex[0].Y;
                    polynode[0].Z = (EleArr + ele - 1)->G.Vertex[0].Z;
                    polynode[1].X = (EleArr + ele - 1)->G.Vertex[1].X;
                    polynode[1].Y = (EleArr + ele - 1)->G.Vertex[1].Y;
                    polynode[1].Z = (EleArr + ele - 1)->G.Vertex[1].Z;
                    polynode[2].X = (EleArr + ele - 1)->G.Vertex[2].X;
                    polynode[2].Y = (EleArr + ele - 1)->G.Vertex[2].Y;
                    polynode[2].Z = (EleArr + ele - 1)->G.Vertex[2].Z;
                    if (nvert == 4) {
                      polynode[3].X = (EleArr + ele - 1)->G.Vertex[3].X;
                      polynode[3].Y = (EleArr + ele - 1)->G.Vertex[3].Y;
                      polynode[3].Z = (EleArr + ele - 1)->G.Vertex[3].Z;
                    }

                    // printf("neBEMChkInPoly for element %d\n", ele);
                    SumOfAngles = neBEMChkInPoly(nvert, polynode, ptintsct);
                    if (fabs(fabs(SumOfAngles) - neBEMtwopi) <= 1.0e-8)
                      InEle = 1;
                    if (debugFn) {
                      // print polynode and InEle
                      printf("Ele: %d\n", ele);
                      printf("ptintsct: %lg, %lg, %lg\n", ptintsct.X,
                             ptintsct.Y, ptintsct.Z);
                      printf("nvert: %d\n", nvert);
                      printf("polynode0: %lg, %lg, %lg\n", polynode[0].X,
                             polynode[0].Y, polynode[0].Z);
                      printf("polynode1: %lg, %lg, %lg\n", polynode[1].X,
                             polynode[1].Y, polynode[1].Z);
                      printf("polynode2: %lg, %lg, %lg\n", polynode[2].X,
                             polynode[2].Y, polynode[2].Z);
                      if (nvert == 4) {
                        printf("polynode3: %lg, %lg, %lg\n", polynode[3].X,
                               polynode[3].Y, polynode[3].Z);
                      }
                      printf("SumOfAngles: %lg, InEle: %d\n", SumOfAngles,
                             InEle);
                      fflush(stdout);
                    }
                    if (InEle) {
                      ptintsct.X = (EleArr + ele - 1)->G.Origin.X;
                      ptintsct.Y = (EleArr + ele - 1)->G.Origin.Y;
                      ptintsct.Z = (EleArr + ele - 1)->G.Origin.Z;
                      // Associate this electron to the identified element
                      EleIntsctd = ele;
                      NbChUpEonEle[ele]++;
                      fprintf(fPtEChUpMap, "%d %lg %lg %lg %d %d %d %d\n",
                              electron, ptintsct.X, ptintsct.Y, ptintsct.Z,
                              prim, InPrim, ele, InEle);

                      if (debugFn) {
                        printf("# electron: %d\n", electron);
                        printf("%lg %lg %lg\n", xlbend, ylbend, zlbend);
                        printf("%lg %lg %lg\n", xend, yend, zend);
                        printf("%lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                               ptintsct.Z);
                        printf("# Associated primitive: %d\n", prim);
                        printf(
                            "# Associated element and origin: %d, %lg, %lg, "
                            "%lg\n",
                            ele, (EleArr + ele - 1)->G.Origin.X,
                            (EleArr + ele - 1)->G.Origin.Y,
                            (EleArr + ele - 1)->G.Origin.Z);
                        printf("#NbChUpEonEle on element: %d\n",
                               NbChUpEonEle[ele]);
                        fprintf(ftmpEF, "#Element: %d\n", ele);
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X,
                                polynode[0].Y, polynode[0].Z);
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[1].X,
                                polynode[1].Y, polynode[1].Z);
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[2].X,
                                polynode[2].Y, polynode[2].Z);
                        if (nvert == 4) {
                          fprintf(ftmpEF, "%lg %lg %lg\n", polynode[3].X,
                                  polynode[3].Y, polynode[3].Z);
                        }
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X,
                                polynode[0].Y, polynode[0].Z);
                        fprintf(ftmpEF, "\n");
                        fflush(stdout);
                      }
                      break;  // desired element has been found!
                    }
                  }  // for all elements on this primitive

                  if (InEle)
                    break;
                  else {
                    neBEMMessage(
                        "Element not identified ... neBEMKnownCharges\n");
                    return -2;
                  }
                }  // if InPrim
              }    // if intersection and no extrasection

              if ((InPrim) && (intersect) && (!extrasect) && (InEle)) {
                // all satisfied
                // do not check any further primitive for this electron
                break;  
              }

              // If, after checking all the primitives, no interstion is found
              // valid
              if (prim ==
                  (NbPrimitives))  // end of the list and no intersection
              {
                int nvert;
                Point3D polynode[4];
                int nearestele = ElementBgn[nearestprim];
                double distele = 1.0e6;
                double mindistele = 1.0e6;  // absurdly high value

                if (debugFn) {
                  printf("prim == (NbPrimitives) ... checking nearest ...\n");
                  printf("nearestprim: %d, mindist: %lg\n", nearestprim,
                         mindist);
                }

                if (mindist <= 10.0e-6) {
                  PrimIntsctd = nearestprim;
                  InPrim = 1;
                } else {
                  InPrim = 0;
                  InEle = 0;
                  break;
                }
                // check all elements
                for (int ele = ElementBgn[nearestprim];  
                     ele <= ElementEnd[nearestprim]; ++ele) {
                  nvert = 0;
                  if ((EleArr + ele - 1)->G.Type == 3) nvert = 3;
                  if ((EleArr + ele - 1)->G.Type == 4) nvert = 4;
                  if (!nvert) {
                    neBEMMessage(
                        "no vertex element! ... neBEMKnownCharges ...\n");
                    return -20;
                  }

                  /*
                  polynode[0].X = (EleArr+ele-1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr+ele-1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr+ele-1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr+ele-1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr+ele-1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr+ele-1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr+ele-1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr+ele-1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr+ele-1)->G.Vertex[2].Z;
                  if (nvert == 4) {
                    polynode[3].X = (EleArr+ele-1)->G.Vertex[3].X;
                    polynode[3].Y = (EleArr+ele-1)->G.Vertex[3].Y;
                    polynode[3].Z = (EleArr+ele-1)->G.Vertex[3].Z;
                  }
                  Vector3D v01, v12, elenorm, unitelenorm;
                  v01.X = polynode[1].X - polynode[0].X;
                  v01.Y = polynode[1].Y - polynode[0].Y;
                  v01.Z = polynode[1].Z - polynode[0].Z;
                  v12.X = polynode[2].X - polynode[1].X;
                  v12.Y = polynode[2].Y - polynode[1].Y;
                  v12.Z = polynode[2].Z - polynode[1].Z;
                  elenorm = Vector3DCrossProduct(&v01, &v12);
                  unitelenorm = UnitVector3D(&elenorm);

                  if ((nvert == 3) || (nvert == 4)) {
                    if (debugFn) { 
                      printf("nearestprim: %d, element: %d\n",
                             nearestprim, ele); 
                      printf("vertex0: %lg, %lg, %lg\n", polynode[0].X,
                             polynode[0].Y, polynode[0].Z); 
                      printf("vertex1: %lg, %lg, %lg\n", polynode[1].X, 
                             polynode[1].Y, polynode[1].Z); 
                      printf("vertex2: %lg, %lg, %lg\n", polynode[2].X,
                             polynode[2].Y, polynode[2].Z); 
                      if (PrimType[prim] == 4) { 
                        printf("vertex3: %lg, %lg, %lg\n", polynode[3].X, 
                               polynode[3].Y, polynode[3].Z);
                      }
                      fprintf(ftmpEF, "#nearestprim: %d, element: %d\n", 
                              nearestprim, ele); 
                      fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X,
                              polynode[0].Y, polynode[0].Z); 
                      fprintf(ftmpEF, "%lg %lg %lg\n", polynode[1].X,
                              polynode[1].Y, polynode[1].Z); 
                      fprintf(ftmpEF, "%lg %lg %lg\n", polynode[2].X,
                              polynode[2].Y, polynode[2].Z); 
                      if (PrimType[prim] == 4) { 
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[3].X, 
                                polynode[3].Y, polynode[3].Z);
                      }
                      fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X, 
                              polynode[0].Y, polynode[0].Z); 
                      fprintf(ftmpEF, "\n"); 
                      fflush(stdout); 
                    } // debugFn
                    // use a, b, c (normal is ai + bj + ck) 
                    // at one of the nodes to get d
                    // ax + by + cz + d = 0 is the equation of the plane
                    double a = unitelenorm.X;
                    double b = unitelenorm.Y;
                    double c = unitelenorm.Z;
                    double d = - a * polynode[0].X - b * polynode[0].Y - c * polynode[0].Z;
                    // distance of the end point to this primitve is
                    distele = (xend * a + yend * b + zend * c + d) /
                              sqrt(a * a + b * b + c * c);
                    distele = fabs(distele); // if only magnitude is required
                    */

                  Vector3D eleOrigin;
                  eleOrigin.X = (EleArr + ele - 1)->G.Origin.X;
                  eleOrigin.Y = (EleArr + ele - 1)->G.Origin.Y;
                  eleOrigin.Z = (EleArr + ele - 1)->G.Origin.Z;
                  distele = (eleOrigin.X - xend) * (eleOrigin.X - xend) +
                            (eleOrigin.Y - yend) * (eleOrigin.Y - yend) +
                            (eleOrigin.Z - zend) * (eleOrigin.Z - zend);
                  distele = pow(distele, 0.5);

                  if (ele == ElementBgn[nearestprim]) {
                    mindistele = distele;
                    nearestele = ele;
                  }
                  if (distele < mindistele) {
                    mindistele = distele;
                    nearestele = ele;
                  }

                  if (debugFn) {
                    // printf("a, b, c, d, dist: %lg, %lg, %lg, %lg, %lg\n",
                    // a, b, c, d,  dist);
                    // printf("vector n: ai + bj + ck\n");
                    // printf("vector v: xgrd, ygrd, zgrd: %lg, %lg, %lg\n",
                    // xgrd, ygrd, zgrd);
                    printf(
                        "distele: %lg, mindistele: %lg,from nearest ele "
                        "origin: %d\n",
                        distele, mindistele, nearestele);
                    fflush(stdout);
                  }

                  // }	// if PrimType is 3 or 4
                }  // for elements in nearestprim

                // if(mindistele <= 10.0e-6)
                // {
                EleIntsctd = nearestele;
                InEle = 1;
                ptintsct.X = (EleArr + EleIntsctd - 1)->G.Origin.X;
                ptintsct.Y = (EleArr + EleIntsctd - 1)->G.Origin.Y;
                ptintsct.Z = (EleArr + EleIntsctd - 1)->G.Origin.Z;
                NbChUpEonEle[EleIntsctd]++;

                fprintf(fPtEChUpMap, "%d %lg %lg %lg %d %d %d %d\n", electron,
                        ptintsct.X, ptintsct.Y, ptintsct.Z, PrimIntsctd, InPrim,
                        EleIntsctd, InEle);
                // }	// if mindistele

                if (debugFn) {
                  printf("# electron: %d\n", electron);
                  printf("%lg %lg %lg\n", xlbend, ylbend, zlbend);
                  printf("%lg %lg %lg\n", xend, yend, zend);
                  printf("%lg, %lg, %lg\n", ptintsct.X, ptintsct.Y, ptintsct.Z);
                  printf("# Associated primitive: %d\n", PrimIntsctd);
                  printf("# Associated element and origin: %d, %lg, %lg, %lg\n",
                         EleIntsctd, (EleArr + EleIntsctd - 1)->G.Origin.X,
                         (EleArr + EleIntsctd - 1)->G.Origin.Y,
                         (EleArr + EleIntsctd - 1)->G.Origin.Z);
                  printf("#NbChUpEonEle on element: %d\n",
                         NbChUpEonEle[EleIntsctd]);
                  fflush(stdout);

                  fprintf(ftmpEF, "#Element: %d\n", EleIntsctd);
                  polynode[0].X = (EleArr + EleIntsctd - 1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr + EleIntsctd - 1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr + EleIntsctd - 1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr + EleIntsctd - 1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr + EleIntsctd - 1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr + EleIntsctd - 1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr + EleIntsctd - 1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr + EleIntsctd - 1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr + EleIntsctd - 1)->G.Vertex[2].Z;
                  if (nvert == 4) {
                    polynode[3].X = (EleArr + EleIntsctd - 1)->G.Vertex[3].X;
                    polynode[3].Y = (EleArr + EleIntsctd - 1)->G.Vertex[3].Y;
                    polynode[3].Z = (EleArr + EleIntsctd - 1)->G.Vertex[3].Z;
                  }
                  fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X, polynode[0].Y,
                          polynode[0].Z);
                  fprintf(ftmpEF, "%lg %lg %lg\n", polynode[1].X, polynode[1].Y,
                          polynode[1].Z);
                  fprintf(ftmpEF, "%lg %lg %lg\n", polynode[2].X, polynode[2].Y,
                          polynode[2].Z);
                  if (nvert == 4) {
                    fprintf(ftmpEF, "%lg %lg %lg\n", polynode[3].X,
                            polynode[3].Y, polynode[3].Z);
                  }
                  fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X, polynode[0].Y,
                          polynode[0].Z);
                  fprintf(ftmpEF, "\n");
                }  // debug
              }    // if prim == NbPrimitives

            }  // for all primitives // just not those on the volume

            if (debugFn)
              printf("writing file for checking electron positions\n");

            if (debugFn)  // check electron positions, volume primitives and
                          // elements
            {
              char elecposdbg[256], enbstr[10];
              sprintf(enbstr, "%d", electron);
              strcpy(elecposdbg, "/tmp/Electron");
              strcat(elecposdbg, enbstr);
              strcat(elecposdbg, ".out");
              FILE *fepd = fopen(elecposdbg, "w");
              if (fepd == NULL) {
                printf(
                    "cannot open writable file to debug electron positions "
                    "...\n");
                printf("returning ...\n");
                return -111;
              }
              // write electron number, end, lbend, volume, primitive, elements,
              // intxn
              fprintf(fepd, "#electron: %d %d\n", enb,
                      electron);  // should print same
              fprintf(fepd, "#last but end position:\n");
              fprintf(fepd, "%lg %lg %lg\n", xlbend, ylbend, zlbend);
              fprintf(fepd, "#end position:\n");
              fprintf(fepd, "%lg %lg %lg\n\n", xend, yend, zend);
              fprintf(fepd, "#intersected primitive number: %d\n", PrimIntsctd);
              if (PrimIntsctd >= 1) {
                fprintf(fepd, "#PrimType: %d\n", PrimType[PrimIntsctd]);
                fprintf(fepd, "#prim vertices:\n");
                fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][0],
                        YVertex[PrimIntsctd][0], ZVertex[PrimIntsctd][0]);
                fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][1],
                        YVertex[PrimIntsctd][1], ZVertex[PrimIntsctd][1]);
                fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][2],
                        YVertex[PrimIntsctd][2], ZVertex[PrimIntsctd][2]);
                if (PrimType[PrimIntsctd] == 4) {
                  fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][3],
                          YVertex[PrimIntsctd][3], ZVertex[PrimIntsctd][3]);
                }
                fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][0],
                        YVertex[PrimIntsctd][0], ZVertex[PrimIntsctd][0]);
                fprintf(fepd, "\n");

                fprintf(fepd, "#ptintsct:\n");
                fprintf(fepd, "%lg %lg %lg\n", ptintsct.X, ptintsct.Y,
                        ptintsct.Z);
                fprintf(fepd, "\n");
              }
              fprintf(fepd, "#intersected element number: %d\n", EleIntsctd);
              if (EleIntsctd >= 1) {
                int gtype = (EleArr + EleIntsctd - 1)->G.Type;
                fprintf(fepd, "#EleType: %d\n", gtype);
                fprintf(fepd, "#element vertices:\n");
                double x0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].X;
                double y0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].Y;
                double z0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].Z;
                double x1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].X;
                double y1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].Y;
                double z1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].Z;
                double x2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].X;
                double y2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].Y;
                double z2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].Z;
                fprintf(fepd, "%lg %lg %lg\n", x0, y0, z0);
                fprintf(fepd, "%lg %lg %lg\n", x1, y1, z1);
                fprintf(fepd, "%lg %lg %lg\n", x2, y2, z2);
                if (gtype == 4) {
                  double x3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].X;
                  double y3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].Y;
                  double z3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].Z;
                  fprintf(fepd, "%lg %lg %lg\n", x3, y3, z3);
                }
                fprintf(fepd, "%lg %lg %lg\n", x0, y0, z0);
                fprintf(fepd, "\n");

                fprintf(fepd, "#ptintsct:\n");
                fprintf(fepd, "%lg %lg %lg\n", ptintsct.X, ptintsct.Y,
                        ptintsct.Z);
                fprintf(fepd, "\n");
              }

              fclose(fepd);
            }  // if 1
            if (debugFn)
              printf("done writing file for checking electron positions\n");
          }  // for all the electrons
          fclose(fPtEChUpMap);
          if (debugFn) printf("done for all electrons\n");

          FILE *fEleEChUpMap = fopen("EleEChUpMap.out", "w");
          if (fEleEChUpMap == NULL) {
            printf("cannot open EleEChUpMap.out file for writing ...\n");
            return 111;
          }
          for (int ele = 1; ele <= NbElements; ++ele) {
            (EleArr + ele - 1)->Assigned +=
                ChUpFactor * Q_E * NbChUpEonEle[ele] / (EleArr + ele - 1)->G.dA;
            fprintf(fEleEChUpMap, "%d %lg %lg %lg %d %lg\n", ele,
                    (EleArr + ele - 1)->G.Origin.X,
                    (EleArr + ele - 1)->G.Origin.Y,
                    (EleArr + ele - 1)->G.Origin.Z, NbChUpEonEle[ele],
                    (EleArr + ele - 1)->Assigned);
          }
          fclose(fEleEChUpMap);

          fclose(ftmpEF);
          free(NbChUpEonEle);
        }  // Calculation for electrons ends

        {  // Calculation for ions
          FILE *fptrChargingUpIFile = fopen(ChargingUpIFile, "r");
          if (fptrChargingUpIFile == NULL) {
            neBEMMessage("ChargingUpI file absent ... returning\n");
            return -10;
          }
          int NbOfI = neBEMGetNbOfLines(ChargingUpIFile);
          if (NbOfI <= 1) {
            neBEMMessage("Too few lines in ChargingUpI ... returning\n");
            fclose(fptrChargingUpIFile);
            return -11;
          } 
          // initialize
          NbChUpIonEle = (int *)malloc((NbElements + 1) * sizeof(int));
          if (NbChUpIonEle == NULL) {
            neBEMMessage("Memory allocation failed ... returning\n");
            fclose(fptrChargingUpIFile);
            return -11;
          }
          for (int ele = 0; ele <= NbElements; ++ele) {  
            // CHECK!!! ele limit starts from 0 but all other from 1 to ...
            NbChUpIonEle[ele] = 0;
          }

          // read the header line
          char header[256];
          fgets(header, 256, fptrChargingUpIFile);

          --NbOfI;  // one line was for the header
          if (debugFn) printf("No. of ions: %d\n", NbOfI);
          const char *tmpIFile = "/tmp/IonTempFile.out";
          FILE *ftmpIF = fopen(tmpIFile, "w");
          if (ftmpIF == NULL) {
            printf("cannot open temporary ion output file ... returning ...\n");
            free(NbChUpEonEle);
            return -100;
          }
          FILE *fPtIChUpMap = fopen("PtIChUpMap.out", "w");
          if (fPtIChUpMap == NULL) {
            printf("cannot open PtIChUpMap.out file for writing ...\n");
            fclose(ftmpIF);
            free(NbChUpEonEle);
            return 110;
          }

          char label;
          int inb, vol;  // label, volume and ion number
          double xlbend, ylbend, zlbend; // lbend == Last But END
          double xend, yend, zend;
          Point3D ptintsct;  // each ion likely to have an intersection point
          for (int ion = 1; ion <= NbOfI; ++ion) {
            fscanf(fptrChargingUpIFile, "%c %d %d %lg %lg %lg %lg %lg %lg\n",
                   &label, &vol, &inb, &xlbend, &xend, &ylbend, &yend, &zlbend,
                   &zend);
            xlbend *= 0.01;
            xend *= 0.01;
            ylbend *= 0.01;
            yend *= 0.01;
            zlbend *= 0.01;
            zend *= 0.01;
            ptintsct.X = 0.0;
            ptintsct.Y = 0.0;
            ptintsct.Z = 0.0;  // initialize

            // find the parametric equation of this last segment
            // if xend < xlbend, swap the directions
            // This has not been mentioned as mandatory in Vince's book
            // "Geometry for Computer Graphics", but implied in the book "A
            // Programmer's Geometry"
            if (xend < xlbend)  // not implemented now
            {
            }
            double lseg = (xend - xlbend) * (xend - xlbend) +
                          (yend - ylbend) * (yend - ylbend) +
                          (zend - zlbend) * (zend - zlbend);
            lseg = sqrt(lseg);
            double xgrd =
                (xend - xlbend) / lseg;  // normalized direction vector
            double ygrd = (yend - ylbend) / lseg;
            double zgrd = (zend - zlbend) / lseg;
            if (debugFn) {
              printf("\nion: %d\n", ion);
              printf("xlbend: %lg, ylbend: %lg, zlbend: %lg\n", xlbend, ylbend,
                     zlbend);
              printf("xend: %lg, yend: %lg, zend: %lg\n", xend, yend, zend);
              printf("xgrd: %lg, ygrd: %lg, zgrd: %lg\n", xgrd, ygrd, zgrd);
              fprintf(ftmpIF, "#e: %d, label: %c, vol: %d\n", ion, label, vol);
              fprintf(ftmpIF, "%lg %lg %lg\n", xlbend, ylbend, zlbend);
              fprintf(ftmpIF, "%lg %lg %lg\n", xend, yend, zend);
              fprintf(ftmpIF, "#xgrd: %lg, ygrd: %lg, zgrd: %lg\n", xgrd, ygrd,
                      zgrd);
              fprintf(ftmpIF, "\n");
            }

            // determine which element gets this electron
            // At present, the logic is as follow:
            // Using the information on the last segment, find out which
            // primitive is pierced by it and at which point From the
            // intersection point, find out the element number on the primitive
            // in a local sense Using the information (start and end elements of
            // a given primitive) identify the element in a global sense. This
            // approach should be lot more efficient than checking intersection
            // element by element.
            // The intersection point is computed following algorithm
            // implemented in a matlab code (plane_imp_line_par_int_3d.m) Also
            // check which primitive in the list is the closet to the end point

            int PrimIntsctd = -1,
                EleIntsctd = -1;   // intersected primitive & element
            int nearestprim = -1;  // absurd value
            double dist = 1.0e6, mindist = 1.0e6;  // absurdly high numbers
            double SumOfAngles;
            // check all primitives
            for (int prim = 1; prim <= NbPrimitives; ++prim) {
              if (InterfaceType[prim] != 4)
                continue;  // primitive not a dielectric

              int intersect = 0, extrasect = 1;  // worst of conditions
              int InPrim = 0, InEle = 0;
              if (debugFn)
                printf("prim: %d, mindist: %lg, nearestprim: %d\n", prim,
                       mindist, nearestprim);

              // get the primitive nodes

              // Use two nodes at a time to get two independent vectors on
              // primitive Get cross-product of these two vector - normal to the
              // plane Note that the normal is already associated with the
              // primitive of type 3 and 4; 2 is wire and does not have any
              // associated normal
              if ((PrimType[prim] == 3) || (PrimType[prim] == 4)) {
                if (debugFn) {
                  printf("prim: %d\n", prim);
                  printf("vertex0: %lg, %lg, %lg\n", XVertex[prim][0],
                         YVertex[prim][0], ZVertex[prim][0]);
                  printf("vertex1: %lg, %lg, %lg\n", XVertex[prim][1],
                         YVertex[prim][1], ZVertex[prim][1]);
                  printf("vertex2: %lg, %lg, %lg\n", XVertex[prim][2],
                         YVertex[prim][2], ZVertex[prim][2]);
                  if (PrimType[prim] == 4) {
                    printf("vertex3: %lg, %lg, %lg\n", XVertex[prim][3],
                           YVertex[prim][3], ZVertex[prim][3]);
                  }
                  fprintf(ftmpIF, "#prim: %d\n", prim);
                  fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][0],
                          YVertex[prim][0], ZVertex[prim][0]);
                  fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][1],
                          YVertex[prim][1], ZVertex[prim][1]);
                  fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][2],
                          YVertex[prim][2], ZVertex[prim][2]);
                  if (PrimType[prim] == 4) {
                    fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][3],
                            YVertex[prim][3], ZVertex[prim][3]);
                  }
                  fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][0],
                          YVertex[prim][0], ZVertex[prim][0]);
                  fprintf(ftmpIF, "\n");
                  fflush(stdout);
                }  // debugFn

                // use a, b, c (normal is ai + bj + ck) at one of the nodes to
                // get d ax + by + cz + d = 0 is the equation of the plane
                double d = -XNorm[prim] * XVertex[prim][0] -
                           YNorm[prim] * YVertex[prim][0] -
                           ZNorm[prim] * ZVertex[prim][0];

                // distance of the end point to this primitve is
                dist =
                    (xend * XNorm[prim] + yend * YNorm[prim] +
                     zend * ZNorm[prim] + d) /
                    sqrt(XNorm[prim] * XNorm[prim] + YNorm[prim] * YNorm[prim] +
                         ZNorm[prim] * ZNorm[prim]);
                dist = fabs(dist);  // if only magnitude is required
                if (prim == 1) {
                  mindist = dist;
                  nearestprim = prim;
                }
                if ((prim == 1) && debugFn)
                  printf(
                      "after prim == 1 check mindist: %lg, nearestprim: %d\n",
                      mindist, nearestprim);

                // Point of intersection
                // Algo as on p62 (pdf 81) of Vince - Geometry for Computer
                // Graphics 1.5.13 Intersection of a line and a plane Algorithm
                // as implemented in plne_imp_line_par_int_3d.m a (nx), b (ny),
                // c (nz), d are a, b, c, d vx, vy, vz are xgrd, ygrd and zgrd
                // tx, ty, tz are xlbend, ylbend, zlbend
                // In the present case, n and v are unit vectors
                double a = XNorm[prim];
                double b = YNorm[prim];
                double c = ZNorm[prim];
                double norm1 = sqrt(a * a + b * b + c * c);
                double norm2 = sqrt(xgrd * xgrd + ygrd * ygrd + zgrd * zgrd);
                double denom =
                    a * xgrd + b * ygrd + c * zgrd;  // (vec)n.(vec)v; if 0, ||
                double tol =
                    1.0e-12;  // CHECK: -8 in original code; sizes small here
                intersect = extrasect = 0;

                if (debugFn) {
                  printf("a, b, c, d, dist: %lg, %lg, %lg, %lg, %lg\n", a, b, c,
                         d, dist);
                  printf("vector n: ai + bj + ck\n");
                  printf("vector v: xgrd, ygrd, zgrd: %lg, %lg, %lg\n", xgrd,
                         ygrd, zgrd);
                  printf("norm1, norm2, (vec n . vec v) denom: %lg, %lg, %lg\n",
                         norm1, norm2, denom);
                  printf("if vec n . vec v == 0, line and plane parallel\n");
                  fflush(stdout);
                }

                if (fabs(denom) < tol * norm1 * norm2) {
                  // line parallel to the plane
                  if (a * xlbend + b * ylbend + c * zlbend + d ==
                      0.0)  // CHECK == for float
                  {
                    intersect = 1;
                    extrasect = 0;
                    ptintsct.X = xlbend;
                    ptintsct.Y = ylbend;
                    ptintsct.Z = zlbend;
                  } else {
                    intersect = 0;
                    extrasect = 0;
                    ptintsct.X = 0.0;  // Wrong to assign 0 values
                    ptintsct.Y =
                        0.0;  // However, they are never going to be used
                    ptintsct.Z = 0.0;  // since intersect is 0
                  }
                  if (debugFn) {
                    printf("line and plane parallel ...\n");
                    printf("intersect: %d, extrasect: %d\n", intersect,
                           extrasect);
                    printf("intersection point: %lg, %lg, %lg\n", ptintsct.X,
                           ptintsct.Y, ptintsct.Z);
                  }       // if line and plane are parallel
                } else {  // if they are not parallel, they must intersect
                  intersect = 1;
                  double t =
                      -(a * xlbend + b * ylbend + c * zlbend + d) / denom;

                  // check whether t is less than the length of the segment
                  // and in the correct direction
                  // If not, then an extrapolated intersection is not of
                  // interest
                  if ((t < 0.0) || (fabs(t) > fabs(lseg))) {
                    extrasect = 1;
                    ptintsct.X = xlbend + t * xgrd;
                    ptintsct.Y = ylbend + t * ygrd;
                    ptintsct.Z = zlbend + t * zgrd;
                  } else {
                    extrasect = 0;
                    ptintsct.X = xlbend + t * xgrd;
                    ptintsct.Y = ylbend + t * ygrd;
                    ptintsct.Z = zlbend + t * zgrd;
                  }
                  if (debugFn) {
                    printf("line and plane NOT parallel ...\n");
                    printf("intersect: %d, extrasect: %d\n", intersect,
                           extrasect);
                    printf("intersection point: %lg, %lg, %lg\n", ptintsct.X,
                           ptintsct.Y, ptintsct.Z);
                    printf("t, lseg: %lg, %lg\n", t, lseg);
                    printf(
                        "for an interesting intersection, lseg > t > 0.0 "
                        "...\n\n");
                    fflush(stdout);
                  }   // must intersect
                }     // if not parallel
              }       // if PrimType is 3 or 4
              else {  // this is a wire primitive - assume no charging up issues
                dist = -1.0;  // an absurd negative distance
                intersect = 0;
                extrasect = 0;
              }  // else PrimType 3 or 4

              if (dist < mindist) {
                mindist = dist;
                nearestprim = prim;
              }

              // implicit assumption: the first primitive that gets pierced by
              // the ray is the one that we want. There can be other primitives
              // that are pierced by the same ray. So, this logic should be
              // refined further
              if ((intersect == 1) && (extrasect == 0)) {
                // check whether the intersection point is within primitive
                // polygon
                int nvert = PrimType[prim];
                Point3D polynode[4];
                polynode[0].X = XVertex[prim][0];
                polynode[0].Y = YVertex[prim][0];
                polynode[0].Z = ZVertex[prim][0];
                polynode[1].X = XVertex[prim][1];
                polynode[1].Y = YVertex[prim][1];
                polynode[1].Z = ZVertex[prim][1];
                polynode[2].X = XVertex[prim][2];
                polynode[2].Y = YVertex[prim][2];
                polynode[2].Z = ZVertex[prim][2];
                if (PrimType[prim] == 4) {
                  polynode[3].X = XVertex[prim][3];
                  polynode[3].Y = YVertex[prim][3];
                  polynode[3].Z = ZVertex[prim][3];
                }
                // printf("neBEMChkInPoly for primitive %d\n", prim);
                SumOfAngles = neBEMChkInPoly(nvert, polynode, ptintsct);
                if (fabs(fabs(SumOfAngles) - neBEMtwopi) <= 1.0e-8) {
                  InPrim = 1;
                  PrimIntsctd = prim;
                }
                if (debugFn) {
                  // print polynode and InPrim
                  printf("Prim: %d\n", prim);
                  printf("ptintsct: %lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                         ptintsct.Z);
                  printf("nvert: %d\n", nvert);
                  printf("polynode0: %lg, %lg, %lg\n", polynode[0].X,
                         polynode[0].Y, polynode[0].Z);
                  printf("polynode1: %lg, %lg, %lg\n", polynode[1].X,
                         polynode[1].Y, polynode[1].Z);
                  printf("polynode2: %lg, %lg, %lg\n", polynode[2].X,
                         polynode[2].Y, polynode[2].Z);
                  if (nvert == 4) {
                    printf("polynode3: %lg, %lg, %lg\n", polynode[3].X,
                           polynode[3].Y, polynode[3].Z);
                  }
                  printf("SumOfAngles: %lg, InPrim: %d\n", SumOfAngles, InPrim);
                  fflush(stdout);
                }
                if (!InPrim) continue;  // check next primitive

                // Once identified, check in which element belonging to this
                // primitive contains the point of intersection
                InEle = 0;
                for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim];
                     ++ele) {
                  nvert = 0;
                  if ((EleArr + ele - 1)->G.Type == 3) nvert = 3;
                  if ((EleArr + ele - 1)->G.Type == 4) nvert = 4;
                  if (!nvert) {
                    neBEMMessage(
                        "no vertex in element! ... neBEMKnownCharges ...\n");
                    if (fPtIChUpMap) fclose(fPtIChUpMap);
                    return -20;
                  }

                  polynode[0].X = (EleArr + ele - 1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr + ele - 1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr + ele - 1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr + ele - 1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr + ele - 1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr + ele - 1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr + ele - 1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr + ele - 1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr + ele - 1)->G.Vertex[2].Z;
                  if (nvert == 4) {
                    polynode[3].X = (EleArr + ele - 1)->G.Vertex[3].X;
                    polynode[3].Y = (EleArr + ele - 1)->G.Vertex[3].Y;
                    polynode[3].Z = (EleArr + ele - 1)->G.Vertex[3].Z;
                  }

                  // printf("neBEMChkInPoly for element %d\n", ele);
                  SumOfAngles = neBEMChkInPoly(nvert, polynode, ptintsct);
                  if (fabs(fabs(SumOfAngles) - neBEMtwopi) <= 1.0e-8) InEle = 1;
                  if (debugFn) {
                    // print polynode and InEle
                    printf("Ele: %d\n", ele);
                    printf("ptintsct: %lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                           ptintsct.Z);
                    printf("nvert: %d\n", nvert);
                    printf("polynode0: %lg, %lg, %lg\n", polynode[0].X,
                           polynode[0].Y, polynode[0].Z);
                    printf("polynode1: %lg, %lg, %lg\n", polynode[1].X,
                           polynode[1].Y, polynode[1].Z);
                    printf("polynode2: %lg, %lg, %lg\n", polynode[2].X,
                           polynode[2].Y, polynode[2].Z);
                    if (nvert == 4) {
                      printf("polynode3: %lg, %lg, %lg\n", polynode[3].X,
                             polynode[3].Y, polynode[3].Z);
                    }
                    printf("SumOfAngles: %lg, InEle: %d\n", SumOfAngles, InEle);
                    fflush(stdout);
                  }
                  if (InEle) {
                    ptintsct.X = (EleArr + ele - 1)->G.Origin.X;
                    ptintsct.Y = (EleArr + ele - 1)->G.Origin.Y;
                    ptintsct.Z = (EleArr + ele - 1)->G.Origin.Z;
                    EleIntsctd = ele;
                    // Associate this electron to the identified element
                    NbChUpIonEle[ele]++;
                    fprintf(fPtIChUpMap, "%d %lg %lg %lg %d %d %d %d\n", ion,
                            ptintsct.X, ptintsct.Y, ptintsct.Z, prim, InPrim,
                            ele, InEle);

                    if (debugFn) {
                      printf("# ion: %d\n", ion);
                      printf("%lg %lg %lg\n", xlbend, ylbend, zlbend);
                      printf("%lg %lg %lg\n", xend, yend, zend);
                      printf("%lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                             ptintsct.Z);
                      printf("# Associated primitive: %d\n", prim);
                      printf(
                          "# Associated element and origin: %d, %lg, %lg, "
                          "%lg\n",
                          ele, (EleArr + ele - 1)->G.Origin.X,
                          (EleArr + ele - 1)->G.Origin.Y,
                          (EleArr + ele - 1)->G.Origin.Z);
                      printf("#NbChUpIonEle on element: %d\n",
                             NbChUpIonEle[ele]);
                      fprintf(ftmpIF, "#Element: %d\n", ele);
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X,
                              polynode[0].Y, polynode[0].Z);
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[1].X,
                              polynode[1].Y, polynode[1].Z);
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[2].X,
                              polynode[2].Y, polynode[2].Z);
                      if (nvert == 4) {
                        fprintf(ftmpIF, "%lg %lg %lg\n", polynode[3].X,
                                polynode[3].Y, polynode[3].Z);
                      }
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X,
                              polynode[0].Y, polynode[0].Z);
                      fprintf(ftmpIF, "\n");
                      fflush(stdout);
                    }
                    break;  // desired element has been found!
                  }         // if InEle
                }           // for all elements on this primitive

                if (InEle) break;
                neBEMMessage(
                    "Element cannot be identified ... neBEMKnownCharges\n");
                return -2;
              }  // if proper intersection and no extrasection

              if ((InPrim) && (intersect) && (!extrasect) && (InEle)) {
                // all satisfied
                // do not check any further primitive for this electron
                break;  
              }

              // If, after checking all the primitives, no interstion is found
              // valid
              if (prim == (NbPrimitives)) {
                // end of the list and no intersection
                int nvert;
                Point3D polynode[4];
                int nearestele = ElementBgn[nearestprim];
                double distele = 1.0e6;
                double mindistele = 1.0e6;  // absurdly high value

                if (debugFn) {
                  printf("prim == (NbPrimitives) ... checking nearest ...\n");
                  printf("nearestprim: %d, mindist: %lg\n", nearestprim,
                         mindist);
                }

                if (mindist <= 10.0e-6) {
                  PrimIntsctd = nearestprim;
                  InPrim = 1;
                } else {
                  InPrim = 0;
                  InEle = 0;
                  break;
                }

                for (int ele = ElementBgn[nearestprim];  // check all elements
                     ele <= ElementEnd[nearestprim]; ++ele) {
                  nvert = 0;
                  if ((EleArr + ele - 1)->G.Type == 3) nvert = 3;
                  if ((EleArr + ele - 1)->G.Type == 4) nvert = 4;
                  if (!nvert) {
                    neBEMMessage(
                        "no vertex element! ... neBEMKnownCharges ...\n");
                    return -20;
                  }

                  /*
                  polynode[0].X = (EleArr+ele-1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr+ele-1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr+ele-1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr+ele-1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr+ele-1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr+ele-1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr+ele-1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr+ele-1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr+ele-1)->G.Vertex[2].Z;
                  if (nvert == 4) {
                    polynode[3].X = (EleArr+ele-1)->G.Vertex[3].X;
                    polynode[3].Y = (EleArr+ele-1)->G.Vertex[3].Y;
                    polynode[3].Z = (EleArr+ele-1)->G.Vertex[3].Z;
                  } 

                  Vector3D v01, v12, elenorm, unitelenorm;
                  v01.X = polynode[1].X - polynode[0].X;
                  v01.Y = polynode[1].Y - polynode[0].Y;
                  v01.Z = polynode[1].Z - polynode[0].Z;
                  v12.X = polynode[2].X - polynode[1].X;
                  v12.Y = polynode[2].Y - polynode[1].Y;
                  v12.Z = polynode[2].Z - polynode[1].Z;
                  elenorm = Vector3DCrossProduct(&v01, &v12);
                  unitelenorm = UnitVector3D(&elenorm);

                  if ((nvert == 3) || (nvert == 4)) {
                    if (debugFn) {
                      printf("nearestprim: %d, element: %d\n",
                             nearestprim, ele); 
                      printf("vertex0: %lg, %lg, %lg\n", polynode[0].X,
                             polynode[0].Y, polynode[0].Z); 
                      printf("vertex1: %lg, %lg, %lg\n", polynode[1].X,
                             polynode[1].Y, polynode[1].Z); 
                      printf("vertex2: %lg, %lg, %lg\n", polynode[2].X, 
                             polynode[2].Y, polynode[2].Z); 
                      if (PrimType[prim] == 4) {
                        printf("vertex3: %lg, %lg, %lg\n", polynode[3].X,
                               polynode[3].Y, polynode[3].Z);
                      }
                      fprintf(ftmpIF, "#nearestprim: %d, element: %d\n", 
                              nearestprim, ele); 
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X,
                              polynode[0].Y, polynode[0].Z); 
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[1].X,
                              polynode[1].Y, polynode[1].Z); 
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[2].X,
                              polynode[2].Y, polynode[2].Z); 
                      if (PrimType[prim] == 4) {
                        fprintf(ftmpIF, "%lg %lg %lg\n", polynode[3].X, 
                                polynode[3].Y, polynode[3].Z);
                      }
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X,
                              polynode[0].Y, polynode[0].Z); 
                      fprintf(ftmpIF, "\n"); 
                      fflush(stdout);
                    } // debugFn

                    // use a, b, c (normal is ai + bj + ck) 
                    // at one of the nodes to get d
                    // ax + by + cz + d = 0 is the equation of the plane
                    double a = unitelenorm.X;
                    double b = unitelenorm.Y;
                    double c = unitelenorm.Z;
                    double d = - a * polynode[0].X - b * polynode[0].Y - c * polynode[0].Z;
                    // distance of the end point to this primitve is
                    distele = (xend * a + yend * b + zend * c + d) /
                              sqrt(a * a + b * b + c * c);
                    distele = fabs(distele); // if only magnitude is required
                  */

                  Vector3D eleOrigin;
                  eleOrigin.X = (EleArr + ele - 1)->G.Origin.X;
                  eleOrigin.Y = (EleArr + ele - 1)->G.Origin.Y;
                  eleOrigin.Z = (EleArr + ele - 1)->G.Origin.Z;
                  distele = (eleOrigin.X - xend) * (eleOrigin.X - xend) +
                            (eleOrigin.Y - yend) * (eleOrigin.Y - yend) +
                            (eleOrigin.Z - zend) * (eleOrigin.Z - zend);
                  distele = pow(distele, 0.5);

                  if (ele == ElementBgn[nearestprim]) {
                    mindistele = distele;
                    nearestele = ele;
                  }
                  if (distele < mindistele) {
                    mindistele = distele;
                    nearestele = ele;
                  }

                  if (debugFn) {
                    // printf("a, b, c, d, dist: %lg, %lg, %lg, %lg, %lg\n",
                    // a, b, c, d,  dist);
                    // printf("vector n: ai + bj + ck\n");
                    // printf("vector v: xgrd, ygrd, zgrd: %lg, %lg, %lg\n",
                    // xgrd, ygrd, zgrd);
                    printf(
                        "distele: %lg, mindist: %lg,  from nearest ele: %d\n",
                        distele, mindistele, nearestele);
                    fflush(stdout);
                  }

                  // }	// if PrimType is 3 or 4
                }  // for elements in nearestprim

                // if(mindistele <= 10.0e-6)
                // {
                EleIntsctd = nearestele;
                InEle = 1;
                ptintsct.X = (EleArr + EleIntsctd - 1)->G.Origin.X;
                ptintsct.Y = (EleArr + EleIntsctd - 1)->G.Origin.Y;
                ptintsct.Z = (EleArr + EleIntsctd - 1)->G.Origin.Z;
                NbChUpIonEle[EleIntsctd]++;

                fprintf(fPtIChUpMap, "%d %lg %lg %lg %d %d %d %d\n", ion,
                        ptintsct.X, ptintsct.Y, ptintsct.Z, PrimIntsctd, InPrim,
                        EleIntsctd, InEle);
                // }

                if (debugFn) {
                  printf("# ion: %d\n", ion);
                  printf("%lg %lg %lg\n", xlbend, ylbend, zlbend);
                  printf("%lg %lg %lg\n", xend, yend, zend);
                  printf("%lg, %lg, %lg\n", ptintsct.X, ptintsct.Y, ptintsct.Z);
                  printf("# Associated primitive: %d\n", PrimIntsctd);
                  printf("# Associated element and origin: %d, %lg, %lg, %lg\n",
                         EleIntsctd, (EleArr + EleIntsctd - 1)->G.Origin.X,
                         (EleArr + EleIntsctd - 1)->G.Origin.Y,
                         (EleArr + EleIntsctd - 1)->G.Origin.Z);
                  printf("#NbChUpIonEle on element: %d\n",
                         NbChUpIonEle[EleIntsctd]);
                  fprintf(ftmpIF, "#Element: %d\n", EleIntsctd);
                  polynode[0].X = (EleArr + EleIntsctd - 1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr + EleIntsctd - 1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr + EleIntsctd - 1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr + EleIntsctd - 1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr + EleIntsctd - 1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr + EleIntsctd - 1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr + EleIntsctd - 1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr + EleIntsctd - 1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr + EleIntsctd - 1)->G.Vertex[2].Z;
                  if (nvert == 4) {
                    polynode[3].X = (EleArr + EleIntsctd - 1)->G.Vertex[3].X;
                    polynode[3].Y = (EleArr + EleIntsctd - 1)->G.Vertex[3].Y;
                    polynode[3].Z = (EleArr + EleIntsctd - 1)->G.Vertex[3].Z;
                  }
                  fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X, polynode[0].Y,
                          polynode[0].Z);
                  fprintf(ftmpIF, "%lg %lg %lg\n", polynode[1].X, polynode[1].Y,
                          polynode[1].Z);
                  fprintf(ftmpIF, "%lg %lg %lg\n", polynode[2].X, polynode[2].Y,
                          polynode[2].Z);
                  if (nvert == 4) {
                    fprintf(ftmpIF, "%lg %lg %lg\n", polynode[3].X,
                            polynode[3].Y, polynode[3].Z);
                  }
                  fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X, polynode[0].Y,
                          polynode[0].Z);
                  fprintf(ftmpIF, "\n");
                  fflush(stdout);
                }  // debug
              }    // if prim == NbPrimitives

            }  // for all primitives // just not those on the volume

            if (debugFn)  // check ion positions, volume primitives and elements
            {
              char ionposdbg[256], inbstr[10];
              sprintf(inbstr, "%d", ion);
              strcpy(ionposdbg, "/tmp/Ion");
              strcat(ionposdbg, inbstr);
              strcat(ionposdbg, ".out");
              FILE *fipd = fopen(ionposdbg, "w");
              if (fipd == NULL) {
                printf(
                    "cannot open writable file to debug ion positions ...\n");
                printf("returning ...\n");
                return -111;
              }
              // write electron number, end, lbend, volume, primitive, elements,
              // intxn
              fprintf(fipd, "#ion: %d %d\n", inb, ion);  // should print same
              fprintf(fipd, "#last but end position:\n");
              fprintf(fipd, "%lg %lg %lg\n", xlbend, ylbend, zlbend);
              fprintf(fipd, "#end position:\n");
              fprintf(fipd, "%lg %lg %lg\n\n", xend, yend, zend);

              fprintf(fipd, "#intersected primitive number: %d\n", PrimIntsctd);
              if (PrimIntsctd >= 1) {
                fprintf(fipd, "#PrimType: %d\n", PrimType[PrimIntsctd]);
                fprintf(fipd, "#prim vertices:\n");
                fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][0],
                        YVertex[PrimIntsctd][0], ZVertex[PrimIntsctd][0]);
                fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][1],
                        YVertex[PrimIntsctd][1], ZVertex[PrimIntsctd][1]);
                fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][2],
                        YVertex[PrimIntsctd][2], ZVertex[PrimIntsctd][2]);
                if (PrimType[PrimIntsctd] == 4) {
                  fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][3],
                          YVertex[PrimIntsctd][3], ZVertex[PrimIntsctd][3]);
                }
                fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][0],
                        YVertex[PrimIntsctd][0], ZVertex[PrimIntsctd][0]);
                fprintf(fipd, "\n");

                fprintf(fipd, "#ptintsct:\n");
                fprintf(fipd, "%lg %lg %lg\n", ptintsct.X, ptintsct.Y,
                        ptintsct.Z);
                fprintf(fipd, "\n");
              }

              fprintf(fipd, "#intersected element number: %d\n", EleIntsctd);
              if (EleIntsctd >= 1) {
                int gtype = (EleArr + EleIntsctd - 1)->G.Type;
                fprintf(fipd, "#EleType: %d\n", gtype);
                fprintf(fipd, "#element vertices:\n");
                double x0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].X;
                double y0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].Y;
                double z0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].Z;
                double x1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].X;
                double y1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].Y;
                double z1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].Z;
                double x2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].X;
                double y2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].Y;
                double z2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].Z;
                fprintf(fipd, "%lg %lg %lg\n", x0, y0, z0);
                fprintf(fipd, "%lg %lg %lg\n", x1, y1, z1);
                fprintf(fipd, "%lg %lg %lg\n", x2, y2, z2);
                if (gtype == 4) {
                  double x3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].X;
                  double y3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].Y;
                  double z3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].Z;
                  fprintf(fipd, "%lg %lg %lg\n", x3, y3, z3);
                }
                fprintf(fipd, "%lg %lg %lg\n", x0, y0, z0);
                fprintf(fipd, "\n");

                fprintf(fipd, "#ptintsct:\n");
                fprintf(fipd, "%lg %lg %lg\n", ptintsct.X, ptintsct.Y,
                        ptintsct.Z);
                fprintf(fipd, "\n");
              }
              fclose(fipd);
            }  // if 1
          }    // for all the ions
          fclose(fPtIChUpMap);

          // This file contains information about number of ions (I)
          // and total (E+I) charge deposition on each element
          FILE *fEleEIChUpMap = fopen("EleE+IChUpMap.out", "w");
          if (fEleEIChUpMap == NULL) {
            printf("cannot open EleE+IChUpMap.out file for writing ...\n");
            return 111;
          }
          for (int ele = 1; ele <= NbElements; ++ele) {
            (EleArr + ele - 1)->Assigned +=
                ChUpFactor * Q_I * NbChUpIonEle[ele] / (EleArr + ele - 1)->G.dA;
            fprintf(fEleEIChUpMap, "%d %lg %lg %lg %d %lg\n", ele,
                    (EleArr + ele - 1)->G.Origin.X,
                    (EleArr + ele - 1)->G.Origin.Y,
                    (EleArr + ele - 1)->G.Origin.Z, NbChUpIonEle[ele],
                    (EleArr + ele - 1)->Assigned);
          }
          fclose(fEleEIChUpMap);

          fclose(ftmpIF);
          free(NbChUpIonEle);
        }  // Calculation for ions ends

      }  // OptChargingUp

      fclose(ChargingUpInpFile);
    }  // else ChargingUpInpFile

    if (debugFn) {
      // print all the elements and their number of charging up e-s and i-s
    }
  }  // charging up parameters set up

  return (0);
}  // neBEMChargingUp ends

// There are several flags associated with this crucial step.
// These flags have a hieracrchy, the first mentioned having the highest
// priority.
// NewModel: 1 implies a fresh calculation.
// NewMesh: 1 implies a new mesh for a new / old device.
// NewBC: 1 implies new RHS for the same LHS; skips matrix inversion.
// NewPP: 1 implies the use of the same solution;
//             skips matrix inversion, as well as, computing the solution.
// It should be noted that a given device can be modeled in different manners.
// The same model for a device can have different discretization,
// same mesh different boundary conditions leading to different solutions
// and the same solution used to carry out different post-processes, we maintain
// four counters as well.
// ModelCntr: keeps track of the model for a given device,
// MeshCntr: keeps track of the mesh for a given model
// BCCntr: keeps track of the association of the boundary condition and
//               its solution. This has to maintained by the user manually and
//               supplied, for example, while carrying out a post-processing
//               for a solution that was computed before.
// PPCntr: numbers different post-processes for a given solution resulting from
//         a given set of preceding conditions.
int neBEMSolve(void) {
  int dbgFn = 0;

  clock_t startSolveClock = clock();

  if (TimeStep < 1) {
    neBEMMessage("neBEMSolve - TimeStep cannot be less than one!;\n");
    neBEMMessage("             Please put TimeStep = 1 for static problems.\n");
  }

  if ((neBEMState == 5) || (neBEMState == 8)) {
    if (neBEMState == 8)  // neBEMState 8 must have inverted flag on
    {  // so it must be a case looking for solution with a NewBC
      if (NewBC == 0) {
        neBEMMessage("neBEMSolve - NewBC zero for neBEMState = 8!");
        neBEMMessage("           - Nothing to be done ... returning.");
        return -1;
      }
    }

    int fstatus;
    if (NewModel) {  // effectively, NewMesh = NewBC = NewPP = 1;
      fstatus = ComputeSolution();
      if (fstatus != 0) {
        neBEMMessage("neBEMSolve - NewModel");
        return -1;
      }
    } else { // NewModel == 0
      if (NewMesh) {
        // effectively, NewBC = NewPP = 1;
        fstatus = ComputeSolution();
        if (fstatus != 0) {
          neBEMMessage("neBEMSolve - NewMesh");
          return -1;
        }
      } else { // NewModel == NewMesh == 0
        if (NewBC) {  // effectively, NewPP = 1;
          fstatus = ComputeSolution();
          if (fstatus != 0) {
            neBEMMessage("neBEMSolve - Failure computing new solution");
            return -1;
          }
        } else { // NewBC == 0
          if (NewPP) {
            fstatus = ReadSolution();
            if (fstatus != 0) {
              neBEMMessage("neBEMSolve - Failure reading solution");
              return (-1);
            }
          } else { // NewPP == 0
            printf("neBEMSolve: Nothing to do ... returning ...\n");
            return (-1);
          }  // NewPP == 0
        }    // NewBC == 0
      }      // NewModel == NewDiscretization == 0
    }        // NewModel == 0

    neBEMState = 9;
  } else {
    printf("neBEMSolve: neBEMSolve can be called only in state 5 / 8 ...\n");
    printf("returning ...\n");
    return (-1);
  }

  if (FailureCntr) {
    printf(
        "neBEMSolve: Approximations were made while computing the influence "
        "coefficients.\n");
    printf("            Please check the \"%s/Isles.log\" file.\n", PPOutDir);
  }

  clock_t stopSolveClock = clock();
  neBEMTimeElapsed(startSolveClock, stopSolveClock);
  printf("to complete solve\n");

  // Prepare voxelized data that will be exported to Garfield++
  if (OptVoxel) {
    clock_t startVoxelClock = clock();

    int fstatus = VoxelFPR();
    if (fstatus != 0) {
      neBEMMessage("neBEMSolve - Failure computing VoxelFPR");
      return -1;
    }

    clock_t stopVoxelClock = clock();
    neBEMTimeElapsed(startVoxelClock, stopVoxelClock);
    printf("to compute VoxelFPR\n");
  }

  // Prepare 3dMap data that will be exported to Garfield++
  if (OptMap) {
    clock_t startMapClock = clock();

    int fstatus = MapFPR();
    if (fstatus != 0) {
      neBEMMessage("neBEMSolve - Failure computing MapFPR");
      return -1;
    }

    clock_t stopMapClock = clock();
    neBEMTimeElapsed(startMapClock, stopMapClock);
    printf("to compute MapFPR\n");
  }

  // allocate memory for potential and field components within FAST volume mesh
  // and compute / read FastVol data
  // Similar allocation, computation and reading may be necessary for the KnCh
  // effects.
  // The other approach could be to create fast volumes that always have both
  // the influences (elements + known charges) added together. This approach
  // seems more managable now and is being followed.
  // Thus, if we want to inspect the effects of elements and known charges
  // separately, we will have to generate one fast volume with OptKnCh = 0,
  // and another with OptKnCh = 1. Subtraction of these two fast volumes will
  // provide us with the effect of KnCh.
  if (OptFastVol) {
    int MaxXCells = BlkNbXCells[1];
    int MaxYCells = BlkNbYCells[1];
    int MaxZCells = BlkNbZCells[1];
    clock_t startFastClock = clock();

    // find maximum number of Xcells etc in all the blocks
    // simplifies memory allocation using nrutils but hogs memory!
    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      if (block == 1) {
        MaxXCells = BlkNbXCells[1];
        MaxYCells = BlkNbYCells[1];
        MaxZCells = BlkNbZCells[1];
      } else {
        if (MaxXCells < BlkNbXCells[block]) MaxXCells = BlkNbXCells[block];
        if (MaxYCells < BlkNbYCells[block]) MaxYCells = BlkNbYCells[block];
        if (MaxZCells < BlkNbZCells[block]) MaxZCells = BlkNbZCells[block];
      }
    }  // loop block for finding maxm cells among all the blocks

    if (dbgFn) {
      printf("OptFastVol: %d\n", OptFastVol);
      printf("NbPtSkip: %d\n", NbPtSkip);
      printf("OptStaggerFastVol: %d\n", OptStaggerFastVol);
      printf("NbStgPtSkip: %d\n", NbStgPtSkip);
      printf("OptReadFastPF: %d\n", OptReadFastPF);
      printf("LX: %le\n", FastVol.LX);
      printf("LY: %le\n", FastVol.LY);
      printf("LZ: %le\n", FastVol.LZ);
      printf("CornerX: %le\n", FastVol.CrnrX);
      printf("CornerY: %le\n", FastVol.CrnrY);
      printf("CornerZ: %le\n", FastVol.CrnrZ);
      printf("YStagger: %le\n", FastVol.YStagger);
      printf("NbOfBlocks: %d\n", FastVol.NbBlocks);
      for (int block = 1; block <= FastVol.NbBlocks; ++block) {
        printf("NbOfXCells[%d]: %d\n", block, BlkNbXCells[block]);
        printf("NbOfYCells[%d]: %d\n", block, BlkNbYCells[block]);
        printf("NbOfZCells[%d]: %d\n", block, BlkNbZCells[block]);
        printf("LZ[%d]: %le\n", block, BlkLZ[block]);
        printf("CornerZ[%d]: %le\n", block, BlkCrnrZ[block]);
      }
      printf("NbOfOmitVols: %d\n", FastVol.NbOmitVols);
      if (FastVol.NbOmitVols) {
        for (int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
          printf("OmitVolLX[%d]: %le\n", omit, OmitVolLX[omit]);
          printf("OmitVolLY[%d]: %le\n", omit, OmitVolLY[omit]);
          printf("OmitVolLZ[%d]: %le\n", omit, OmitVolLZ[omit]);
          printf("OmitCrnrX[%d]: %le\n", omit, OmitVolCrnrX[omit]);
          printf("OmitCrnrY[%d]: %le\n", omit, OmitVolCrnrY[omit]);
          printf("OmitCrnrZ[%d]: %le\n", omit, OmitVolCrnrZ[omit]);
        }
      }
      printf("MaxXCells, MaxYCells, MaxZCells: %d, %d, %d\n", MaxXCells,
             MaxYCells, MaxZCells);
    }  // dbgFn

    /* Memory wastage!!! Improve as soon as possible. */
    FastPot = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1, MaxYCells + 1,
                       1, MaxZCells + 1);
    FastFX = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1, MaxYCells + 1,
                      1, MaxZCells + 1);
    FastFY = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1, MaxYCells + 1,
                      1, MaxZCells + 1);
    FastFZ = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1, MaxYCells + 1,
                      1, MaxZCells + 1);

    if (OptStaggerFastVol) {
      /* Memory wastage!!! Improve as soon as possible. */
      FastStgPot = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1,
                            MaxYCells + 1, 1, MaxZCells + 1);
      FastStgFX = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);
      FastStgFY = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);
      FastStgFZ = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);
    }  // if OptStaggerFastVol

    if ((OptCreateFastPF) && (!OptReadFastPF))  // reading overrides creation
    {
      int fstatus = FastVolPF();
      if (fstatus != 0) {
        neBEMMessage("neBEMSolve - Failure computing FastVolPF");
        return -1;
      }
    }  // if OptCreateFastPF

    if (OptReadFastPF) {
      int nbXCells, nbYCells, nbZCells;
      int tmpblk;
      double xpt, ypt, zpt;

      char FastVolPFFile[256];
      strcpy(FastVolPFFile, BCOutDir);
      strcat(FastVolPFFile, "/FastVolPF.out");
      FILE* fFastVolPF = fopen(FastVolPFFile, "r");
      if (fFastVolPF == NULL) {
        neBEMMessage("in neBEMSolve - FastVolPFFile");
        return -1;
      }

      fscanf(fFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

      for (int block = 1; block <= FastVol.NbBlocks; ++block) {
        nbXCells = BlkNbXCells[block];
        nbYCells = BlkNbYCells[block];
        nbZCells = BlkNbZCells[block];

        for (int i = 1; i <= nbXCells + 1; ++i) {
          for (int j = 1; j <= nbYCells + 1; ++j) {
            for (int k = 1; k <= nbZCells + 1; ++k) {
              fscanf(fFastVolPF, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                     &tmpblk, &xpt, &ypt, &zpt, &FastPot[block][i][j][k],
                     &FastFX[block][i][j][k], &FastFY[block][i][j][k],
                     &FastFZ[block][i][j][k]);
            }  // loop k
          }    // loop j
        }      // loop i
      }        // loop block
      fclose(fFastVolPF);

      if (OptStaggerFastVol) {
        char FastStgVolPFFile[256];
        strcpy(FastStgVolPFFile, BCOutDir);
        strcat(FastStgVolPFFile, "/FastStgVolPF.out");
        FILE *fFastStgVolPF = fopen(FastStgVolPFFile, "r");
        if (fFastStgVolPF == NULL) {
          neBEMMessage("in neBEMSolve - FastStgVolPFFile");
          return -1;
        }

        fscanf(fFastStgVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

        for (int block = 1; block <= FastVol.NbBlocks; ++block) {
          nbXCells = BlkNbXCells[block];
          nbYCells = BlkNbYCells[block];
          nbZCells = BlkNbZCells[block];

          for (int i = 1; i <= nbXCells + 1; ++i) {
            for (int j = 1; j <= nbYCells + 1; ++j) {
              for (int k = 1; k <= nbZCells + 1; ++k) {
                fscanf(fFastStgVolPF, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                       &tmpblk, &xpt, &ypt, &zpt, &FastStgPot[block][i][j][k],
                       &FastStgFX[block][i][j][k], &FastStgFY[block][i][j][k],
                       &FastStgFZ[block][i][j][k]);
              }  // loop k
            }    // loop j
          }      // loop i
        }        // loop block
        fclose(fFastStgVolPF);
      }  // if OptStaggerFastVol
    }    // if OptReadFastPF

    clock_t stopFastClock = clock();
    neBEMTimeElapsed(startFastClock, stopFastClock);
    printf("to compute / read FastVolPF\n");
  }  // if OptFastVol

  if (OptWtFldFastVol)  // allocate memory for weighting field fast volume
                        // variables
  {
    int MaxXCells = WtFldBlkNbXCells[1];
    int MaxYCells = WtFldBlkNbYCells[1];
    int MaxZCells = WtFldBlkNbZCells[1];
    clock_t startFastClock = clock();

    // find maximum number of Xcells etc in all the blocks
    // simplifies memory allocation using nrutils but hogs memory!
    for (int block = 1; block <= WtFldFastVol.NbBlocks; ++block) {
      if (block == 1) {
        MaxXCells = WtFldBlkNbXCells[1];
        MaxYCells = WtFldBlkNbYCells[1];
        MaxZCells = WtFldBlkNbZCells[1];
      } else {
        if (MaxXCells < WtFldBlkNbXCells[block])
          MaxXCells = WtFldBlkNbXCells[block];
        if (MaxYCells < WtFldBlkNbYCells[block])
          MaxYCells = WtFldBlkNbYCells[block];
        if (MaxZCells < WtFldBlkNbZCells[block])
          MaxZCells = WtFldBlkNbZCells[block];
      }
    }  // loop block for finding maxm cells among all the blocks

    if (dbgFn) {
      printf("OptWtFldFastVol: %d\n", OptWtFldFastVol);
      printf("WtFldNbPtSkip: %d\n", WtFldNbPtSkip);
      printf("OptWtFldStaggerFastVol: %d\n", OptWtFldStaggerFastVol);
      printf("WtFldNbStgPtSkip: %d\n", WtFldNbStgPtSkip);
      printf("OptWtFldReadFastPF: %d\n", OptWtFldReadFastPF);
      printf("LX: %le\n", WtFldFastVol.LX);
      printf("LY: %le\n", WtFldFastVol.LY);
      printf("LZ: %le\n", WtFldFastVol.LZ);
      printf("CornerX: %le\n", WtFldFastVol.CrnrX);
      printf("CornerY: %le\n", WtFldFastVol.CrnrY);
      printf("CornerZ: %le\n", WtFldFastVol.CrnrZ);
      printf("YStagger: %le\n", WtFldFastVol.YStagger);
      printf("NbOfBlocks: %d\n", WtFldFastVol.NbBlocks);
      for (int block = 1; block <= WtFldFastVol.NbBlocks; ++block) {
        printf("NbOfXCells[%d]: %d\n", block, WtFldBlkNbXCells[block]);
        printf("NbOfYCells[%d]: %d\n", block, WtFldBlkNbYCells[block]);
        printf("NbOfZCells[%d]: %d\n", block, WtFldBlkNbZCells[block]);
        printf("LZ[%d]: %le\n", block, WtFldBlkLZ[block]);
        printf("CornerZ[%d]: %le\n", block, WtFldBlkCrnrZ[block]);
      }
      printf("NbOfOmitVols: %d\n", WtFldFastVol.NbOmitVols);
      if (WtFldFastVol.NbOmitVols) {
        for (int omit = 1; omit <= WtFldFastVol.NbOmitVols; ++omit) {
          printf("OmitVolLX[%d]: %le\n", omit, WtFldOmitVolLX[omit]);
          printf("OmitVolLY[%d]: %le\n", omit, WtFldOmitVolLY[omit]);
          printf("OmitVolLZ[%d]: %le\n", omit, WtFldOmitVolLZ[omit]);
          printf("OmitCrnrX[%d]: %le\n", omit, WtFldOmitVolCrnrX[omit]);
          printf("OmitCrnrY[%d]: %le\n", omit, WtFldOmitVolCrnrY[omit]);
          printf("OmitCrnrZ[%d]: %le\n", omit, WtFldOmitVolCrnrZ[omit]);
        }
      }
      printf("MaxXCells, MaxYCells, MaxZCells: %d, %d, %d\n", MaxXCells,
             MaxYCells, MaxZCells);
    }  // dbgFn

    /* Memory wastage!!! Improve as soon as possible. */
    WtFldFastPot = d4tensor(1, WtFldFastVol.NbBlocks, 1, MaxXCells + 1, 1,
                            MaxYCells + 1, 1, MaxZCells + 1);
    WtFldFastFX = d4tensor(1, WtFldFastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);
    WtFldFastFY = d4tensor(1, WtFldFastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);
    WtFldFastFZ = d4tensor(1, WtFldFastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);

    if (OptWtFldStaggerFastVol) {
      /* Memory wastage!!! Improve as soon as possible. */
      WtFldFastStgPot = d4tensor(1, WtFldFastVol.NbBlocks, 1, MaxXCells + 1, 1,
                                 MaxYCells + 1, 1, MaxZCells + 1);
      WtFldFastStgFX = d4tensor(1, WtFldFastVol.NbBlocks, 1, MaxXCells + 1, 1,
                                MaxYCells + 1, 1, MaxZCells + 1);
      WtFldFastStgFY = d4tensor(1, WtFldFastVol.NbBlocks, 1, MaxXCells + 1, 1,
                                MaxYCells + 1, 1, MaxZCells + 1);
      WtFldFastStgFZ = d4tensor(1, WtFldFastVol.NbBlocks, 1, MaxXCells + 1, 1,
                                MaxYCells + 1, 1, MaxZCells + 1);
    }  // if OptWtFldStaggerFastVol

    if ((OptWtFldCreateFastPF) && (!OptWtFldReadFastPF))  // reading overrides
    {
      // Computing weighting field fast volume has not been implemented
      neBEMMessage(
          "neBEMSolve - Failure computing WtFldFastVolPF: not implemented");
      return -1;
    }  // if OptWtFldCreateFastPF

    if (OptWtFldReadFastPF)  // reading option overrides creation
    {
      int nbXCells, nbYCells, nbZCells;
      int tmpblk;
      double xpt, ypt, zpt;

      char FastVolPFFile[256];
      strcpy(FastVolPFFile, BCOutDir);
      strcat(FastVolPFFile, "/WtFldFastVolPF.out");
      FILE *fFastVolPF = fopen(FastVolPFFile, "r");
      if (fFastVolPF == NULL) {
        neBEMMessage("in neBEMSolve - WtFldFastVolPFFile");
        return -1;
      }

      fscanf(fFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

      for (int block = 1; block <= WtFldFastVol.NbBlocks; ++block) {
        nbXCells = WtFldBlkNbXCells[block];
        nbYCells = WtFldBlkNbYCells[block];
        nbZCells = WtFldBlkNbZCells[block];

        for (int i = 1; i <= nbXCells + 1; ++i) {
          for (int j = 1; j <= nbYCells + 1; ++j) {
            for (int k = 1; k <= nbZCells + 1; ++k) {
              fscanf(fFastVolPF, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                     &tmpblk, &xpt, &ypt, &zpt, &WtFldFastPot[block][i][j][k],
                     &WtFldFastFX[block][i][j][k], &WtFldFastFY[block][i][j][k],
                     &WtFldFastFZ[block][i][j][k]);
            }  // loop k
          }    // loop j
        }      // loop i
      }        // loop block
      fclose(fFastVolPF);

      if (OptWtFldStaggerFastVol) {
        char FastStgVolPFFile[256];
        FILE *fFastStgVolPF;
        strcpy(FastStgVolPFFile, BCOutDir);
        strcat(FastStgVolPFFile, "/WtFldFastStgVolPF.out");
        fFastStgVolPF = fopen(FastStgVolPFFile, "r");

        if (fFastStgVolPF == NULL) {
          neBEMMessage("in neBEMSolve - WtFldFastStgVolPFFile");
          return -1;
        }

        fscanf(fFastStgVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

        for (int block = 1; block <= WtFldFastVol.NbBlocks; ++block) {
          nbXCells = WtFldBlkNbXCells[block];
          nbYCells = WtFldBlkNbYCells[block];
          nbZCells = WtFldBlkNbZCells[block];

          for (int i = 1; i <= nbXCells + 1; ++i) {
            for (int j = 1; j <= nbYCells + 1; ++j) {
              for (int k = 1; k <= nbZCells + 1; ++k) {
                fscanf(fFastStgVolPF, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                       &tmpblk, &xpt, &ypt, &zpt,
                       &WtFldFastStgPot[block][i][j][k],
                       &WtFldFastStgFX[block][i][j][k],
                       &WtFldFastStgFY[block][i][j][k],
                       &WtFldFastStgFZ[block][i][j][k]);
              }  // loop k
            }    // loop j
          }      // loop i
        }        // loop block
        fclose(fFastStgVolPF);
      }  // if OptWtFldStaggerFastVol
    }    // if OptWtFldReadFastPF

    clock_t stopFastClock = clock();
    neBEMTimeElapsed(startFastClock, stopFastClock);
    printf("to compute / read FastVolPF\n");
  }  // if OptWtFldFastVol

  return (0);
}  // neBEMSolve ends

// Get potential and field at a given point
int neBEMField(Point3D *point, double *potential, Vector3D *field) {
  if (neBEMState < 9) {
    printf("neBEMField cannot be called before reaching state 9.\n");
    return (-1);
  }

  // printf("neBEMField called %8d times", ++neBEMFieldCallCntr);

  double Pot;
  int fstatus;
  if (OptFastVol)  // Note: this is not the Create or Read option
  {
    fstatus = FastPFAtPoint(point, &Pot, field);
    if (fstatus != 0) {
      neBEMMessage("neBEMField - FastPFAtPoint");
      return -1;
    }
  } else {
    fstatus = PFAtPoint(point, &Pot, field);
    if (fstatus != 0) {
      neBEMMessage("neBEMField - PFAtPoint");
      return -1;
    }
  }

  *potential = Pot;

  // printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");

  return (0);
}  // neBEMField ends

// Actual preparation of the weighting field.
// The return value identifies the weighting field. Error: id < 0.
// The list contains all primitives that are part of this particular
// read-out group. These can come from several volumes, but it is
// not anticipated that only some of the primitives of one volume
// are listed.
// Weighitng field boundary conditions are useful only when the inverted
// matrix is available.
// This state is assigned either after element discretization has been
// completed or in a condition when we are looking for modifying only the
// boundary condition for a device having same geometry (hence, the same
// inverted influence coefficient matrix)
int neBEMPrepareWeightingField(int nprim, int primlist[]) {
  static int IdWtField = 0;

  if (neBEMState < 7) {
    printf(
        "neBEMPrepareWeightingField: Weighting computations only meaningful "
        "beyond neBEMState 7 ...\n");
    return -1;
  }

  // Find first free slot
  const int MaxWtField = 100; // used again while deallocating these memories
  if (WtFieldChDen == NULL)
    WtFieldChDen = (double **)malloc(MaxWtField * sizeof(double *));
  if (AvWtChDen == NULL)
    AvWtChDen = (double **)malloc(MaxWtField * sizeof(double *));

  ++IdWtField;
  if (IdWtField >= MaxWtField) {
    printf(
        "neBEMPrepareWeightingField: reached MaxWtField weighting fields.\n");
    return -1;
  }

  // Allocate a new column to store this solution set
  WtFieldChDen[IdWtField] = (double *)malloc((NbElements + 2) * sizeof(double));
  AvWtChDen[IdWtField] = (double *)malloc((NbPrimitives + 2) * sizeof(double));

  int fstatus =
      WeightingFieldSolution(nprim, primlist, WtFieldChDen[IdWtField]);
  if (!fstatus)  // estimate primitive related average wt field charge densities
  {
    // OMPCheck - may be parallelized
    for (int prim = 1; prim <= NbPrimitives; ++prim) {
      double area = 0.0;  // need area of the primitive as well!
      AvWtChDen[IdWtField][prim] = 0.0;

      for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim]; ++ele) {
        area += (EleArr + ele - 1)->G.dA;
        AvWtChDen[IdWtField][prim] +=
            WtFieldChDen[IdWtField][ele] * (EleArr + ele - 1)->G.dA;
      }

      AvWtChDen[IdWtField][prim] /= area;
    }
  }
  if (fstatus != 0) {
    neBEMMessage("neBEMPrepareWeightingField - WeightingFieldSolution");
    return -1;
  }

  return IdWtField;
}  // neBEMPrepareWeightingField ends

// Deallocates memory reserved for a weighting field
void neBEMDeleteWeightingField(int IdWtField) { 
  free(WtFieldChDen[IdWtField]); 
  free(AvWtChDen[IdWtField]);
}

// Deallocates all memory reserved for all weighting fields
void neBEMDeleteAllWeightingFields(void) {
  const int MaxWtField = 100;	// being used while allocating memory
  for (int id = 1; id < MaxWtField; ++id)	{ // count from 1
    free(WtFieldChDen[id]);
    free(AvWtChDen[id]);
  }
  free(WtFieldChDen); 
  free(AvWtChDen);
}

// Get weighting field (potential also) at a specific point
// returns DBL_MAX as the value of potential when something goes wrong.
double neBEMWeightingField(Point3D *point, Vector3D *field, int IdWtField) {
  double potential;

  if (neBEMState < 9) {
    printf("neBEMWeightingField cannot be called before reaching state 9.\n");
    return (-1);
  }

  if (OptFixedWtField) { 
    // minimum computation, too restricted!
    potential = FixedWtPotential;
    field->X = FixedWtFieldX;
    field->Y = FixedWtFieldY;
    field->Z = FixedWtFieldZ;
  } else if (OptWtFldFastVol) {
    // bit more computation, lot more flexibility
    // Note: this is not the Creat or Read option
    int fstatus = WtFldFastPFAtPoint(point, &potential, field);
    if (fstatus != 0) {
      neBEMMessage("neBEMWeightingField - WtFldFastPFAtPoint");
      return DBL_MAX;
    }
  } else {
    int fstatus = WtPFAtPoint(point, &potential, field, IdWtField);
    if (fstatus != 0) {
      neBEMMessage("neBEMWeightingField - WtPFAtPoint");
      return DBL_MAX;
    }
  }

  return potential;
}  // neBEMWeightingField ends

double neBEMVolumeCharge(int volume) {
  // Initialise the sum
  double sumcharge = 0.0;

  // Loop over all elements
  for (int elem = 1; elem <= NbElements; ++elem) {
    // Find the primitive number for the element
    int prim = (EleArr + elem - 1)->PrimitiveNb;

    // Find out to which volume this belongs to
    int volref = VolRef1[prim];

    // Skip the rest if the volume is not right
    if (volref + 1 != volume) {
      continue;
    }

    // Calculate the periodic weight of the primitive
    int rptCnt = 0;
    if (PeriodicInX[prim] || PeriodicInY[prim] || PeriodicInZ[prim]) {
      for (int xrpt = -PeriodicInX[prim]; xrpt <= PeriodicInX[prim]; ++xrpt)
        for (int yrpt = -PeriodicInY[prim]; yrpt <= PeriodicInY[prim]; ++yrpt)
          for (int zrpt = -PeriodicInZ[prim]; zrpt <= PeriodicInZ[prim];
               ++zrpt) {
            if ((xrpt == 0) && (yrpt == 0) && (zrpt == 0))
              continue;
            else
              ++rptCnt;
          }
    } else {
      rptCnt = 1;
    }

    // Add the charge
    // printf("Element: %d, volume: %d, charge: %g\n", elem, volref,
    //    (EleArr+elem-1)->Solution * (EleArr+elem-1)->G.dA);
    sumcharge +=
        rptCnt * (EleArr + elem - 1)->Solution * (EleArr + elem - 1)->G.dA;
  }  // loop over elements

  // Return the result
  // printf("Charge on volume %d: %g C\n", volume, sumcharge);
  return sumcharge;
}  // end of neBEMVolumeCharge

int neBEMEnd(void) {
  fprintf(fIsles,
          "IslesCntr: %d, ExactCntr: %d, FailureCntr: %d, ApproxCntr: %d\n",
          IslesCntr, ExactCntr, FailureCntr, ApproxCntr);
  fclose(fIsles);
  fIsles = NULL;
  printf("neBEM ends ... bye!\n");

  return 0;
}  // neBEMEnd ends

// In a given problem, the DeviceOutDir is the one that is unique for a given
// device. Within it, there are several sub-directories, one related to each
// device counter; within each device counter directory, there can be
// several sub-directories related to each Mesh specification;
// within each mesh sub-dir, there can be several sub-directories, one for each
// set of boundary conditions and finally, for each boundary condition,
// several related to each post-processing counter.
int CreateDirStr(void) {
  char strModelCntr[10], strMeshCntr[10], strBCCntr[10], strPPCntr[10];
  int CreateOrUseDir(char[]);
  int CreateDirOrQuit(char[]);

  sprintf(strModelCntr, "/Model%d", ModelCntr);
  sprintf(strMeshCntr, "/M%d", MeshCntr);
  sprintf(strBCCntr, "/BC%d", BCCntr);
  sprintf(strPPCntr, "/PP%d", PPCntr);

  strcpy(ModelOutDir, DeviceOutDir);
  strcat(ModelOutDir, strModelCntr);
  strcpy(NativeOutDir, ModelOutDir);
  strcat(NativeOutDir, "/neBEMNatives/");
  strcpy(NativePrimDir, NativeOutDir);
  strcat(NativePrimDir, "Primitives/");
  strcpy(MeshOutDir, ModelOutDir);
  strcat(MeshOutDir, strMeshCntr);
  strcpy(BCOutDir, MeshOutDir);
  strcat(BCOutDir, strBCCntr);
  strcpy(PPOutDir, BCOutDir);
  strcat(PPOutDir, strPPCntr);

  // Create DeviceOutDir, if necessary
  int fstatus = CreateOrUseDir(DeviceOutDir);
  if (fstatus != 0) {
    neBEMMessage("CreateDirStr - CreateOrUseDir");
    return -1;
  }

  if (NewModel) {
    // create ModelOutDir
    if (OptReuseDir) {
      fstatus = CreateOrUseDir(ModelOutDir);
      fstatus = CreateOrUseDir(NativeOutDir);
      fstatus = CreateOrUseDir(NativePrimDir);
    } else {
      fstatus = CreateDirOrQuit(ModelOutDir);
      fstatus = CreateDirOrQuit(NativeOutDir);
      fstatus = CreateDirOrQuit(NativePrimDir);
    }
    if (fstatus != 0) {
      neBEMMessage("CreateDirStr - ModelOutDir");
      return -1;
    }
  }

  if (NewMesh) {
    // create MeshOutDir
    if (OptReuseDir)
      fstatus = CreateOrUseDir(MeshOutDir);
    else
      fstatus = CreateDirOrQuit(MeshOutDir);
    if (fstatus != 0) {
      neBEMMessage("CreateDirStr - MeshOutDir");
      return -1;
    }
  }

  if (NewBC) {
    // create BCOutDir
    if (OptReuseDir)
      fstatus = CreateOrUseDir(BCOutDir);
    else
      fstatus = CreateDirOrQuit(BCOutDir);
    if (fstatus != 0) {
      neBEMMessage("CreateDirStr - BCOutDir");
      return -1;
    }
  }

  if (NewPP) {
    // create PPOutDir
    if (OptReuseDir)
      fstatus = CreateOrUseDir(PPOutDir);
    else
      fstatus = CreateDirOrQuit(PPOutDir);
    if (fstatus != 0) {
      neBEMMessage("CreateDirStr - PPOutDir");
      return -1;
    }
  }

  // Create other relevant sub-directories
  char subdir[256];

  strcpy(subdir, ModelOutDir);
  strcat(subdir, "/Primitives/");
  if (OptReuseDir)
    fstatus = CreateOrUseDir(subdir);
  else
    fstatus = CreateDirOrQuit(subdir);

  strcpy(subdir, MeshOutDir);
  strcat(subdir, "/Elements/");
  if (OptReuseDir)
    fstatus = CreateOrUseDir(subdir);
  else
    fstatus = CreateDirOrQuit(subdir);

  strcpy(subdir, MeshOutDir);
  strcat(subdir, "/GViewDir/");
  if (OptReuseDir)
    fstatus = CreateOrUseDir(subdir);
  else
    fstatus = CreateDirOrQuit(subdir);

  return (0);
}  // CreateDirStr ends

// Create or use directory dirname
int CreateOrUseDir(char dirname[]) {
  char dirstr[256];
  struct stat st;

  if (stat(dirname, &st) == 0)  // feel safe to use an existing directory
  {
    printf("Previous %s exists ... using the existing directory ... \n",
           dirname);
  } else {
    sprintf(dirstr, "mkdir -p %s", dirname);
    if (system(dirstr))  // returns 0 if successful
    {
      printf("Cannot create dirname %s ... returning ...\n", dirname);
      return (-1);
    }
  }

  return (0);
}  // CreateOrUseDir ends

// Create directory dirname or quit reporting failure
int CreateDirOrQuit(char dirname[]) {
  char dirstr[256];
  struct stat st;

  if (stat(dirname, &st) == 0)  // not safe to use an existing directory
  {
    printf("Previous %s exists ... please check inputs and counters ... \n",
           dirname);
    return (-1);
  } else {
    sprintf(dirstr, "mkdir -p %s", dirname);
    if (system(dirstr))  // returns 0 if successful
    {
      printf("Cannot create dirname %s ... returning ...\n", dirname);
      return (-1);
    }
  }

  return (0);
}  // CreateOrQuitDir ends

int neBEMMessage(const char *message) {
  fprintf(stdout, "neBEMMessage: %s\n", message);

  return 0;
}  // neBEMMessage ends

int WritePrimitives(void) {
  char PrimitiveFile[256];

  strcpy(PrimitiveFile, ModelOutDir);
  strcat(PrimitiveFile, "/Primitives/StorePrims.out");

  FILE *fStrPrm = fopen(PrimitiveFile, "w");
  if (fStrPrm == NULL) {
    neBEMMessage("WritePrimitives - Could not create file to store primitives");
    return -1;
  }

  fprintf(fStrPrm, "%d %d\n", NbVolumes, VolMax);
  fprintf(fStrPrm, "%d\n", NbPrimitives);
  fprintf(fStrPrm, "%d\n", MaxNbVertices);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    fprintf(fStrPrm, "%d\n", PrimType[prim]);
    fprintf(fStrPrm, "%d\n", InterfaceType[prim]);
    fprintf(fStrPrm, "%d\n", NbVertices[prim]);

    for (int vert = 0; vert < NbVertices[prim]; ++vert) {
      fprintf(fStrPrm, "%le %le %le\n", XVertex[prim][vert],
              YVertex[prim][vert], ZVertex[prim][vert]);
    }  // vert loop

    fprintf(fStrPrm, "%le %le %le\n", XNorm[prim], YNorm[prim], ZNorm[prim]);
    fprintf(fStrPrm, "%le\n", Radius[prim]);

    fprintf(fStrPrm, "%le %le %le %le %le\n", Epsilon1[prim], Epsilon2[prim],
            Lambda[prim], ApplPot[prim], ApplCh[prim]);

    fprintf(fStrPrm, "%d %d\n", VolRef1[prim], VolRef2[prim]);

    fprintf(fStrPrm, "%d %d %d\n", PeriodicTypeX[prim], PeriodicTypeY[prim],
            PeriodicTypeZ[prim]);
    fprintf(fStrPrm, "%d %d %d\n", PeriodicInX[prim], PeriodicInY[prim],
            PeriodicInZ[prim]);
    fprintf(fStrPrm, "%le %le %le\n", XPeriod[prim], YPeriod[prim],
            ZPeriod[prim]);
    fprintf(fStrPrm, "%le %le %le\n", MirrorDistXFromOrigin[prim],
            MirrorDistYFromOrigin[prim], MirrorDistZFromOrigin[prim]);
  }  // prim loop

  fclose(fStrPrm);

  return 0;
}  // WritePrimitives ends

int WriteElements(void) {
  char ElementFile[256];

  strcpy(ElementFile, MeshOutDir);
  strcat(ElementFile, "/Elements/StoreElems.out");

  FILE *fStrEle = fopen(ElementFile, "w");
  if (fStrEle == NULL) {
    neBEMMessage("WriteElements - Could not create file to store elements");
    return -1;
  }

  fprintf(fStrEle, "%d %d\n", NbSurfs, NbWires);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    fprintf(fStrEle, "%d %d\n", NbSurfSegX[prim], NbSurfSegZ[prim]);
    fprintf(fStrEle, "%d\n", NbWireSeg[prim]);
  }

  fprintf(fStrEle, "%d\n", NbElements);

  for (int ele = 1; ele <= NbElements; ++ele) {
    fprintf(fStrEle, "%d %d %d %d %d\n", (EleArr + ele - 1)->DeviceNb,
            (EleArr + ele - 1)->ComponentNb, (EleArr + ele - 1)->PrimitiveNb,
            (EleArr + ele - 1)->InterfaceId, (EleArr + ele - 1)->Id);
    fprintf(fStrEle, "%d %le %le %le %le %le %le\n", (EleArr + ele - 1)->G.Type,
            (EleArr + ele - 1)->G.Origin.X, (EleArr + ele - 1)->G.Origin.Y,
            (EleArr + ele - 1)->G.Origin.Z, (EleArr + ele - 1)->G.LX,
            (EleArr + ele - 1)->G.LZ, (EleArr + ele - 1)->G.dA);
    fprintf(fStrEle, "%le %le %le\n", (EleArr + ele - 1)->G.DC.XUnit.X,
            (EleArr + ele - 1)->G.DC.XUnit.Y, (EleArr + ele - 1)->G.DC.XUnit.Z);
    fprintf(fStrEle, "%le %le %le\n", (EleArr + ele - 1)->G.DC.YUnit.X,
            (EleArr + ele - 1)->G.DC.YUnit.Y, (EleArr + ele - 1)->G.DC.YUnit.Z);
    fprintf(fStrEle, "%le %le %le\n", (EleArr + ele - 1)->G.DC.ZUnit.X,
            (EleArr + ele - 1)->G.DC.ZUnit.Y, (EleArr + ele - 1)->G.DC.ZUnit.Z);
    fprintf(fStrEle, "%d %le\n", (EleArr + ele - 1)->E.Type,
            (EleArr + ele - 1)->E.Lambda);
    fprintf(fStrEle, "%d %le %le %le %le\n", (EleArr + ele - 1)->BC.NbOfBCs,
            (EleArr + ele - 1)->BC.CollPt.X, (EleArr + ele - 1)->BC.CollPt.Y,
            (EleArr + ele - 1)->BC.CollPt.Z, (EleArr + ele - 1)->BC.Value);
    fprintf(fStrEle, "%le %le\n", (EleArr + ele - 1)->Solution,
            (EleArr + ele - 1)->Assigned);
  }

  fprintf(fStrEle, "%d %d %d %d\n", NbPointsKnCh, NbLinesKnCh, NbAreasKnCh,
          NbVolumesKnCh);

  for (int pt = 1; pt <= NbPointsKnCh; ++pt) {
    fprintf(fStrEle, "%d %le\n", (PointKnChArr + pt - 1)->Nb,
            (PointKnChArr + pt - 1)->Assigned);
    fprintf(fStrEle, "%le %le %le\n", (PointKnChArr + pt - 1)->P.X,
            (PointKnChArr + pt - 1)->P.Y, (PointKnChArr + pt - 1)->P.Z);
  }

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    fprintf(fStrEle, "%d %le %le\n", (LineKnChArr + line - 1)->Nb,
            (LineKnChArr + line - 1)->Radius,
            (LineKnChArr + line - 1)->Assigned);
    fprintf(fStrEle, "%le %le %le\n", (LineKnChArr + line - 1)->Start.X,
            (LineKnChArr + line - 1)->Start.Y,
            (LineKnChArr + line - 1)->Start.Z);
    fprintf(fStrEle, "%le %le %le\n", (LineKnChArr + line - 1)->Stop.X,
            (LineKnChArr + line - 1)->Stop.Y, (LineKnChArr + line - 1)->Stop.Z);
  }

  for (int area = 1; area <= NbAreasKnCh; ++area) {
    fprintf(fStrEle, "%d %d %le\n", (AreaKnChArr + area - 1)->Nb,
            (AreaKnChArr + area - 1)->NbVertices,
            (AreaKnChArr + area - 1)->Assigned);
    for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices; ++vert) {
      fprintf(fStrEle, "%le %le %le\n",
              (AreaKnChArr + area - 1)->Vertex[vert].X,
              (AreaKnChArr + area - 1)->Vertex[vert].Y,
              (AreaKnChArr + area - 1)->Vertex[vert].Z);
    }
  }

  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    fprintf(fStrEle, "%d %d %le\n", (VolumeKnChArr + vol - 1)->Nb,
            (VolumeKnChArr + vol - 1)->NbVertices,
            (VolumeKnChArr + vol - 1)->Assigned);
    for (int vert = 1; vert <= (VolumeKnChArr + vol - 1)->NbVertices; ++vert) {
      fprintf(fStrEle, "%le %le %le\n",
              (VolumeKnChArr + vol - 1)->Vertex[vert].X,
              (VolumeKnChArr + vol - 1)->Vertex[vert].Y,
              (VolumeKnChArr + vol - 1)->Vertex[vert].Z);
    }
  }

  fclose(fStrEle);

  return 0;
}  // WriteElements ends

int ReadPrimitives(void) {
  int dbgFn = 0;

  char PrimitiveFile[256];

  strcpy(PrimitiveFile, ModelOutDir);
  strcat(PrimitiveFile, "/Primitives/StorePrims.out");

  FILE *fStrPrm;
  fStrPrm = fopen(PrimitiveFile, "r");
  if (fStrPrm == NULL) {
    neBEMMessage("ReadPrimitives - Could not open file to read primitives");
    return -1;
  }

  fscanf(fStrPrm, "%d %d\n", &NbVolumes, &VolMax);
  fscanf(fStrPrm, "%d\n", &NbPrimitives);
  fscanf(fStrPrm, "%d\n", &MaxNbVertices);

  // assign neBEMState and allocate memory
  neBEMState = 2;
  PrimType = ivector(1, NbPrimitives);
  NbVertices = ivector(1, NbPrimitives);
  XVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  YVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  ZVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  XNorm = dvector(1, NbPrimitives);
  YNorm = dvector(1, NbPrimitives);
  ZNorm = dvector(1, NbPrimitives);
  Radius = dvector(1, NbPrimitives);  // can lead to a little memory misuse
  VolRef1 = ivector(1, NbPrimitives);
  VolRef2 = ivector(1, NbPrimitives);
  NbSurfSegX = ivector(1, NbPrimitives);
  NbSurfSegZ = ivector(1, NbPrimitives);
  NbWireSeg = ivector(1, NbPrimitives);  // little memory misuse
  InterfaceType = ivector(1, NbPrimitives);
  Lambda = dvector(1, NbPrimitives);
  ApplPot = dvector(1, NbPrimitives);
  ApplCh = dvector(1, NbPrimitives);
  PeriodicTypeX = ivector(1, NbPrimitives);
  PeriodicTypeY = ivector(1, NbPrimitives);
  PeriodicTypeZ = ivector(1, NbPrimitives);
  PeriodicInX = ivector(1, NbPrimitives);
  PeriodicInY = ivector(1, NbPrimitives);
  PeriodicInZ = ivector(1, NbPrimitives);
  XPeriod = dvector(1, NbPrimitives);
  YPeriod = dvector(1, NbPrimitives);
  ZPeriod = dvector(1, NbPrimitives);
  MirrorDistXFromOrigin = dvector(1, NbPrimitives);
  MirrorDistYFromOrigin = dvector(1, NbPrimitives);
  MirrorDistZFromOrigin = dvector(1, NbPrimitives);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    fscanf(fStrPrm, "%d\n", &PrimType[prim]);
    fscanf(fStrPrm, "%d\n", &InterfaceType[prim]);
    fscanf(fStrPrm, "%d\n", &NbVertices[prim]);

    for (int vert = 0; vert < NbVertices[prim]; ++vert) {
      fscanf(fStrPrm, "%le %le %le\n", &XVertex[prim][vert],
             &YVertex[prim][vert], &ZVertex[prim][vert]);
    }  // vert loop

    fscanf(fStrPrm, "%le %le %le\n", &XNorm[prim], &YNorm[prim], &ZNorm[prim]);
    fscanf(fStrPrm, "%le\n", &Radius[prim]);

    fscanf(fStrPrm, "%le %le %le %le %le\n", &Epsilon1[prim], &Epsilon2[prim],
           &Lambda[prim], &ApplPot[prim], &ApplCh[prim]);

    fscanf(fStrPrm, "%d %d\n", &VolRef1[prim], &VolRef2[prim]);

    fscanf(fStrPrm, "%d %d %d\n", &PeriodicTypeX[prim], &PeriodicTypeY[prim],
           &PeriodicTypeZ[prim]);
    fscanf(fStrPrm, "%d %d %d\n", &PeriodicInX[prim], &PeriodicInY[prim],
           &PeriodicInZ[prim]);
    fscanf(fStrPrm, "%le %le %le\n", &XPeriod[prim], &YPeriod[prim],
           &ZPeriod[prim]);
    fscanf(fStrPrm, "%le %le %le\n", &MirrorDistXFromOrigin[prim],
           &MirrorDistYFromOrigin[prim], &MirrorDistZFromOrigin[prim]);
  }  // prim loop

  volRef = ivector(0, VolMax);
  volShape = ivector(0, VolMax);
  volMaterial = ivector(0, VolMax);
  volEpsilon = dvector(0, VolMax);
  volPotential = dvector(0, VolMax);
  volCharge = dvector(0, VolMax);
  volBoundaryType = ivector(0, VolMax);
  for (int volref = 0; volref <= VolMax; ++volref) {
    neBEMVolumeDescription(volref, &volShape[volref], &volMaterial[volref],
                           &volEpsilon[volref], &volPotential[volref],
                           &volCharge[volref], &volBoundaryType[volref]);
    if (dbgFn) {
      printf("volref: %d\n", volref);
      printf("shape: %d,  material: %d\n", volShape[volref],
             volMaterial[volref]);
      printf("eps: %lg,  pot: %lg\n", volEpsilon[volref], volPotential[volref]);
      printf("q: %lg,  type: %d\n", volCharge[volref], volBoundaryType[volref]);
    }
  }  // volume loop

  fclose(fStrPrm);

  return 0;
}  // ReadPrimitives ends

int ReadElements(void) {
  char ElementFile[256];

  strcpy(ElementFile, MeshOutDir);
  strcat(ElementFile, "/Elements/StoreElems.out");

  FILE *fStrEle = fopen(ElementFile, "r");
  if (fStrEle == NULL) {
    neBEMMessage("ReadElements - Could not open file to read elements");
    return -1;
  }

  fscanf(fStrEle, "%d %d\n", &NbSurfs, &NbWires);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    fscanf(fStrEle, "%d %d\n", &NbSurfSegX[prim], &NbSurfSegZ[prim]);
    fscanf(fStrEle, "%d\n", &NbWireSeg[prim]);
  }

  fscanf(fStrEle, "%d\n", &NbElements);

  if (neBEMState == 3) {
    if (EleArr) {
      Element *tmp = (Element *)realloc(EleArr, NbElements * sizeof(Element));
      if (tmp != NULL) {
        EleArr = tmp;
        EleCntr = 0;
      } else {
        free(EleArr);
        printf("neBEMDiscretize: Re-allocating EleArr failed.\n");
        return (1);
      }
      printf("neBEMDiscretize: Re-allocated EleArr.\n");
    }  // if EleArr => re-allocation
    else {
      EleArr = (Element *)malloc(NbElements * sizeof(Element));
      if (EleArr == NULL) {
        neBEMMessage("neBEMDiscretize - EleArr malloc");
        return -1;
      }
    }  // else EleArr => fresh allocation
  } else {
    neBEMMessage("neBEMDiscretize - EleArr malloc; neBEMState mismatch!");
    return -1;
  }  // else neBEMState == 3

  for (int ele = 1; ele <= NbElements; ++ele) {
    fscanf(fStrEle, "%hd %d %d %d %d\n", &(EleArr + ele - 1)->DeviceNb,
           &(EleArr + ele - 1)->ComponentNb, &(EleArr + ele - 1)->PrimitiveNb,
           &(EleArr + ele - 1)->InterfaceId, &(EleArr + ele - 1)->Id);
    fscanf(fStrEle, "%hd %le %le %le %le %le %le\n",
           &(EleArr + ele - 1)->G.Type, &(EleArr + ele - 1)->G.Origin.X,
           &(EleArr + ele - 1)->G.Origin.Y, &(EleArr + ele - 1)->G.Origin.Z,
           &(EleArr + ele - 1)->G.LX, &(EleArr + ele - 1)->G.LZ,
           &(EleArr + ele - 1)->G.dA);
    fscanf(fStrEle, "%le %le %le\n", &(EleArr + ele - 1)->G.DC.XUnit.X,
           &(EleArr + ele - 1)->G.DC.XUnit.Y,
           &(EleArr + ele - 1)->G.DC.XUnit.Z);
    fscanf(fStrEle, "%le %le %le\n", &(EleArr + ele - 1)->G.DC.YUnit.X,
           &(EleArr + ele - 1)->G.DC.YUnit.Y,
           &(EleArr + ele - 1)->G.DC.YUnit.Z);
    fscanf(fStrEle, "%le %le %le\n", &(EleArr + ele - 1)->G.DC.ZUnit.X,
           &(EleArr + ele - 1)->G.DC.ZUnit.Y,
           &(EleArr + ele - 1)->G.DC.ZUnit.Z);
    fscanf(fStrEle, "%hd %le\n", &(EleArr + ele - 1)->E.Type,
           &(EleArr + ele - 1)->E.Lambda);
    fscanf(fStrEle, "%hd %le %le %le %le\n", &(EleArr + ele - 1)->BC.NbOfBCs,
           &(EleArr + ele - 1)->BC.CollPt.X, &(EleArr + ele - 1)->BC.CollPt.Y,
           &(EleArr + ele - 1)->BC.CollPt.Z, &(EleArr + ele - 1)->BC.Value);
    fscanf(fStrEle, "%le %le\n", &(EleArr + ele - 1)->Solution,
           &(EleArr + ele - 1)->Assigned);
  }

  fscanf(fStrEle, "%d %d %d %d\n", &NbPointsKnCh, &NbLinesKnCh, &NbAreasKnCh,
         &NbVolumesKnCh);

  for (int pt = 1; pt <= NbPointsKnCh; ++pt) {
    fscanf(fStrEle, "%d %le\n", &(PointKnChArr + pt - 1)->Nb,
           &(PointKnChArr + pt - 1)->Assigned);
    fscanf(fStrEle, "%le %le %le\n", &(PointKnChArr + pt - 1)->P.X,
           &(PointKnChArr + pt - 1)->P.Y, &(PointKnChArr + pt - 1)->P.Z);
  }

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    fscanf(fStrEle, "%d %le %le\n", &(LineKnChArr + line - 1)->Nb,
           &(LineKnChArr + line - 1)->Radius,
           &(LineKnChArr + line - 1)->Assigned);
    fscanf(fStrEle, "%le %le %le\n", &(LineKnChArr + line - 1)->Start.X,
           &(LineKnChArr + line - 1)->Start.Y,
           &(LineKnChArr + line - 1)->Start.Z);
    fscanf(fStrEle, "%le %le %le\n", &(LineKnChArr + line - 1)->Stop.X,
           &(LineKnChArr + line - 1)->Stop.Y,
           &(LineKnChArr + line - 1)->Stop.Z);
  }

  for (int area = 1; area <= NbAreasKnCh; ++area) {
    fscanf(fStrEle, "%d %d %le\n", &(AreaKnChArr + area - 1)->Nb,
           &(AreaKnChArr + area - 1)->NbVertices,
           &(AreaKnChArr + area - 1)->Assigned);
    for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices; ++vert) {
      fscanf(fStrEle, "%le %le %le\n",
             &(AreaKnChArr + area - 1)->Vertex[vert].X,
             &(AreaKnChArr + area - 1)->Vertex[vert].Y,
             &(AreaKnChArr + area - 1)->Vertex[vert].Z);
    }
  }

  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    fscanf(fStrEle, "%d %d %le\n", &(VolumeKnChArr + vol - 1)->Nb,
           &(VolumeKnChArr + vol - 1)->NbVertices,
           &(VolumeKnChArr + vol - 1)->Assigned);
    for (int vert = 1; vert <= (VolumeKnChArr + vol - 1)->NbVertices; ++vert) {
      fscanf(fStrEle, "%le %le %le\n",
             &(VolumeKnChArr + vol - 1)->Vertex[vert].X,
             &(VolumeKnChArr + vol - 1)->Vertex[vert].Y,
             &(VolumeKnChArr + vol - 1)->Vertex[vert].Z);
    }
  }

  fclose(fStrEle);

  return 0;
}  // ReadElements ends

// returns the time elapsed between start and stop
// Interface.h
double neBEMTimeElapsed(clock_t startClock, clock_t stopClock) {
  double elapsedTime;

  elapsedTime = ((double)(stopClock - startClock)) / CLOCKS_PER_SEC;
  printf("neBEMV%s TimeElapsed ===> %lg seconds ", neBEMVersion, elapsedTime);

  return (elapsedTime);
}

// returns number of lines in file fname
// Interface.h
int neBEMGetNbOfLines(const char fname[]) {
  unsigned int number_of_lines = 0;

  FILE *infile = fopen(fname, "r");

  int ch;
  while (EOF != (ch = getc(infile)))
    if ('\n' == ch) ++number_of_lines;

  return number_of_lines;
}

// check whether a 3D point is within a 3D polygon
// Downloads/ComNum/Geometry/Determining if a point lies on the interior of a
// polygon.htm
double neBEMChkInPoly(int n, Point3D *p, Point3D q) {
  double anglesum = 0.0;
  // double theta = 0.0;
  Point3D p1, p2;

  // printf("In neBEMChkInPoly ... \n");
  // printf("n: %d\n", n);

  for (int i = 0; i < n; i++) {
    // printf("i: %d\n", i);
    p1.X = p[i].X - q.X;
    p1.Y = p[i].Y - q.Y;
    p1.Z = p[i].Z - q.Z;
    // printf("p[%d]: %lg %lg %lg\n", i, p[i].X, p[i].Y, p[i].Z);
    // printf("q: %lg %lg %lg\n", q.X, q.Y, q.Z);

    if (i < n - 1) {
      p2.X = p[i + 1].X - q.X;
      p2.Y = p[i + 1].Y - q.Y;
      p2.Z = p[i + 1].Z - q.Z;
    } else {
      p2.X = p[0].X - q.X;
      p2.Y = p[0].Y - q.Y;
      p2.Z = p[0].Z - q.Z;
    }

    double m1 = MODULUS(p1);
    double m2 = MODULUS(p2);
    // printf("m1: %lg, m2: %lg, m1*m2: %lg", m1, m2, m1*m2);

    if (m1 * m2 <= 1.0e-12) {
      // vetors of 1 micron - we may need to reduce further
      return (neBEMtwopi); /* We are on a node, consider this inside */
    }
    double costheta = (p1.X * p2.X + p1.Y * p2.Y + p1.Z * p2.Z) / (m1 * m2);

    /*
    double oldtheta = theta;
    theta = acos(costheta);
    // printf("n: %d, i: %d, theta: %lg\n", n, i, neBEMrtod*theta);
    if (Sign(theta) != Sign(oldtheta)) {
      // polygon either non-covex, or the point is outside the polygon 
      return(0.0);  // absurd value implying outside polygon
    }
    */
    anglesum += acos(costheta);
    // printf("n: %d, i: %d, anglesum: %lg %lg\n", n, i, anglesum, neBEMrtod*anglesum);
  }

  return (anglesum);
}  // neBEMChkInPoly

#ifdef __cplusplus
} // namespace
#endif
