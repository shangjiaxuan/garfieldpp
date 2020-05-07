#include <stdio.h>

#include "neBEMInterface.h"
#include "NR.h"
#include "Vector.h"
#include "neBEM.h"

#include "Garfield/ComponentNeBem3d.hh"

/// Assign default values to some of the important global variables.
int neBEMSetDefaults(void) {

  neBEMState = 0;

  NbVolumes = 2;
  VolMax = 2;

  NbPrimitives = 1;
  MaxNbVertices = 4;

  NbSurfs = 1;
  NbWires = 0;

  MinNbElementsOnLength = 1;
  MaxNbElementsOnLength = 100;
  ElementLengthRqstd = 100.0e-6;

  NewModel = 1;
  NewMesh = 1;
  NewBC = 1;
  NewPP = 1;
  ModelCntr = 1;
  MeshCntr = 1;
  BCCntr = 1;
  PPCntr = 1;

  DebugLevel = 0;

  LengthScale = 1.0;

  TimeStep = 1;

  strcpy(DeviceOutDir, "Outputs");

  OptDeviceFile = 0;
  strcpy(DeviceInputFile, "");

  OptGnuplot = 0;
  OptGnuplotPrimitives = 0;
  OptGnuplotElements = 0;
  OptPrimitiveFiles = 0;
  OptElementFiles = 0;

  OptPrintPrimaryDetails = 0;
  OptPrintVolumeDetails = 0;
  OptPrintVertexAndNormal = 0;

  OptGnuplot = 1;
  OptGnuplotPrimitives = 1;
  OptGnuplotElements = 1;
  OptPrimitiveFiles = 1;
  OptElementFiles = 1;

  // Same directory can be used over and over again - BEWARE!
  OptReuseDir = 1;  

  // Matrix inversion procedure => 0: LU, 1: SVD, 2:GSL
  OptInvMatProc = 0;  

  OptValidateSolution = 0;
  OptForceValidation = 0;
  OptStorePrimitives = 1;
  OptStoreElements = 1;
  OptStoreInflMatrix = 0;
  OptStoreInvMatrix = 1;
  OptFormattedFile = 1;
  OptUnformattedFile = 0;
  OptRepeatLHMatrix = 0;

  OptSystemChargeZero = 1;

  return 0;
}  // neBEMSetDefaults ends

int ReadInitFile(char filename[]) {
  FILE* finit = fopen(filename, "r");
  if (finit == NULL) {
    neBEMMessage("ReadInitFile - fail to open init file");
    return -1;
  }

  fscanf(finit, "MinNbElementsOnLength: %d\n", &MinNbElementsOnLength);
  fscanf(finit, "MaxNbElementsOnLength: %d\n", &MaxNbElementsOnLength);
  fscanf(finit, "ElementLengthRqstd: %le\n", &ElementLengthRqstd);
  fscanf(finit, "LengthScale: %le\n", &LengthScale);

  fscanf(finit, "DebugLevel: %d\n", &DebugLevel);

  fscanf(finit, "NewModel: %d\n", &NewModel);
  fscanf(finit, "NewMesh: %d\n", &NewMesh);
  fscanf(finit, "NewBC: %d\n", &NewBC);
  fscanf(finit, "NewPP: %d\n", &NewPP);
  fscanf(finit, "ModelCntr: %d\n", &ModelCntr);
  fscanf(finit, "MeshCntr: %d\n", &MeshCntr);
  fscanf(finit, "BCCntr: %d\n", &BCCntr);
  fscanf(finit, "PPCntr: %d\n", &PPCntr);

  fscanf(finit, "DeviceOutDir: %s\n", DeviceOutDir);

  fscanf(finit, "OptDeviceFile: %d\n", &OptDeviceFile);
  fscanf(finit, "DeviceInputFile: %s\n", DeviceInputFile);

  fscanf(finit, "OptPrintPrimaryDetails: %d\n", &OptPrintPrimaryDetails);
  fscanf(finit, "OptPrintVolumeDetails: %d\n", &OptPrintVolumeDetails);
  fscanf(finit, "OptPrintVertexAndNormal: %d\n", &OptPrintVertexAndNormal);

  fscanf(finit, "OptGnuplot: %d\n", &OptGnuplot);
  fscanf(finit, "OptGnuplotPrimitives: %d\n", &OptGnuplotPrimitives);
  fscanf(finit, "OptGnuplotElements: %d\n", &OptGnuplotElements);

  fscanf(finit, "OptPrimitiveFiles: %d\n", &OptPrimitiveFiles);
  fscanf(finit, "OptElementFiles: %d\n", &OptElementFiles);

  fscanf(finit, "OptReuseDir: %d\n", &OptReuseDir);

  fscanf(finit, "OptInvMatProc: %d\n", &OptInvMatProc);

  fscanf(finit, "OptValidateSolution: %d\n", &OptValidateSolution);
  fscanf(finit, "OptForceValidation: %d\n", &OptForceValidation);
  fscanf(finit, "OptStorePrimitives: %d\n", &OptStorePrimitives);
  fscanf(finit, "OptStoreElements: %d\n", &OptStoreElements);
  fscanf(finit, "OptStoreInflMatrix: %d\n", &OptStoreInflMatrix);
  fscanf(finit, "OptStoreInvMatrix: %d\n", &OptStoreInvMatrix);
  fscanf(finit, "OptFormattedFile: %d\n", &OptFormattedFile);
  fscanf(finit, "OptUnformattedFile: %d\n", &OptUnformattedFile);
  fscanf(finit, "OptRepeatLHMatrix: %d\n", &OptRepeatLHMatrix);

  fscanf(finit, "OptSystemChargeZero: %d\n", &OptSystemChargeZero);

  fscanf(finit, "PrimAfter: %d\n", &PrimAfter);

  fclose(finit);

  fprintf(stdout, "MinNbElementsOnLength: %d\n", MinNbElementsOnLength);
  fprintf(stdout, "MaxNbElementsOnLength: %d\n", MaxNbElementsOnLength);
  fprintf(stdout, "ElementLengthRqstd: %le\n", ElementLengthRqstd);
  fprintf(stdout, "LengthScale: %le\n", LengthScale);

  fprintf(stdout, "NewModel: %d\n", NewModel);
  fprintf(stdout, "NewMesh: %d\n", NewMesh);
  fprintf(stdout, "NewBC: %d\n", NewBC);
  fprintf(stdout, "NewPP: %d\n", NewPP);
  fprintf(stdout, "ModelCntr: %d\n", ModelCntr);
  fprintf(stdout, "MeshCntr: %d\n", MeshCntr);
  fprintf(stdout, "BCCntr: %d\n", BCCntr);
  fprintf(stdout, "PPCntr: %d\n", PPCntr);

  fprintf(stdout, "DeviceOutDir: %s\n", DeviceOutDir);

  fprintf(stdout, "OptDeviceFile: %d\n", OptDeviceFile);
  fprintf(stdout, "DeviceInputFile: %s\n", DeviceInputFile);

  fprintf(stdout, "OptPrintPrimaryDetails: %d\n", OptPrintPrimaryDetails);
  fprintf(stdout, "OptPrintVolumeDetails: %d\n", OptPrintVolumeDetails);
  fprintf(stdout, "OptPrintVertexAndNormal: %d\n", OptPrintVertexAndNormal);

  fprintf(stdout, "OptGnuplot: %d\n", OptGnuplot);
  fprintf(stdout, "OptGnuplotPrimitives: %d\n", OptGnuplotPrimitives);
  fprintf(stdout, "OptGnuplotElements: %d\n", OptGnuplotElements);

  fprintf(stdout, "OptPrimitiveFiles: %d\n", OptPrimitiveFiles);
  fprintf(stdout, "OptElementFiles: %d\n", OptElementFiles);

  fprintf(stdout, "OptReuseDir: %d\n", OptReuseDir);

  fprintf(stdout, "OptValidateSolution: %d\n", OptValidateSolution);
  fprintf(stdout, "OptStorePrimitives: %d\n", OptStorePrimitives);
  fprintf(stdout, "OptStoreElements: %d\n", OptStoreElements);
  fprintf(stdout, "OptStoreInflMatrix: %d\n", OptStoreInflMatrix);
  fprintf(stdout, "OptStoreInvMatrix: %d\n", OptStoreInvMatrix);
  fprintf(stdout, "OptFormattedFile: %d\n", OptFormattedFile);
  fprintf(stdout, "OptUnformattedFile: %d\n", OptUnformattedFile);
  fprintf(stdout, "OptRepeatLHMatrix: %d\n", OptRepeatLHMatrix);

  fprintf(stdout, "OptSystemChargeZero: %d\n", OptSystemChargeZero);

  fprintf(stdout, "PrimAfter: %d\n", PrimAfter);

  return 0;
}

/// Do-nothing function (no file inputs).
int neBEMGetInputsFromFiles(void) { return 0; } 

/// Return the number of primitives.
int neBEMGetNbPrimitives() {
  if (!Garfield::gComponentNeBem3d) return 0;
  return Garfield::gComponentNeBem3d->GetNumberOfPrimitives();
}

/// Return one primitive at a time.
int neBEMGetPrimitive(int prim, int* nvertex, double xvert[], double yvert[],
                      double zvert[], double* xnorm, double* ynorm,
                      double* znorm, int* volref1, int* volref2) {

  if (!Garfield::gComponentNeBem3d) return -1;
  if (prim < 1) return -1;
  double a = 0., b = 0., c = 0.;
  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  int vol1 = 0, vol2 = 0;
  if (!Garfield::gComponentNeBem3d->GetPrimitive(prim - 1, a, b, c, 
                                                 xv, yv, zv, vol1, vol2)) {
    return -1;
  }
  const size_t nv = xv.size();
  *nvertex = nv;
  for (size_t i = 0; i < nv; ++i) {
    // Convert from cm to m.
    xvert[i] = 0.01 * xv[i];
    yvert[i] = 0.01 * yv[i];
    zvert[i] = 0.01 * zv[i];
  }
  *xnorm = a;
  *ynorm = b;
  *znorm = c;
  if (nv == 2) *xnorm *= 0.01;
  *volref1 = vol1;
  *volref2 = vol2;
  return 0;
}

/// Return information about periodicities.
/// ix: 0 = no periodicity in x,
///     1 = simple periodicity of length sx,
///     2 = mirror periodicity of length sx,
///     3 = axial periodicity around x with sector angle sx,
///     4 = rotational symmetry around x
/// jx: number of periodic copies internally in neBEM
int neBEMGetPeriodicities(int /*prim*/, int* ix, int* jx, double* sx, int* iy,
                          int* jy, double* sy, int* iz, int* jz, double* sz) {
  // Not implemented yet. Assume no periodicity for the time being.
  *ix = 0;
  *iy = 0;
  *iz = 0;
  *jx = 0;
  *jy = 0;
  *jz = 0;

  // Convert from cm to m.
  *sx *= 0.01;
  *sy *= 0.01;
  *sz *= 0.01;

  return 0;
}

/// Return information about mirror periodicity.
/// ix: 0 = no mirror in x,
///     1 = producing charge density of opposite sign (infinite cond plane)
///     2 = producing charge density of same sign
/// jx: not used at present
/// sx: x distance of the mirror from origin
int neBEMGetMirror(int /*prim*/, int* ix, int* jx, double* sx, int* iy, int* jy,
                   double* sy, int* iz, int* jz, double* sz) {

  // Not implemented yet. Assume no mirror periodicity for the time being.
  *ix = *iy = *iz = 0;
  *jx = *jy = *jz = 0;
  // Mirror assumed to be passing through the origin.
  *sx = *sy = *sz = 0.;  

  // Convert from cm to m.
  *sx *= 0.01;
  *sy *= 0.01;
  *sz *= 0.01;

  return 0;
}

/// Return infinite (bounding) conductors if any.
/// ixmin=0: lower x-plane does not exist
/// ixmin=1: lower x-plane does exist
/// cxmin: coordinate of lower x-plane
/// vxmin: potential of lower x-plane
/// Similar for ixmax, iymin, iymax, izmin, izmax
int neBEMGetBoundingPlanes(int* ixmin, double* cxmin, double* vxmin, int* ixmax,
                           double* cxmax, double* vxmax, int* iymin,
                           double* cymin, double* vymin, int* iymax,
                           double* cymax, double* vymax, int* izmin,
                           double* czmin, double* vzmin, int* izmax,
                           double* czmax, double* vzmax) {
  // Not implemented yet. Assume no planes for the time being.
  *ixmin = *ixmax = 0;
  *iymin = *iymax = 0;
  *izmin = *izmax = 0;
  *vxmin = *vxmax = 0.;
  *vymin = *vymax = 0.;
  *vzmin = *vzmax = 0.;
  // Convert from cm to m.
  *cxmin *= 0.01;
  *cxmax *= 0.01;
  *cymin *= 0.01;
  *cymax *= 0.01;
  *czmin *= 0.01;
  *czmax *= 0.01;

  return 0;
}

/// Return information about a volume.
int neBEMVolumeDescription(int vol, int* shape, int* material,
                           double* epsilon, double* potential, double* charge,
                           int* boundarytype) {
  
  if (!Garfield::gComponentNeBem3d) return -1;
  if (!Garfield::gComponentNeBem3d->GetVolume(vol, *shape, *material, 
                                              *epsilon, *potential, *charge,
                                              *boundarytype)) {
    return false;
  }
  return true; 
}

/// Return the volume in which a point is located.
int neBEMVolumePoint(double x, double y, double z) {
  
  if (!Garfield::gComponentNeBem3d) return -1;
  return Garfield::gComponentNeBem3d->GetVolume(100. * x, 100. * y, 100. * z);
}

/// Return the primitives for a volume.
/// TODO! Do we need this?
void neBEMVolumePrimitives(int /*vol*/, int* /*nprim*/, int /*primlist*/[]) {

}

