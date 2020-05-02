// Gets number of primitives (polygonal surfaces (2D) and wires (1D)) and
// their details
// Assumptions:
//	MaxNbVertices: 4 (maximum nb of vertices a polygon can have)

#include <stdio.h>
#include <sys/stat.h>

#include "Interface.h"
#include "NR.h"
#include "Vector.h"
#include "neBEM.h"

// -------------------------------------------------------------------

void bemini_(int* ifail) {
  // Garfield calls this procedure when initialising
  *ifail = neBEMInitialize();
}

// -------------------------------------------------------------------

// Assign default values to some of the important global variables.
// These variables can also be read in from a configuration file.
int neBEMSetDefaults(void) {
  int ReadInitFile(char filename[]);

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

  char* homedir;
  homedir = getenv("HOME");
  if (homedir == NULL) neBEMMessage("getenv homedir");
  printf("user homedir is %s\n", homedir);

  strcpy(DeviceOutDir, "Outputs");
  // strcpy(DeviceOutDir, homedir);
  // strcat(DeviceOutDir, "/neBEMWD/Outputs/Garfield/V");
  // strcat(DeviceOutDir, neBEMVersion);
  // strcat(DeviceOutDir, "/");

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

  OptReuseDir = 1;  // Same directory can be used over and over again - BEWARE!

  OptInvMatProc = 0;  // matrix inversion procedure => 0: LU, 1: SVD, 2:GSL

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

  // The above parameters can be overridden by an `init' file.
  // No initialization file usually for Garfield. However, it may be allowed
  // if felt necessary.
  /*
  char initFileName[256];
  strcpy(initFileName, homedir);
  strcat(initFileName, "/neBEMWD/InitFiles/");
  strcat(initFileName, "Garfield2neBEM.init");

  struct stat st;
  if(stat(initFileName, &st) == 0)
          {
          printf("Found an init file ... updating parameters.\n");
          int fstatus = ReadInitFile(initFileName);
          if(fstatus != 0)neBEMMessage("neBEMSetDefaults - ReadInitFile");
          }
  else
          {
          printf("No init file ... continuing with the parameter values.\n");
          }
  */

  return (0);
}  // neBEMSetDefaults ends

// -------------------------------------------------------------------
void bempar_(int* bemmin, int* bemmax, double* bemtarget, int* newflag,
             int* storeinverted, int* debuglevel, int* inversionmethod) {
  // Set the user options
  MinNbElementsOnLength = *bemmin;
  MaxNbElementsOnLength = *bemmax;
  ElementLengthRqstd = (*bemtarget) * 0.01;

  // New model / reuse existing model flag
  if (*newflag) {
    NewModel = 1;
    NewMesh = 1;
    NewBC = 1;
    NewPP = 1;
  } else {
    NewModel = 0;
    NewMesh = 0;
    NewBC = 1;  // Changed by Supratik
    NewPP = 1;
  }

  // Pass debug level
  DebugLevel = *debuglevel;

  // Store inverted matrix or not
  if (*storeinverted) {
    OptStoreInvMatrix = 1;
    OptFormattedFile = 1;
  } else {
    OptStoreInvMatrix = 0;
    OptFormattedFile = 0;
  }

  // Matrix inversion method
  OptInvMatProc = *inversionmethod;

  // Debugging
  printf("bempar: Discretisation target size: %g\n", ElementLengthRqstd);
  printf("        Element count range: %d to %d\n", MinNbElementsOnLength,
         MaxNbElementsOnLength);
  printf("        Model = %d, Mesh = %d, BC = %d, PP = %d\n", NewModel, NewMesh,
         NewBC, NewPP);
  printf("        Store inverted = %d, formatted file = %d\n",
         OptStoreInvMatrix, OptFormattedFile);
  printf("        Inversion method = %d (0 = LU, 1 = SVD)\n", OptInvMatProc);
}

// -------------------------------------------------------------------
int ReadInitFile(char filename[]) {
  FILE* finit;
  finit = fopen(filename, "r");
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

// -------------------------------------------------------------------

// Dummy function doing nothing since, while using Garfield, we do not
// have file inputs.
int neBEMGetInputsFromFiles(void) { return (0); }

// -------------------------------------------------------------------

void bemgeo_(int* ifail) {
  // Garfield calls this procedure to ask neBEM to start to read
  *ifail = neBEMReadGeometry();
}

// -------------------------------------------------------------------

// Garfield procedure for counting the primitives
// TODO!
// extern int bemnpr_();

int neBEMGetNbPrimitives() {
  // Number of panels
  // TODO!
  // return bemnpr_();
  return 0;
}

// -------------------------------------------------------------------

// Garfield procedure that returns one primitive at the time
// TODO!
/*
extern void bempri_(int* elem, int* nvertex,
                   double* xvert, double* yvert, double* zvert,
                   double* xnorm, double* ynorm, double* znorm,
                   int* volref1, int* volref2, int* ifail);
*/
int neBEMGetPrimitive(int prim, int* nvertex, double xvert[], double yvert[],
                      double zvert[], double* xnorm, double* ynorm,
                      double* znorm, int* volref1, int* volref2) {
  int ifail, ee = prim, ivertex;
  // why +1? - changed by supratik here and in neBEMInterface call
  // TODO!
  // bempri_(&ee, nvertex, xvert, yvert, zvert, xnorm, ynorm, znorm,
  //         volref1, volref2, &ifail);
  // Convert from cm to m.
  for (ivertex = 0; ivertex < *nvertex; ivertex++) {
    xvert[ivertex] *= 0.01;
    yvert[ivertex] *= 0.01;
    zvert[ivertex] *= 0.01;
  }
  if (*nvertex == 2) {
    *xnorm *= 0.01;
  }

  // Adjust the volume numbering
  (*volref1)--;
  (*volref2)--;

  // Return error status.
  return ifail;
}

// -------------------------------------------------------------------

// Garfield procedure returning the periodicities
// ix: 0 = no periodicity in x,
//     1 = simple periodicity of length sx,
//     2 = mirror periodicity of length sx,
//     3 = axial periodicity around x with sector angle sx,
//     4 = rotational symmetry around x
// jx: number of periodic copies internally in neBEM
// TODO!
/*
extern int bemper_(int* ivol,
                   int* ix, int* jx, double* sx,
                   int* iy, int* jy, double* sy,
                   int* iz, int* jz, double* sz);
*/

int neBEMGetPeriodicities(int prim, int* ix, int* jx, double* sx, int* iy,
                          int* jy, double* sy, int* iz, int* jz, double* sz) {
  // Obtain the periodicities
  // TODO!
  // int ifail = bemper_(&prim, ix, jx, sx, iy, jy, sy, iz, jz, sz);
  int ifail = 0;

  // Convert from cm to m.
  *sx *= 0.01;
  *sy *= 0.01;
  *sz *= 0.01;

  // Return error status.
  return ifail;
}

// -------------------------------------------------------------------

// Garfield procedure returning the mirror information
// ix: 0 = no mirror in x,
//     1 = producing charge density of oppostie sign (infinite cond plane)
//     2 = producing charge density of same sign
// jx: not used at present
// sx: x distance of the mirror from origin
// TODO!
/*
extern int bemper_(int* ivol,
                   int* ix, int* jx, double* sx,
                   int* iy, int* jy, double* sy,
                   int* iz, int* jz, double* sz);
*/

int neBEMGetMirror(int prim, int* ix, int* jx, double* sx, int* iy, int* jy,
                   double* sy, int* iz, int* jz, double* sz) {
  // Obtain mirror information
  // int ifail = bemper_(&prim, ix, jx, sx, iy, jy, sy, iz, jz, sz);
  int ifail = 0;

  *jx = *jy = *jz = 0;    // not used at present
  *sx = *sy = *sz = 0.0;  // mirror assumed to be passing through the origin

  // Only value 2 implies mirror periodicity
  // At present MirrorTypeX etc is always 2. When Mirrors of other kinds will be
  // allowed, this will need modification.
  if ((*ix != 2) && (*iy != 2) && (*iz != 2)) {
    *ix = *iy = *iz = 0;
  }

  // Only one reflection is allowed at present
  if (*ix == 2) {
    *iy = *iz = 0;
  } else if (*iy == 2)
    *iz = 0;

  // Convert from cm to m.
  *sx *= 0.01;
  *sy *= 0.01;
  *sz *= 0.01;

  // Return error status.
  return ifail;
}

// -------------------------------------------------------------------

// Garfield procedure returning the infinite (bounding) conductors
// ixmin=0: lower x-plane does not exist
// ixmin=1: lower x-plane does exist
// cxmin: coordinate of lower x-plane
// vxmin: potential of lower x-plane
// Similar for ixmax, iymin, iymax, izmin, izmax
/*
extern int bempla_(int* ixmin, double* cxmin, double* vxmin,
                   int* ixmax, double* cxmax, double* vxmax,
                   int* iymin, double* cymin, double* vymin,
                   int* iymax, double* cymax, double* vymax,
                   int* izmin, double* czmin, double* vzmin,
                   int* izmax, double* czmax, double* vzmax);
*/

int neBEMGetBoundingPlanes(int* ixmin, double* cxmin, double* vxmin, int* ixmax,
                           double* cxmax, double* vxmax, int* iymin,
                           double* cymin, double* vymin, int* iymax,
                           double* cymax, double* vymax, int* izmin,
                           double* czmin, double* vzmin, int* izmax,
                           double* czmax, double* vzmax) {
  // Obtain the periodicities
  // TODO!
  // int ifail = bempla_(ixmin, cxmin, vxmin, ixmax, cxmax, vxmax,
  //                     iymin, cymin, vymin, iymax, cymax, vymax,
  //                     izmin, czmin, vzmin, izmax, czmax, vzmax);
  int ifail = 0;

  // Convert from cm to m.
  *cxmin *= 0.01;
  *cxmax *= 0.01;
  *cymin *= 0.01;
  *cymax *= 0.01;
  *czmin *= 0.01;
  *czmax *= 0.01;

  // Supratik: tore the values somewhere

  // Return error status.
  return ifail;
}

// -------------------------------------------------------------------

// Garfield procedure that works out the target size for elements
// TODO!
// extern void bemszd_(int* prim, double* size, int* ifail);

void bemdis_(int* ifail) {
  int** elementNbs;
  int ifail1, nb1, nb2;
  double size, l1, l2;

  elementNbs = imatrix(1, NbPrimitives, 1, 2);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    // Try and identify the element
    // TODO!
    // bemszd_(&prim, &size, &ifail1);
    if (size <= 0 || ifail1 != 0) {
      printf(
          "bemdis: Primitive %d discretisation size set to target element size "
          "(%g m).\n",
          prim, ElementLengthRqstd);
      size = ElementLengthRqstd;
    } else {
      size *= 0.01;
    }

    // Work out the element dimensions (from AnalyzeSurface)
    l1 = (XVertex[prim][0] - XVertex[prim][1]) *
             (XVertex[prim][0] - XVertex[prim][1]) +
         (YVertex[prim][0] - YVertex[prim][1]) *
             (YVertex[prim][0] - YVertex[prim][1]) +
         (ZVertex[prim][0] - ZVertex[prim][1]) *
             (ZVertex[prim][0] - ZVertex[prim][1]);
    nb1 = (int)(sqrt(fabs(l1)) / size);
    // Truncate to desired range
    if (nb1 < MinNbElementsOnLength) {
      // printf("bemdis: Primitive %d discretisation 1 increased from %d to
      // minimum allowed (%d).\n", prim, nb1, MinNbElementsOnLength);
      nb1 = MinNbElementsOnLength;
    } else if (nb1 > MaxNbElementsOnLength) {
      printf(
          "bemdis: Primitive %d discretisation 1 decreased from %d to maximum "
          "allowed (%d).\n",
          prim, nb1, MaxNbElementsOnLength);
      nb1 = MaxNbElementsOnLength;
    }
    if (NbVertices[prim] > 2) {
      l2 = (XVertex[prim][2] - XVertex[prim][1]) *
               (XVertex[prim][2] - XVertex[prim][1]) +
           (YVertex[prim][2] - YVertex[prim][1]) *
               (YVertex[prim][2] - YVertex[prim][1]) +
           (ZVertex[prim][2] - ZVertex[prim][1]) *
               (ZVertex[prim][2] - ZVertex[prim][1]);
      nb2 = (int)(sqrt(fabs(l2)) / size);
      if (nb2 < MinNbElementsOnLength) {
        // printf("bemdis: Primitive %d discretisation 2 increased from %d to
        // minimum allowed (%d).\n", prim, nb2, MinNbElementsOnLength);
        nb2 = MinNbElementsOnLength;
      } else if (nb2 > MaxNbElementsOnLength) {
        printf(
            "bemdis: Primitive %d discretisation 2 decreased from %d to "
            "maximum allowed (%d).\n",
            prim, nb2, MaxNbElementsOnLength);
        nb2 = MaxNbElementsOnLength;
      }
    } else {
      l2 = 0;
      nb2 = 0;
    }
    // printf("bemdis: Primitive %d discretisation size = %g, primitive size l1
    // = %g, l2 = %g, n1 = %d, n2 = %d.\n", prim, size, sqrt(l1), sqrt(l2), nb1,
    // nb2);

    // Assign
    elementNbs[prim][1] = nb1;
    elementNbs[prim][2] = nb2;
  }
  *ifail = neBEMDiscretize(elementNbs);
}

// -------------------------------------------------------------------

void bembnd_(int* ifail) {
  // Obtain boundary conditions
  *ifail = neBEMBoundaryConditions();
}

// -------------------------------------------------------------------

// Garfield procedure providing information on a volume
// TODO
/*
extern void bemvol_(int* volref, int* shape,
                    int* material, double* epsilon,
                    double* potential, double* charge,
                    int* boundarytype,
                    int* ifail);
*/

int neBEMVolumeDescription(int volref, int* shape, int* material,
                           double* epsilon, double* potential, double* charge,
                           int* boundarytype) {
  int ifail = 0, vv = volref + 1;
  // TODO!
  // bemvol_(&vv, shape, material, epsilon, potential, charge, boundarytype,
  // &ifail);
  return ifail;
}

// -------------------------------------------------------------------

void bemsol_(int* ifail) {
  // Garfield calls this procedure to ask neBEM to solve or to retrieve
  *ifail = neBEMSolve();
  // no check here? what is the value of ifail?
  // if(fstatus != 0)
  // 	{
  // 	neBEMMessage("bemsol_ - bemsol");
  // 	return -1;
  // 	}
}

// -------------------------------------------------------------------

// Garfield procedure that tells in which volume a point is located
// TODO!
// extern void celvpt_(double* x, double* y, double* z, int* volref);

int neBEMVolumePoint(double x, double y, double z) {
  // Work out in which volume the point is located
  int volref;
  double xx = x * 100.0, yy = y * 100.0, zz = z * 100.0;
  // celvpt_(&xx, &yy, &zz, &volref);
  return volref - 1;
}

// -------------------------------------------------------------------

// Garfield procedure returning the primitives for a volume
// TODO!
// extern void bemvpr_(int* volume, int* nprim, int* primlist);
void neBEMVolumePrimitives(int volume, int* nprim, int primlist[]) {
  // bemvpr_(&volume, nprim, primlist);
}

// -------------------------------------------------------------------

// Garfield procedure returning the primitives with the same boundaries
// TODO!
// extern void bemsbc_(int* prim1, int* nprim, int* primlist);
void neBEMSameBoundaryConditions(int prim1, int* nprim, int primlist[]) {
  // bemsbc_(&prim1, nprim, primlist);
}

// -------------------------------------------------------------------

// Weighting field: prepare, calculate, delete

void bemrqw_(int* nprim, int* primlist, int* id) {
  int idwf = neBEMPrepareWeightingField((*nprim), primlist);
  *id = idwf;
}

void bemdlw_(int* id) { neBEMDeleteWeightingField(*id); }

// -------------------------------------------------------------------

void ionbem_(float* x, float* y, float* z, float* ex, float* ey, float* ez,
             int* id) {
  // Check location
  int ivol = neBEMVolumePoint(*x / 100, *y / 100, *z / 100);
  if (ivol > -1) {
    *ex = 0.0;
    *ey = 0.0;
    *ez = 0.0;
    return;
  }

  // Computes weighting field id at position (x,y,z)
  Point3D fieldpt;
  fieldpt.X = (double)*x / 100.0;
  fieldpt.Y = (double)*y / 100.0;
  fieldpt.Z = (double)*z / 100.0;

  // Declare the output
  Vector3D field;

  // Compute the field
  int fstatus = neBEMWeightingField(&fieldpt, &field, *id);
  if (fstatus != 0) {
    printf("ionbem: Status from neBEMWeightingField = %d\n", fstatus);
    return;
  }

  // Convert units and precision
  *ex = (float)(field.X) / 100.0;
  *ey = (float)(field.Y) / 100.0;
  *ez = (float)(field.Z) / 100.0;
}

// -------------------------------------------------------------------

void efcbem_(float* x, float* y, float* z, float* ex, float* ey, float* ez,
             float* v, int* opt, int* loc) {
  // Default value
  *ex = 0;
  *ey = 0;
  *ez = 0;
  *v = 0;
  *loc = 0;

  // Check location
  int ivol = neBEMVolumePoint(*x / 100, *y / 100, *z / 100);
  if (ivol > -1) {
    const int MXWIRE = 2000;
    *loc = 2 * MXWIRE + ivol + 1;

    int shape, material, boundarytype;
    double epsilon, potential, charge;
    neBEMVolumeDescription(ivol, &shape, &material, &epsilon, &potential,
                           &charge, &boundarytype);
    if (boundarytype == 1) {
      *v = potential;
      return;
    }
  }

  // Construct a point
  Point3D fieldpt;
  fieldpt.X = (double)*x / 100.0;
  fieldpt.Y = (double)*y / 100.0;
  fieldpt.Z = (double)*z / 100.0;
  //  printf("efcbem:\tcalled for (%g,%g,%g) m\n",fieldpt.X, fieldpt.Y,
  //  fieldpt.Z);

  // Declare the output
  double potential;
  Vector3D field;

  // Compute the field
  int fstatus = neBEMField(&fieldpt, &potential, &field);
  if (fstatus != 0) {
    printf("efcbem: Status from neBEMField = %d\n", fstatus);
    *loc = -10;
    return;
  }
  //  printf("\tField = (%g,%g,%g) V/m, V = %g V\n", field.X, field.Y, field.Z,
  //  potential);
  // Convert units and precision
  *ex = (float)(field.X) / 100.0;
  *ey = (float)(field.Y) / 100.0;
  *ez = (float)(field.Z) / 100.0;
  *v = (float)potential;
}

// -------------------------------------------------------------------

void bemvoq_(int* volume, double* q) { *q = neBEMVolumeCharge(*volume); }

// -------------------------------------------------------------------

void bemend_() {
  // Garfield calls this procedure when terminating

  int fstatus = neBEMEnd();
  if (fstatus != 0) {
    neBEMMessage("bemend_ - bemend");
    // return -1;
  }
}

// -------------------------------------------------------------------

int neBEMGetState(void) { return neBEMState; }

// -------------------------------------------------------------------

void bemqst_(int* status) { *status = neBEMGetState(); }
