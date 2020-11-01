/*
(c) 2005, Supratik Mukhopadhayay, Nayana Majumdar
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "neBEMInterface.h"
#include "Isles.h"
#include "NR.h"
#include "Vector.h"
#include "neBEM.h"

#ifdef __cplusplus
namespace {
  static constexpr double InvFourPiEps0 = 1. / MyFACTOR;
}

namespace neBEM {
#endif

// Weighting field function (WtPFAtPoint) has not been modified!
// Should be merged with PFAtPoint function.
// Check the notes written ahead of the weighting field function.

// Potential per unit charge density on an element
double GetPotential(int ele, Point3D *localP) {
  double value;

  switch ((EleArr + ele - 1)->G.Type) {
    case 4:  // rectangular element
      value = RecPot(ele, localP);
      break;
    case 3:  // triangular element
      value = TriPot(ele, localP);
      break;
    case 2:  // linear (wire) element
      value = WirePot(ele, localP);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends

  return (value);
}  // end of GetPotential

// Potential due to unit charge density on this element
double RecPot(int ele, Point3D *localP) {
  if (DebugLevel == 301) {
    printf("In RecPot ...\n");
  }

  double Pot;
  Vector3D Field;
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  double diag = sqrt(a * a + b * b);  // diagonal of the element

  // distance of field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  if (dist >= FarField * diag)  // all are distances and, hence, +ve
  {
    double dA = a * b;
    Pot = dA / dist;
  } else {
    // normalize distances by `a' while sending - likely to improve accuracy
    int fstatus =
        ExactRecSurf(xpt / a, ypt / a, zpt / a, -0.5, -(b / a) / 2.0,
                     0.5, (b / a) / 2.0, &Pot, &Field);
    if (fstatus)  // non-zero
    {
      printf("problem in computing Potential of rectangular element ... \n");
      printf("a: %lg, b: %lg, X: %lg, Y: %lg, Z: %lg\n", a, b, xpt, ypt, zpt);
      // printf("returning ...\n");
      // return -1; void function at present
    }
    Pot *= a;  // rescale Potential - cannot be done outside because of the `if'
  }

#ifdef __cplusplus
  return Pot * InvFourPiEps0;
#else
  return (Pot / MyFACTOR);
#endif
}  // end of RecPot

// Potential due to unit charge density on a triangular element
double TriPot(int ele, Point3D *localP) {
  if (DebugLevel == 301) {
    printf("In TriPot ...\n");
  }

  double Pot;
  Vector3D Field;
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  // distance of field point from element centroid
  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  // largest side (hypotenuse) of the element
  double diag = sqrt(a * a + b * b);  
  
  const double xm = xpt - a / 3.;
  const double zm = zpt - b / 3.;
  double dist = sqrt(xm * xm + ypt * ypt + zm * zm);

  if (dist >= FarField * diag) {
    double dA = 0.5 * a * b;  // area of the triangular element
    Pot = dA / dist;
  } else {
    int fstatus = ExactTriSurf(b / a, xpt / a, ypt / a, zpt / a, &Pot, &Field);
    if (fstatus) { // non-zero
      printf("problem in computing Potential of triangular element ... \n");
      printf("a: %lg, b: %lg, X: %lg, Y: %lg, Z: %lg\n", a, b, xpt, ypt, zpt);
      // printf("returning ...\n");
      // return -1; void function at present
    }
    Pot *= a;  // rescale Potential
  }

#ifdef __cplusplus
  return Pot * InvFourPiEps0;
#else
  return (Pot / MyFACTOR);
#endif
}  // end of TriPot

// Potential due to unit charge density on this element
// arguments not normalized yet
double WirePot(int ele, Point3D *localP) {
  if (DebugLevel == 301) {
    printf("In WirePot ...\n");
  }

  double Pot;
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double rW = (EleArr + ele - 1)->G.LX;
  double lW = (EleArr + ele - 1)->G.LZ;

  // field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  if (dist >= FarField * lW)  // all are distances and, hence, +ve
  {
    double dA = 2.0 * ST_PI * rW * lW;
    // Pot = ApproxP_W(rW, lW, X, Y, Z, 1);
    Pot = dA / dist;
  } else if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST) &&
             (fabs(zpt) < MINDIST)) {
    Pot = ExactCentroidalP_W(rW, lW);
  } else if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST)) {
    Pot = ExactAxialP_W(rW, lW, localP->Z);
  } else {
    Pot = ExactThinP_W(rW, lW, xpt, ypt, zpt);
  }

#ifdef __cplusplus
  return Pot * InvFourPiEps0;
#else
  return (Pot / MyFACTOR);
#endif
}  // end of WirePot

// Flux per unit charge density on an element returned as globalF
// in the global coordiante system
void GetFluxGCS(int ele, Point3D *localP, Vector3D *globalF) {
  Vector3D localF;

  switch ((EleArr + ele - 1)->G.Type) {
    case 4:  // rectangular element
      RecFlux(ele, localP, &localF);
      break;
    case 3:  // triangular element
      TriFlux(ele, localP, &localF);
      break;
    case 2:  // linear (wire) element
      WireFlux(ele, localP, &localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends

  (*globalF) = RotateVector3D(&localF, &(EleArr + ele - 1)->G.DC, local2global);
}  // end of GetFluxGCS

// Flux per unit charge density on an element returned as localF
// in the local coordiante system
void GetFlux(int ele, Point3D *localP, Vector3D *localF) {
  switch ((EleArr + ele - 1)->G.Type) {
    case 4:  // rectangular element
      RecFlux(ele, localP, localF);
      break;
    case 3:  // triangular element
      TriFlux(ele, localP, localF);
      break;
    case 2:  // linear (wire) element
      WireFlux(ele, localP, localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends
}  // end of GetFlux

// local coord flux localF per unit charge density on one rectangular element
void RecFlux(int ele, Point3D *localP, Vector3D *localF) {
  if (DebugLevel == 301) {
    printf("In RecFlux ...\n");
  }
  
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  double diag = sqrt(a * a + b * b);  // diagonal of the element

  // distance of field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  // no re-scaling necessary for `E'
  if (dist >= FarField * diag) {
    const double f = a * b / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    double Pot;
    int fstatus =
        ExactRecSurf(xpt / a, ypt / a, zpt / a, -0.5, -(b / a) / 2.0,
                     0.5, (b / a) / 2.0, &Pot, localF);
    if (fstatus) { // non-zero
      printf("problem in computing flux of rectangular element ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
  }

#ifdef __cplusplus
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
}  // end of RecFlux

// local coord flux per unit charge density on a triangluar element
void TriFlux(int ele, Point3D *localP, Vector3D *localF) {
  if (DebugLevel == 301) {
    printf("In TriFlux ...\n");
  }

  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  double diag = sqrt(a * a + b * b);  // diagonal of the element

  // printf("In TriFlux\n");
  // printf("a: %lg, b: %lg, X: %lg, Y: %lg, Z: %lg\n", a, b, X, Y, Z);

  // distance of field point from element centroid
  const double xm = xpt - a / 3.;
  const double zm = zpt - b / 3.;
  double dist = sqrt(xm * xm + ypt * ypt + zm * zm);

  if (dist >= FarField * diag) {
    const double f = 0.5 * a * b / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    double Pot;
    int fstatus = ExactTriSurf(b / a, xpt / a, ypt / a, zpt / a, &Pot, localF);
    // fstatus = ApproxTriSurf(b/a, X/a, Y/a, Z/a, 5000, 5000, &Pot, &Flux);
    if (fstatus) { // non-zero
      printf("problem in computing flux of triangular element ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
  }

#ifdef __cplusplus
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
  // printf("Pot: %lg, Ex: %lg, Ey: %lg, Ez: %lg\n",
  // Pot, localF.X, localF.Y, localF.Z);
  // printf("Out of TriFlux\n");
}  // end of TriFlux

// local coord flux per unit charge density on a wire element
void WireFlux(int ele, Point3D *localP, Vector3D *localF) {
  if (DebugLevel == 301) {
    printf("In WireFlux ...\n");
  }

  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double rW = (EleArr + ele - 1)->G.LX;
  double lW = (EleArr + ele - 1)->G.LZ;

  // field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  if (dist >= FarField * lW) {
    const double f = 2.0 * ST_PI * rW * lW / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST)) {
      localF->X = localF->Y = 0.0;
    } else {
      localF->X = ExactThinFX_W(rW, lW, xpt, ypt, zpt);

      localF->Y = ExactThinFY_W(rW, lW, xpt, ypt, zpt);
    }

    // Ez
    localF->Z = ExactThinFZ_W(rW, lW, xpt, ypt, zpt);
  }

#ifdef __cplusplus
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
}  // end of WireFlux

/* PFAtPoint without multi-threading (borrowed from 1.8.15)
// do not erase this redundant function - the multi-threaded version is quite
// complex and this may work as a fall-back option.
// Gives three components of the total Potential and flux in the global
// coordinate system due to all the elements
int PFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  double xfld = globalP->X;
  double yfld = globalP->Y;
  double zfld = globalP->Z;
  Point3D fldpt;
  fldpt.X = xfld; fldpt.Y = yfld; fldpt.Z = zfld;
  double TransformationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0}};

  // Compute Potential and field at different locations
  *Potential = globalF->X = globalF->Y = globalF->Z = 0.0;

  // Effects due to base primitives and their repetitions are considered in the
  // local coordinate system of the primitive (or element), while effects due to
  // mirror elements and their repetitions are considered in the global
  // coordinate system (GCS). This works because the direction cosines of a
  // primitive (and its elements) and those of its repetitions are the same.
  // As a result, we can do just one transformation from local to global at the
  // end of calculations related to a primitive. This can save substantial
  // computation if a discretized version of the primitive is being used since
  // we avoid one unnecessary transformation for each element that comprises a
  // primitive.
  // Begin with primitive description of the device
  double tmpPot;
  Vector3D tmpF;

  for(unsigned int primsrc = 1; primsrc <= NbPrimitives; ++primsrc) {
    double xpsrc = PrimOriginX[primsrc];
    double ypsrc = PrimOriginY[primsrc];
    double zpsrc = PrimOriginZ[primsrc];

    Vector3D localF;
    localF.X = localF.Y = localF.Z = 0.0;

    // Set up transform matrix for this primitive, which is also the same for
    // all the elements belonging to this primitive
    TransformationMatrix[0][0] = PrimDC[primsrc].XUnit.X;
    TransformationMatrix[0][1] = PrimDC[primsrc].XUnit.Y;
    TransformationMatrix[0][2] = PrimDC[primsrc].XUnit.Z;
    TransformationMatrix[1][0] = PrimDC[primsrc].YUnit.X;
    TransformationMatrix[1][1] = PrimDC[primsrc].YUnit.Y;
    TransformationMatrix[1][2] = PrimDC[primsrc].YUnit.Z;
    TransformationMatrix[2][0] = PrimDC[primsrc].ZUnit.X;
    TransformationMatrix[2][1] = PrimDC[primsrc].ZUnit.Y;
    TransformationMatrix[2][2] = PrimDC[primsrc].ZUnit.Z;

    // The total influence is due to primitives on the basic device and due to
    // virtual primitives arising out of repetition, reflection etc and not
    // residing on the basic device

    { // basic device
      Point3D localPP; // point primitive
      // point translated to the ECS origin, but axes direction global
      { // Rotate point3D from global to local system
        double InitialVector[3] = {xfld - xpsrc, yfld - ypsrc, zfld - zpsrc};
        double FinalVector[3] = {0., 0., 0.};
        for (int i = 0; i < 3; ++i) {
          for (int j = 0 ; j < 4; ++j) {
            FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
          }
        }
        localPP.X = FinalVector[0];
        localPP.Y = FinalVector[1];
        localPP.Z = FinalVector[2];
      } // Point3D rotated

      // evaluate possibility whether primitive influence is accurate enough
      // This could be based on localPP and the subtended solid angle
      int PrimOK = 0;

      if (PrimOK) { // if 1, then only primitive influence will be considered
        // Potential and flux (local system) due to base primitive
        GetPrimPF(primsrc, &localPP, &tmpPot, &tmpF);
        (*Potential) += (AvChDen[primsrc]+AvAsgndChDen[primsrc]) * tmpPot;
        localF.X += (AvChDen[primsrc]+AvAsgndChDen[primsrc]) * tmpF.X;
        localF.Y += (AvChDen[primsrc]+AvAsgndChDen[primsrc]) * tmpF.Y;
        localF.Z += (AvChDen[primsrc]+AvAsgndChDen[primsrc]) * tmpF.Z;
        if (DebugLevel == 301) {
          printf("PFAtPoint base primitive =>\n");
          printf("primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n",
                  primsrc, localPP.X, localPP.Y, localPP.Z);
          printf("primsrc: %d, Pot: %lg, Fx: %lg, Fx: %lg, Fz: %lg\n",
                  primsrc, tmpPot, tmpF.X, tmpF.Y, tmpF.Z);
          printf("primsrc: %d, SumPot: %lg, SumFx: %lg, SumFy: %lg, SumFz: %lg\n",
                  primsrc, *Potential, localF.X, localF.Y, localF.Z);
        }
      } else {
        // element influence

        for (unsigned int ele = ElementBgn[primsrc]; ele <= ElementEnd[primsrc]; ++ele) {
          double xsrc = (EleArr+ele-1)->G.Origin.X;
          double ysrc = (EleArr+ele-1)->G.Origin.Y; 
          double zsrc = (EleArr+ele-1)->G.Origin.Z;
          // Rotate point3D from global to local system; 
          // matrix as for primitive 
          double InitialVector[3] = {xfld - xsrc, yfld - ysrc, zfld - zsrc};
          double FinalVector[3] = {0., 0., 0.};
          for (int i = 0; i < 3; ++i) {
            for (int j = 0 ; j < 3; ++j) {
              FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
            }
          }
          Point3D localPE; // point element
          localPE.X = FinalVector[0];
          localPE.Y = FinalVector[1];
          localPE.Z = FinalVector[2];
          // Potential and flux (local system)
          GetPF(ele, &localPE, &tmpPot, &tmpF);
          (*Potential) += ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned) * tmpPot;
          localF.X += ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned) * tmpF.X;
          localF.Y += ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned) * tmpF.Y;
          localF.Z += ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned) * tmpF.Z;
          if (DebugLevel == 301) {
            printf("PFAtPoint base primitive =>\n");
            printf("ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n", 
                   ele, localPE.X, localPE.Y, localPE.Z); 
            printf("ele: %d, Pot: %lg, Fx: %lg, Fx: %lg, Fz: %lg\n", 
                   ele, tmpPot, tmpF.X, tmpF.Y, tmpF.Z); 
            printf("ele: %d, SumPot: %lg, SumFx: %lg, SumFy: %lg, SumFz: %lg\n", 
                   ele, *Potential, localF.X, localF.Y, localF.Z);
          }
        } // for all the elements on this primsrc primitive 
      } // else elements influence 
    } // basic device ends

    if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] || MirrorTypeZ[primsrc]) {
      // Mirror effect of base primitives
      printf("Mirror may not be correctly implemented ...\n");
      exit(0);
    } // Mirror effect ends

    // Flux due to repeated primitives
    if ((PeriodicTypeX[primsrc] == 1) || (PeriodicTypeY[primsrc] == 1) || 
        (PeriodicTypeZ[primsrc] == 1)) {
      if (PeriodicInX[primsrc] || PeriodicInY[primsrc] || PeriodicInZ[primsrc]) {
        for (int xrpt = -PeriodicInX[primsrc]; xrpt <= PeriodicInX[primsrc]; ++xrpt) {
          double XPOfRpt = xpsrc + XPeriod[primsrc] * (double)xrpt;
          for (int yrpt = -PeriodicInY[primsrc]; yrpt <= PeriodicInY[primsrc]; ++yrpt) {
            double YPOfRpt = ypsrc + YPeriod[primsrc] * (double)yrpt;
            for (int zrpt = -PeriodicInZ[primsrc]; zrpt <= PeriodicInZ[primsrc]; ++zrpt) {
              double ZPOfRpt = zpsrc + ZPeriod[primsrc] * (double)zrpt;
              if ((xrpt == 0) && (yrpt == 0) && (zrpt == 0)) {
                continue; // this is the base device
              }
              { // basic primitive repeated 
                Point3D localPPR;
                {
                  double InitialVector[3] = {xfld - XPOfRpt, yfld - YPOfRpt, zfld - ZPOfRpt};
                  double FinalVector[3] = {0., 0., 0.};
                  for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                      FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
                    }
                  }
                  localPPR.X = FinalVector[0];
                  localPPR.Y = FinalVector[1];
                  localPPR.Z = FinalVector[2];
                } // Point3D rotated
                int PrimOK = 0;
                // consider primitive representation accurate enough if it is
                // repeated and beyond PrimAfter repetitions. 
                if (PrimAfter == 0) {
                  // If PrimAfter is zero, PrimOK is always zero
                  PrimOK = 0;
                  // and the following is not evaluated at all.
                } else if (abs(xrpt) > PrimAfter && abs(yrpt) > PrimAfter) {
                  PrimOK = 1;
                }
                if (PrimOK) {
                  // use primitive representation
                  // Potential and flux (local system) due to repeated primitive 
                  GetPrimPF(primsrc, &localPPR, &tmpPot, &tmpF);
                  (*Potential) += (AvChDen[primsrc]+AvAsgndChDen[primsrc]) * tmpPot;
                  localF.X += (AvChDen[primsrc]+AvAsgndChDen[primsrc]) * tmpF.X; 
                  localF.Y += (AvChDen[primsrc]+AvAsgndChDen[primsrc]) * tmpF.Y; 
                  localF.Z += (AvChDen[primsrc]+AvAsgndChDen[primsrc]) * tmpF.Z; 
                  if (DebugLevel == 301) {
                    printf("primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal: %lg\n", 
                           primsrc, localPPR.X, localPPR.Y, localPPR.Z);
                    printf("primsrc: %d, Pot: %lg, Fx: %lg, Fy: %lg, Fz: %lg\n", 
                           primsrc,
                           tmpPot * (AvChDen[primsrc]+AvAsgndChDen[primsrc]),
                           tmpF.X * (AvChDen[primsrc]+AvAsgndChDen[primsrc]),
                           tmpF.Y * (AvChDen[primsrc]+AvAsgndChDen[primsrc]),
                           tmpF.Z * (AvChDen[primsrc]+AvAsgndChDen[primsrc]));
                    printf("primsrc: %d, SumPot: %lg, SumFx: %lg, SumFy: %lg, SumFz: %lg\n", 
                           primsrc, *Potential, localF.X, localF.Y, localF.Z);
                  }
                } else {
                  // use discretized representation of a repeated primitive
                  Point3D localPER; // point element repeated
                  for (unsigned int ele = ElementBgn[primsrc]; ele <= ElementEnd[primsrc]; ++ele) {
                    double xsrc = (EleArr+ele-1)->G.Origin.X;
                    double ysrc = (EleArr+ele-1)->G.Origin.Y; 
                    double zsrc = (EleArr+ele-1)->G.Origin.Z;
                                                                
                    double XOfRpt = xsrc + XPeriod[primsrc] * (double)xrpt; 
                    double YOfRpt = ysrc + YPeriod[primsrc] * (double)yrpt;
                    double ZOfRpt = zsrc + ZPeriod[primsrc] * (double)zrpt;
                                                                
                    Rotate from global to local system  
                    double InitialVector[3] = {xfld - XOfRpt, yfld - YOfRpt, zfld - ZOfRpt}; 
                    double FinalVector[3] = {0., 0., 0.};
                    for (int i = 0; i < 3; ++i) {
                      for (int j = 0; j < 3; ++j) {
                        FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
                      }
                    }
                    localPER.X = FinalVector[0]; 
                    localPER.Y = FinalVector[1]; 
                    localPER.Z = FinalVector[2];

                    GetPF(ele, &localPER, &tmpPot, &tmpF);
                    *Potential) += ((EleArr+ele-1)->Solution + (EleArr+ele-1)->Assigned) * tmpPot;
                    localF.X += ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned) * tmpF.X;
                    localF.Y += ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned) * tmpF.Y;
                    localF.Z += ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned) * tmpF.Z;
                    if (DebugLevel == 301) {
                      printf("primsrc: %d, ele: %d, xlocal: %lg, ylocal: %lg, zlocal: %lg\n", 
                             primsrc, ele, localPER.X, localPER.Y, localPER.Z); 
                      printf("primsrc: %d, Pot: %lg, Fx: %lg, Fy: %lg, Fz: %lg\n", primsrc 
                             tmpPot * ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned),
                             tmpF.X * ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned),
                             tmpF.Y * ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned),
                             tmpF.Z * ((EleArr+ele-1)->Solution+(EleArr+ele-1)->Assigned));
                      printf("primsrc: %d, SumPot: %lg, SumFx: %lg, SumFy: %lg, SumFz: %lg\n", 
                             primsrc, *Potential, localF.X, localF.Y, localF.Z);
                    }
                  } // for all elements on this primitive
                }
              } // repetition of basic primitive
            } // for zrpt
          } // for yrpt
        } // for xrpt
      } // PeriodicInX || PeriodicInY || PeriodicInZ
    } // PeriodicType == 1

    tmpF = RotateVector3D(&localF, &PrimDC[primsrc], local2global);
    globalF->X += tmpF.X;
    globalF->Y += tmpF.Y;
    globalF->Z += tmpF.Z;
  } // for all primitives: basic device, mirror reflections and repetitions

  // ExactPointP and ExactPointF should also have an ExactPointPF
  // Similarly for area and volume element related functions
  // since there is no intermediate function that interfaces ExactPointP etc
  // division by MyFACTOR is necessary
  for (unsigned int point = 1; point <= NbPtsKnCh; ++point) {
    tmpPot = ExactPointP(&(PtKnChArr+point-1)->P, globalP);
    (*Potential) += (PtKnChArr+point-1)->Assigned * tmpPot / MyFACTOR;
    ExactPointF(&(PtKnChArr+point-1)->P, globalP, &tmpF);
    globalF->X += (PtKnChArr+point-1)->Assigned * tmpF.X / MyFACTOR;
    globalF->Y += (PtKnChArr+point-1)->Assigned * tmpF.Y / MyFACTOR;
    globalF->Z += (PtKnChArr+point-1)->Assigned * tmpF.Z / MyFACTOR;
  }

  for (unsigned int line = 1; line <= NbLinesKnCh; ++line) {
    (*Potential) += 0.0;
    globalF->X += 0.0;
    globalF->Y += 0.0;
    globalF->Z += 0.0;
  }

  for (unsigned int area = 1; area <= NbAreasKnCh; ++area) {
    (*Potential) += 0.0;
    globalF->X += 0.0;
    globalF->Y += 0.0;
    globalF->Z += 0.0;
  }

  for (unsigned int vol = 1; vol <= NbVolsKnCh; ++vol) {
    (*Potential) += 0.0;
    globalF->X += 0.0;
    globalF->Y += 0.0;
    globalF->Z += 0.0;
  }

  (*Potential) += VSystemChargeZero; // respect total system charge constraint
  return(0);
} // end of PFAtPoint
// do not erase this redundant function.
PFAtPoint borrowed from V1.8.15 ends here */

// Gives three components of the total Potential and flux in the global
// coordinate system due to all the interface elements and all known charges.
// It may be interesting to inspect the relative influence of these two factors.
int PFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  double ElePot;
  Vector3D EleglobalF;
  int fstatus = ElePFAtPoint(globalP, &ElePot, &EleglobalF);
  if (fstatus) {
    printf(
        "Problem in ElePFAtPoint being called from PFAtPoint ... returning\n");
    return (-1);
  }
  *Potential = ElePot;
  globalF->X = EleglobalF.X;
  globalF->Y = EleglobalF.Y;
  globalF->Z = EleglobalF.Z;

  if (OptKnCh) {
    double KnChPot;
    Vector3D KnChglobalF;
    fstatus = KnChPFAtPoint(globalP, &KnChPot, &KnChglobalF);
    if (fstatus) {
      printf(
          "Problem in KnChPFAtPoint being called from PFAtPoint ... "
          "returning\n");
      return (-1);
    }
    *Potential += KnChPot;
    globalF->X += KnChglobalF.X;
    globalF->Y += KnChglobalF.Y;
    globalF->Z += KnChglobalF.Z;
  }

  return 0;
}  // PFAtPoint ends

// Gives three components of the total Potential and flux in the global
// coordinate system only due to all the interface elements.
// Multi-threading implemented in the following routine
int ElePFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;

  const double xfld = globalP->X;
  const double yfld = globalP->Y;
  const double zfld = globalP->Z;

  // Compute Potential and field at different locations
  *Potential = globalF->X = globalF->Y = globalF->Z = 0.0;

  // Effects due to base primitives and their repetitions are considered in the
  // local coordinate system of the primitive (or element), while effects due to
  // mirror elements and their repetitions are considered in the global
  // coordinate system (GCS). This works because the direction cosines of a
  // primitive (and its elements) and those of its repetitions are the same.
  // As a result, we can do just one transformation from local to global at the
  // end of calculations related to a primitive. This can save substantial
  // computation if a discretized version of the primitive is being used since
  // we avoid one unnecessary transformation for each element that comprises a
  // primitive.
  // Begin with primitive description of the device

  // Scope in OpenMP: Variables in the global data space are accessible to all
  // threads, while variables in a thread's private space is accessible to the
  // thread only (there are several variations - copying outside region etc)
  // Field point remains the same - kept outside private
  // source point changes with change in primitive - private
  // TransformationMatrix changes - kept within private (Note: matrices with
  // fixed dimensions can be maintained, but those with dynamic allocation
  // can not).
  double *pPot = dvector(1, NbPrimitives);
  // Field components in LCS for a primitive and its other incarnations.
  double *plFx = dvector(1, NbPrimitives);  
  double *plFy = dvector(1, NbPrimitives);
  double *plFz = dvector(1, NbPrimitives);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    pPot[prim] = plFx[prim] = plFy[prim] = plFz[prim] = 0.0;
  }

#ifdef _OPENMP
  int tid = 0, nthreads = 1;
  #pragma omp parallel private(tid, nthreads)
#endif
  {

#ifdef _OPENMP
    if (dbgFn) {
      tid = omp_get_thread_num();
      if (tid == 0) {
        nthreads = omp_get_num_threads();
        printf("PFAtPoint computation with %d threads\n", nthreads);
      }
    }
#endif
// by default, nested parallelization is off in C
#ifdef _OPENMP
#pragma omp for 
#endif
    for (int primsrc = 1; primsrc <= NbPrimitives; ++primsrc) {
      if (dbgFn) {
        printf("Evaluating effect of primsrc %d using on %lg, %lg, %lg\n",
               primsrc, xfld, yfld, zfld);
        fflush(stdout);
      }

      const double xpsrc = PrimOriginX[primsrc];
      const double ypsrc = PrimOriginY[primsrc];
      const double zpsrc = PrimOriginZ[primsrc];

      // Field in the local frame.
      double lFx = 0.;
      double lFy = 0.;
      double lFz = 0.;

      // Set up transform matrix for this primitive, which is also the same 
      // for all the elements belonging to this primitive.
      double TransformationMatrix[3][3];
      TransformationMatrix[0][0] = PrimDC[primsrc].XUnit.X;
      TransformationMatrix[0][1] = PrimDC[primsrc].XUnit.Y;
      TransformationMatrix[0][2] = PrimDC[primsrc].XUnit.Z;
      TransformationMatrix[1][0] = PrimDC[primsrc].YUnit.X;
      TransformationMatrix[1][1] = PrimDC[primsrc].YUnit.Y;
      TransformationMatrix[1][2] = PrimDC[primsrc].YUnit.Z;
      TransformationMatrix[2][0] = PrimDC[primsrc].ZUnit.X;
      TransformationMatrix[2][1] = PrimDC[primsrc].ZUnit.Y;
      TransformationMatrix[2][2] = PrimDC[primsrc].ZUnit.Z;

      // The total influence is due to primitives on the basic device and due to
      // virtual primitives arising out of repetition, reflection etc and not
      // residing on the basic device

      // basic primitive

      // Evaluate possibility whether primitive influence is accurate enough.
      // This could be based on localPP and the subtended solid angle.
      // If 1, then only primitive influence will be considered.
      int PrimOK = 0;

      if (PrimOK) {
        // Only primitive influence will be considered
        // Potential and flux (local system) due to base primitive
        double tmpPot;
        Vector3D tmpF;
        // Rotate point from global to local system
        double InitialVector[3] = {xfld - xpsrc, yfld - ypsrc, zfld - zpsrc};
        double FinalVector[3] = {0., 0., 0.};
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
          }
        }
        Point3D localPP;
        localPP.X = FinalVector[0];
        localPP.Y = FinalVector[1];
        localPP.Z = FinalVector[2];
        GetPrimPF(primsrc, &localPP, &tmpPot, &tmpF);
        const double qpr = AvChDen[primsrc] + AvAsgndChDen[primsrc];
        pPot[primsrc] += qpr * tmpPot;
        lFx += qpr * tmpF.X;
        lFy += qpr * tmpF.Y;
        lFz += qpr * tmpF.Z;
        // if(DebugLevel == 301)
        if (dbgFn) {
          printf("PFAtPoint base primitive =>\n");
          printf("primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n",
                 primsrc, localPP.X, localPP.Y, localPP.Z);
          printf("primsrc: %d, Pot: %lg, Fx: %lg, Fx: %lg, Fz: %lg\n",
                 primsrc, tmpPot, tmpF.X, tmpF.Y, tmpF.Z);
          printf("primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
                 primsrc, pPot[primsrc], lFx, lFy, lFz);
          fflush(stdout);
          // exit(-1);
        }
      } else {
        // Need to consider element influence.
        double tPot;
        Vector3D tF;
        double ePot = 0.;
        Vector3D eF;
        eF.X = 0.0;
        eF.Y = 0.0;
        eF.Z = 0.0;
        const int eleMin = ElementBgn[primsrc];
        const int eleMax = ElementEnd[primsrc];
        for (int ele = eleMin; ele <= eleMax; ++ele) {
          const double xsrc = (EleArr + ele - 1)->G.Origin.X;
          const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
          const double zsrc = (EleArr + ele - 1)->G.Origin.Z;
          // Rotate from global to local system; matrix as for primitive
          double InitialVector[3] = {xfld - xsrc, yfld - ysrc, zfld - zsrc};
          double FinalVector[3] = {0., 0., 0.};
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
            }
          }
          Point3D localPE;
          localPE.X = FinalVector[0];
          localPE.Y = FinalVector[1];
          localPE.Z = FinalVector[2];

          // Potential and flux (local system) due to base primitive
          GetPF(ele, &localPE, &tPot, &tF);
          const double qel = (EleArr + ele - 1)->Solution + (EleArr + ele - 1)->Assigned;
          ePot += qel * tPot;
          eF.X += qel * tF.X;
          eF.Y += qel * tF.Y;
          eF.Z += qel * tF.Z;
          // if(DebugLevel == 301)
          if (dbgFn) {
            printf("PFAtPoint base primitive:%d\n", primsrc);
            printf("ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n", ele,
                   localPE.X, localPE.Y, localPE.Z);
            printf(
                "ele: %d, tPot: %lg, tFx: %lg, tFy: %lg, tFz: %lg, Solution: "
                "%g\n",
                ele, tPot, tF.X, tF.Y, tF.Z, qel);
            printf("ele: %d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n", ele,
                   ePot, eF.X, eF.Y, eF.Z);
            fflush(stdout);
          }
        }  // for all the elements on this primsrc primitive

        pPot[primsrc] += ePot;
        lFx += eF.X;
        lFy += eF.Y;
        lFz += eF.Z;
        if (dbgFn) {
          printf(
              "prim%d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n",
              primsrc, ePot, eF.X, eF.Y, eF.Z);
          printf("prim%d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n", primsrc,
                 pPot[primsrc], lFx, lFy, lFz);
          fflush(stdout);
        }
      }  // else elements influence

      // if(DebugLevel == 301)
      if (dbgFn) {
        printf("basic primitive\n");
        printf("primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
               primsrc, pPot[primsrc], lFx, lFy, lFz);
        fflush(stdout);
      }

      if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
          MirrorTypeZ[primsrc]) {  // Mirror effect of base primitives
        printf("Mirror may not be correctly implemented ...\n");
        exit(0);
      } // Mirror effect ends

      // Flux due to repeated primitives
      if ((PeriodicTypeX[primsrc] == 1) || (PeriodicTypeY[primsrc] == 1) ||
          (PeriodicTypeZ[primsrc] == 1)) {
        const int perx = PeriodicInX[primsrc];
        const int pery = PeriodicInY[primsrc];
        const int perz = PeriodicInZ[primsrc];
        if (perx || pery || perz) {
          for (int xrpt = -perx; xrpt <= perx; ++xrpt) {
            const double xShift = XPeriod[primsrc] * (double)xrpt;
            const double XPOfRpt = xpsrc + xShift;
            for (int yrpt = -pery; yrpt <= pery; ++yrpt) {
              const double yShift = YPeriod[primsrc] * (double)yrpt;
              const double YPOfRpt = ypsrc + yShift;
              for (int zrpt = -perz; zrpt <= perz; ++zrpt) {
                const double zShift = ZPeriod[primsrc] * (double)zrpt;
                const double ZPOfRpt = zpsrc + zShift;
                // Skip the basic primitive.
                if ((xrpt == 0) && (yrpt == 0) && (zrpt == 0)) continue;

                // Basic primitive repeated
                int PrimOK = 0;

                // consider primitive representation accurate enough if it is
                // repeated and beyond PrimAfter repetitions.
                if (PrimAfter == 0) {
                  // If PrimAfter is zero, PrimOK is always zero
                  PrimOK = 0;
                } else if ((abs(xrpt) > PrimAfter) && (abs(yrpt) > PrimAfter)) {
                  PrimOK = 1;
                }
                if (PrimOK) {
                  // Use primitive representation
                  // Rotate point from global to local system
                  double InitialVector[3] = {xfld - XPOfRpt, yfld - YPOfRpt, zfld - ZPOfRpt};
                  double FinalVector[3] = {0., 0., 0.};
                  for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                      FinalVector[i] +=
                          TransformationMatrix[i][j] * InitialVector[j];
                    }
                  }
                  Point3D localPPR;
                  localPPR.X = FinalVector[0];
                  localPPR.Y = FinalVector[1];
                  localPPR.Z = FinalVector[2];
                  // Potential and flux (local system) due to repeated
                  // primitive
                  double tmpPot;
                  Vector3D tmpF;
                  GetPrimPF(primsrc, &localPPR, &tmpPot, &tmpF);
                  const double qpr = AvChDen[primsrc] + AvAsgndChDen[primsrc];
                  pPot[primsrc] += qpr * tmpPot;
                  lFx += qpr * tmpF.X;
                  lFy += qpr * tmpF.Y;
                  lFz += qpr * tmpF.Z;
                  // if(DebugLevel == 301)
                  if (dbgFn) {
                    printf(
                        "primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal: "
                        "%lg\n",
                        primsrc, localPPR.X, localPPR.Y, localPPR.Z);
                    printf(
                        "primsrc: %d, Pot: %lg, Fx: %lg, Fy: %lg, Fz: %lg\n",
                        primsrc, tmpPot * qpr, tmpF.X * qpr, tmpF.Y * qpr,
                        tmpF.Z * qpr);
                    printf(
                        "primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: "
                        "%lg\n",
                        primsrc, pPot[primsrc], lFx, lFy, lFz);
                    fflush(stdout);
                  }
                } else {
                  // Use discretized representation of a repeated primitive
                  double tPot;
                  Vector3D tF;
                  double erPot = 0.0;
                  Vector3D erF;
                  erF.X = 0.0;
                  erF.Y = 0.0;
                  erF.Z = 0.0;
                  const int eleMin = ElementBgn[primsrc];
                  const int eleMax = ElementEnd[primsrc];
                  for (int ele = eleMin; ele <= eleMax; ++ele) {
                    const double xrsrc = (EleArr + ele - 1)->G.Origin.X;
                    const double yrsrc = (EleArr + ele - 1)->G.Origin.Y;
                    const double zrsrc = (EleArr + ele - 1)->G.Origin.Z;

                    const double XEOfRpt = xrsrc + xShift;
                    const double YEOfRpt = yrsrc + yShift;
                    const double ZEOfRpt = zrsrc + zShift;
                    // Rotate from global to local system
                    double InitialVector[3] = {xfld - XEOfRpt, yfld - YEOfRpt, zfld - ZEOfRpt};
                    double FinalVector[3] = {0., 0., 0.};
                    for (int i = 0; i < 3; ++i) {
                      for (int j = 0; j < 3; ++j) {
                        FinalVector[i] +=
                            TransformationMatrix[i][j] * InitialVector[j];
                      }
                    }
                    Point3D localPER;
                    localPER.X = FinalVector[0];
                    localPER.Y = FinalVector[1];
                    localPER.Z = FinalVector[2];

                    // Allowed, because all the local coordinates have the
                    // same orientations. Only the origins are mutually
                    // displaced along a line.
                    GetPF(ele, &localPER, &tPot, &tF);
                    const double qel = (EleArr + ele - 1)->Solution + (EleArr + ele - 1)->Assigned;
                    erPot += qel * tPot;
                    erF.X += qel * tF.X;
                    erF.Y += qel * tF.Y;
                    erF.Z += qel * tF.Z;
                    // if(DebugLevel == 301)
                    if (dbgFn) {
                      printf("PFAtPoint base primitive:%d\n", primsrc);
                      printf("ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n",
                             ele, localPER.X, localPER.Y, localPER.Z);
                      printf(
                          "ele: %d, tPot: %lg, tFx: %lg, tFy: %lg, tFz: %lg, "
                          "Solution: %g\n",
                          ele, tPot, tF.X, tF.Y, tF.Z, qel);
                      printf(
                          "ele: %d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n",
                          ele, erPot, erF.X, erF.Y, erF.Z);
                      fflush(stdout);
                    }
                  }  // for all the elements on this primsrc repeated
                     // primitive

                  pPot[primsrc] += erPot;
                  lFx += erF.X;
                  lFy += erF.Y;
                  lFz += erF.Z;
                }  // else discretized representation of this primitive

                // if(DebugLevel == 301)
                if (dbgFn) {
                  printf("basic repeated xrpt: %d. yrpt: %d, zrpt: %d\n",
                         xrpt, yrpt, zrpt);
                  printf(
                      "primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
                      primsrc, pPot[primsrc], lFx, lFy, lFz);
                  fflush(stdout);
                }

                if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
                    MirrorTypeZ[primsrc]) {  // Mirror effect of repeated
                                             // primitives - not parallelized
                  printf(
                      "Mirror not correctly implemented in this version of "
                      "neBEM ...\n");
                  exit(0);

                  double tmpPot;
                  Vector3D tmpF;
                  Point3D srcptp;
                  Point3D localPPRM;  // point primitive repeated mirrored
                  DirnCosn3D DirCos;

                  srcptp.X = XPOfRpt;
                  srcptp.Y = YPOfRpt;
                  srcptp.Z = ZPOfRpt;

                  Point3D fldpt;
                  fldpt.X = xfld;
                  fldpt.Y = yfld;
                  fldpt.Z = zfld;
                  if (MirrorTypeX[primsrc]) {
                    MirrorTypeY[primsrc] = 0;
                    MirrorTypeZ[primsrc] = 0;
                  }
                  if (MirrorTypeY[primsrc]) MirrorTypeZ[primsrc] = 0;

                  if (MirrorTypeX[primsrc]) {
                    localPPRM = ReflectPrimitiveOnMirror(
                        'X', primsrc, srcptp, fldpt,
                        MirrorDistXFromOrigin[primsrc], &DirCos);

                    // check whether primitive description is good enough
                    int PrimOK = 0;
                    if (PrimOK) {
                      GetPrimPFGCS(primsrc, &localPPRM, &tmpPot, &tmpF,
                                   &DirCos);
                      const double qpr = AvChDen[primsrc] + AvAsgndChDen[primsrc];
                      if (MirrorTypeX[primsrc] == 1) {
                        // opposite charge density
                        pPot[primsrc] -= qpr * tmpPot;
                        lFx -= qpr * tmpF.X;
                        lFy -= qpr * tmpF.Y;
                        lFz -= qpr * tmpF.Z;
                      } else if (MirrorTypeX[primsrc] == 2) {
                        // same charge density
                        pPot[primsrc] += qpr * tmpPot;
                        lFx += qpr * tmpF.X;
                        lFy += qpr * tmpF.Y;
                        lFz += qpr * tmpF.Z;
                      }
                    } else {              // consider element representation
                      Point3D localPERM;  // point element repeated mirrored
                      Point3D srcpte;

                      const int eleMin = ElementBgn[primsrc];
                      const int eleMax = ElementEnd[primsrc];
                      for (int ele = eleMin; ele <= eleMax; ++ele) {
                        const double xsrc = (EleArr + ele - 1)->G.Origin.X;
                        const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
                        const double zsrc = (EleArr + ele - 1)->G.Origin.Z;

                        const double XEOfRpt = xsrc + xShift;
                        const double YEOfRpt = ysrc + yShift;
                        const double ZEOfRpt = zsrc + zShift;

                        srcpte.X = XEOfRpt;
                        srcpte.Y = YEOfRpt;
                        srcpte.Z = ZEOfRpt;

                        localPERM = ReflectOnMirror(
                            'X', ele, srcpte, fldpt,
                            MirrorDistXFromOrigin[primsrc], &DirCos);
                        GetPFGCS(ele, &localPERM, &tmpPot, &tmpF, &DirCos);  // force?
                        const double qel = (EleArr + ele - 1)->Solution + (EleArr + ele - 1)->Assigned;
                        if (MirrorTypeX[primsrc] == 1) {
                          // opposite charge density
                          pPot[primsrc] -= qel * tmpPot;
                          lFx -= qel * tmpF.X;
                          lFy -= qel * tmpF.Y;
                          lFz -= qel * tmpF.Z;
                        } else if (MirrorTypeX[primsrc] == 2) {
                          // same charge density
                          pPot[primsrc] += qel * tmpPot;
                          lFx += qel * tmpF.X;
                          lFy += qel * tmpF.Y;
                          lFz += qel * tmpF.Z;
                        }
                      }  // loop for all elements on the primsrc primitive
                    }    // else element representation
                  }      // MirrorTypeX

                  if (MirrorTypeY[primsrc]) {
                    localPPRM = ReflectOnMirror('Y', primsrc, srcptp, fldpt,
                                                MirrorDistYFromOrigin[primsrc],
                                                &DirCos);

                    // check whether primitive description is good enough
                    int PrimOK = 0;
                    if (PrimOK) {
                      GetPrimPFGCS(primsrc, &localPPRM, &tmpPot, &tmpF,
                                   &DirCos);
                      const double qpr = AvChDen[primsrc] + AvAsgndChDen[primsrc];
                      if (MirrorTypeY[primsrc] == 1) {
                        // opposite charge density
                        pPot[primsrc] -= qpr * tmpPot;
                        lFx -= qpr * tmpF.X;
                        lFy -= qpr * tmpF.Y;
                        lFz -= qpr * tmpF.Z;
                      } else if (MirrorTypeY[primsrc] == 2) {
                        // same charge density
                        pPot[primsrc] += qpr * tmpPot;
                        lFx += qpr * tmpF.X;
                        lFy += qpr * tmpF.Y;
                        lFz += qpr * tmpF.Z;
                      }
                    } else {  // consider element representation
                      Point3D localPERM;
                      Point3D srcpte;

                      const int eleMin = ElementBgn[primsrc];
                      const int eleMax = ElementEnd[primsrc];
                      for (int ele = eleMin; ele <= eleMax; ++ele) {
                        const double xsrc = (EleArr + ele - 1)->G.Origin.X;
                        const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
                        const double zsrc = (EleArr + ele - 1)->G.Origin.Z;

                        const double XEOfRpt = xsrc + xShift;
                        const double YEOfRpt = ysrc + yShift;
                        const double ZEOfRpt = zsrc + zShift;

                        srcpte.X = XEOfRpt;
                        srcpte.Y = YEOfRpt;
                        srcpte.Z = ZEOfRpt;

                        localPERM = ReflectOnMirror(
                            'Y', ele, srcpte, fldpt,
                            MirrorDistYFromOrigin[primsrc], &DirCos);
                        GetPFGCS(ele, &localPERM, &tmpPot, &tmpF, &DirCos);
                        const double qel = (EleArr + ele - 1)->Solution + (EleArr + ele - 1)->Assigned;
                        if (MirrorTypeY[primsrc] == 1) {
                          // opposite charge density
                          pPot[primsrc] -= qel * tmpPot;
                          lFx -= qel * tmpF.X;
                          lFy -= qel * tmpF.Y;
                          lFz -= qel * tmpF.Z;
                        } else if (MirrorTypeY[primsrc] == 2) {
                          // same charge density
                          pPot[primsrc] += qel * tmpPot;
                          lFx += qel * tmpF.X;
                          lFy += qel * tmpF.Y;
                          lFz += qel * tmpF.Z;
                        }
                      }  // loop for all elements on the primsrc primitive
                    }    // else element representations
                  }      // MirrorTypeY

                  if (MirrorTypeZ[primsrc]) {
                    localPPRM = ReflectOnMirror('Z', primsrc, srcptp, fldpt,
                                                MirrorDistZFromOrigin[primsrc],
                                                &DirCos);

                    // check whether primitive description is good enough
                    int PrimOK = 0;
                    if (PrimOK) {
                      GetPrimPFGCS(primsrc, &localPPRM, &tmpPot, &tmpF,
                                   &DirCos);
                      const double qpr = AvChDen[primsrc] + AvAsgndChDen[primsrc];
                      if (MirrorTypeZ[primsrc] == 1) {
                        // opposite charge density
                        pPot[primsrc] -= qpr * tmpPot;
                        lFx -= qpr * tmpF.X;
                        lFy -= qpr * tmpF.Y;
                        lFz -= qpr * tmpF.Z;
                      } else if (MirrorTypeZ[primsrc] == 2) {
                        // same charge density
                        pPot[primsrc] += qpr * tmpPot;
                        lFx += qpr * tmpF.X;
                        lFy += qpr * tmpF.Y;
                        lFz += qpr * tmpF.Z;
                      }
                    } else {
                      // elements to be considered
                      Point3D localPERM;
                      Point3D srcpte;

                      const int eleMin = ElementBgn[primsrc];
                      const int eleMax = ElementEnd[primsrc];
                      for (int ele = eleMin; ele <= eleMax; ++ele) {
                        const double xsrc = (EleArr + ele - 1)->G.Origin.X;
                        const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
                        const double zsrc = (EleArr + ele - 1)->G.Origin.Z;

                        const double XEOfRpt = xsrc + xShift;
                        const double YEOfRpt = ysrc + yShift;
                        const double ZEOfRpt = zsrc + zShift;

                        srcpte.X = XEOfRpt;
                        srcpte.Y = YEOfRpt;
                        srcpte.Z = ZEOfRpt;

                        localPERM = ReflectOnMirror(
                            'Z', ele, srcpte, fldpt,
                            MirrorDistZFromOrigin[primsrc], &DirCos);
                        GetPFGCS(ele, &localPERM, &tmpPot, &tmpF, &DirCos);
                        const double qel = (EleArr + ele - 1)->Solution + (EleArr + ele - 1)->Assigned;
                        if (MirrorTypeZ[primsrc] == 1) {
                          // opposite charge density
                          pPot[primsrc] -= qel * tmpPot;
                          lFx -= qel * tmpF.X;
                          lFy -= qel * tmpF.Y;
                          lFz -= qel * tmpF.Z;
                        } else if (MirrorTypeZ[primsrc] == 2) {
                          // same charge density
                          pPot[primsrc] += qel * tmpPot;
                          lFx += qel * tmpF.X;
                          lFy += qel * tmpF.Y;
                          lFz += qel * tmpF.Z;
                        }
                      }  // loop for all elements on the primsrc primitive
                    }    // else consider element representation
                  }      // MirrorTypeZ
                }        // Mirror effect for repeated primitives ends

              }  // for zrpt
            }    // for yrpt
          }      // for xrpt
        }        // PeriodicInX || PeriodicInY || PeriodicInZ
      }          // PeriodicType == 1
      Vector3D localF;
      localF.X = lFx;
      localF.Y = lFy;
      localF.Z = lFz;
      Vector3D tmpF = RotateVector3D(&localF, &PrimDC[primsrc], local2global);
      plFx[primsrc] = tmpF.X;
      plFy[primsrc] = tmpF.Y;
      plFz[primsrc] = tmpF.Z;
    }  // for all primitives: basic device, mirror reflections and repetitions
  }    // pragma omp parallel

  double totPot = 0.0;
  Vector3D totF;
  totF.X = totF.Y = totF.Z = 0.0;
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    totPot += pPot[prim];
    totF.X += plFx[prim];
    totF.Y += plFy[prim];
    totF.Z += plFz[prim];
  }

  // This is done at the end of the function - before freeing memory
  *Potential = totPot;
  globalF->X = totF.X;
  globalF->Y = totF.Y;
  globalF->Z = totF.Z;

  (*Potential) += VSystemChargeZero;  // respect total system charge constraint

  if (dbgFn) {
    printf("Final values due to all primitives: ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n");	// refer, do not
    // uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", xfld, yfld, zfld,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  free_dvector(pPot, 1, NbPrimitives);
  free_dvector(plFx, 1, NbPrimitives);
  free_dvector(plFy, 1, NbPrimitives);
  free_dvector(plFz, 1, NbPrimitives);

  return (0);
}  // end of ElePFAtPoint

// Evaluate effects due to known charge distributions
// Since there is no intermediate function that interfaces PointKnChPF etc,
// division by MyFACTOR is necessary.
// CHECK OpenMP / GPU possibilities:
// Do parallelize before using these known charges - points or distributions
int KnChPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;
  double tmpPot;
  Point3D tmpPt;
  tmpPt.X = globalP->X;
  tmpPt.Y = globalP->Y;
  tmpPt.Z = globalP->Z;
  Vector3D tmpF;

  for (int point = 1; point <= NbPointsKnCh; ++point) {
    tmpPot = PointKnChPF((PointKnChArr + point - 1)->P, tmpPt, &tmpF);
    (*Potential) += (PointKnChArr + point - 1)->Assigned * tmpPot / MyFACTOR;
    globalF->X += (PointKnChArr + point - 1)->Assigned * tmpF.X / MyFACTOR;
    globalF->Y += (PointKnChArr + point - 1)->Assigned * tmpF.Y / MyFACTOR;
    globalF->Z += (PointKnChArr + point - 1)->Assigned * tmpF.Z / MyFACTOR;
  }

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    tmpPot = LineKnChPF((LineKnChArr + line - 1)->Start,
                        (LineKnChArr + line - 1)->Stop,
                        (LineKnChArr + line - 1)->Radius, tmpPt, &tmpF);
    (*Potential) += (LineKnChArr + line - 1)->Assigned * tmpPot / MyFACTOR;
    globalF->X += (LineKnChArr + line - 1)->Assigned * tmpF.X / MyFACTOR;
    globalF->Y += (LineKnChArr + line - 1)->Assigned * tmpF.Y / MyFACTOR;
    globalF->Z += (LineKnChArr + line - 1)->Assigned * tmpF.Z / MyFACTOR;
  }

  for (int area = 1; area <= NbAreasKnCh; ++area) {
    tmpPot = AreaKnChPF((AreaKnChArr + area - 1)->NbVertices,
                        ((AreaKnChArr + area - 1)->Vertex), tmpPt, &tmpF);
    (*Potential) += (AreaKnChArr + area - 1)->Assigned * tmpPot / MyFACTOR;
    globalF->X += (AreaKnChArr + area - 1)->Assigned * tmpF.X / MyFACTOR;
    globalF->Y += (AreaKnChArr + area - 1)->Assigned * tmpF.Y / MyFACTOR;
    globalF->Z += (AreaKnChArr + area - 1)->Assigned * tmpF.Z / MyFACTOR;
  }

  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    tmpPot = VolumeKnChPF((VolumeKnChArr + vol - 1)->NbVertices,
                          ((VolumeKnChArr + vol - 1)->Vertex), tmpPt, &tmpF);
    (*Potential) += (VolumeKnChArr + vol - 1)->Assigned * tmpPot / MyFACTOR;
    globalF->X += (VolumeKnChArr + vol - 1)->Assigned * tmpF.X / MyFACTOR;
    globalF->Y += (VolumeKnChArr + vol - 1)->Assigned * tmpF.Y / MyFACTOR;
    globalF->Z += (VolumeKnChArr + vol - 1)->Assigned * tmpF.Z / MyFACTOR;
  }

  if (dbgFn) {
    printf("Final values due to all known charges: ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n");	// refer, do not
    // uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", tmpPt.X, tmpPt.Y, tmpPt.Z,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  return 0;
}  // KnChPFAtPoint ends

// Compute voxelized data for export to Garfield++
int VoxelFPR(void) {
  int dbgFn = 0;
  int fstatus;

  int nbXCells;
  int nbYCells;
  int nbZCells;
  double startX;
  double startY;
  double startZ;
  double delX;
  double delY;
  double delZ;

  printf("\nPotential and field computation for voxelized data export\n");

  char VoxelFile[256];
  FILE *fVoxel;
  strcpy(VoxelFile, BCOutDir);
  strcat(VoxelFile, "/VoxelFPR.out");
  fVoxel = fopen(VoxelFile, "w");
  if (fVoxel == NULL) {
    neBEMMessage("VoxelFPR - VoxelFile");
    return -1;
  }
  fprintf(
      fVoxel,
      "# X(cm)\tY(cm)\tZ(cm)\tFX(V/cm)\tFY(V/cm)\tFZ(V/cm)\tPot(V)\tRegion\n");

  if (dbgFn) {
    printf("VoxelFPR.out created ...\n");
    fflush(stdout);
  }

  nbXCells = Voxel.NbXCells;
  nbYCells = Voxel.NbYCells;
  nbZCells = Voxel.NbZCells;
  startX = Voxel.Xmin;
  startY = Voxel.Ymin;
  startZ = Voxel.Zmin;
  delX = (Voxel.Xmax - Voxel.Xmin) / nbXCells;
  delY = (Voxel.Ymax - Voxel.Ymin) / nbYCells;
  delZ = (Voxel.Zmax - Voxel.Zmin) / nbZCells;

  int ivol;  // relates XYZ position to volume number
  double *VoxelFX, *VoxelFY, *VoxelFZ, *VoxelP;
  VoxelFX = dvector(0, nbZCells + 1);
  VoxelFY = dvector(0, nbZCells + 1);
  VoxelFZ = dvector(0, nbZCells + 1);
  VoxelP = dvector(0, nbZCells + 1);

  if (dbgFn) {
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    fflush(stdout);
  }

  for (int i = 1; i <= nbXCells + 1; ++i) {
    for (int j = 1; j <= nbYCells + 1; ++j) {
      if (dbgFn) {
        printf("VoxelFPR => i: %4d, j: %4d", i, j);
        fflush(stdout);
      }

      Point3D point;
#ifdef _OPENMP
      int nthreads = 1, tid = 0;
      #pragma omp parallel private(nthreads, tid)
#endif
      {
#ifdef _OPENMP
        if (dbgFn) {
          tid = omp_get_thread_num();
          if (tid == 0) {
            nthreads = omp_get_num_threads();
            printf("Starting voxel computation with %d threads\n", nthreads);
          }
        }
#endif
        int k;
        Vector3D field;
        double potential;
#ifdef _OPENMP
        #pragma omp for private(k, point, potential, field)
#endif
        for (k = 1; k <= nbZCells + 1; ++k) {
          potential = 0.0;
          field.X = field.Y = field.Z = 0.0;

          point.X = startX + (i - 1) * delX;  // all 3 components need to be
          point.Y = startY + (j - 1) * delY;  // evaluated after pragma omp
          point.Z = startZ + (k - 1) * delZ;

          if (dbgFn) {
            printf("i, j, k: %d, %d, %d\n", i, j, k);
            printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                   point.X / LengthScale, point.Y / LengthScale,
                   point.Z / LengthScale);
            fflush(stdout);
          }

          fstatus = PFAtPoint(&point, &potential, &field);
          if (fstatus != 0) {
            neBEMMessage("wrong PFAtPoint return value in VoxelFPR\n");
            // return -1;
          }

          if (dbgFn) {
            printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                   point.X / LengthScale, point.Y / LengthScale,
                   point.Z / LengthScale, field.X, field.Y, field.Z,
                   potential / LengthScale);
            fflush(stdout);
          }

          VoxelFX[k] = field.X;
          VoxelFY[k] = field.Y;
          VoxelFZ[k] = field.Z;
          VoxelP[k] = potential;
        }  // loop k
      }    // pragma omp parallel

      for (int k = 1; k <= nbZCells + 1; ++k)  // file output
      {
        point.X = startX + (i - 1) * delX;
        point.Y = startY + (j - 1) * delY;
        point.Z = startZ + (k - 1) * delZ;

        ivol = neBEMVolumePoint(point.X, point.Y, point.Z);
        /*
        volMaterial[ivol];	// region linked to material
        neBEMVolumeDescription(ivol, &vshp, &vmat, &veps, &vpot, &vq, &vtype);
if(dbgFn)
{
printf("volref: %d\n", ivol);
printf("shape: %d,  material: %d\n", volShape[ivol], volMaterial[ivol]);
printf("eps: %d,  pot: %d\n", volEpsilon[ivol], volPotential[ivol]);
printf("q: %d,  type: %d\n", volCharge[ivol], volBoundaryType[ivol]);
printf("shape: %d,  material: %d\n", vshp, vmat);
printf("eps: %d,  pot: %d\n", veps, vpot);
printf("q: %d,  type: %d\n", vq, vtype);
}
        */

        fprintf(fVoxel,
                "%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%4d\n",
                100.0 * point.X / LengthScale, 100.0 * point.Y / LengthScale,
                100.0 * point.Z / LengthScale, VoxelFX[k] / 100.0,
                VoxelFY[k] / 100.0, VoxelFZ[k] / 100.0, VoxelP[k] / LengthScale,
                ivol + 1);
      }
      fflush(fVoxel);  // file output over

      // printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    }  // loop j
  }    // loop i

  fclose(fVoxel);

  free_dvector(VoxelFX, 0, nbZCells + 1);
  free_dvector(VoxelFY, 0, nbZCells + 1);
  free_dvector(VoxelFZ, 0, nbZCells + 1);
  free_dvector(VoxelP, 0, nbZCells + 1);

  return 0;
}  // end of VoxelFPR

// Compute 3dMap data for export to Garfield++.
// The 3dMap data for weighting field can be generated using the same function.
// After creation, the exported file has to be properly
// named and then imported by ComponentNeBem3dMap.
int MapFPR(void) {
  int dbgFn = 0;
  int fstatus;

  int nbXCells;
  int nbYCells;
  int nbZCells;
  double startX;
  double startY;
  double startZ;
  double delX;
  double delY;
  double delZ;

  printf("\nPotential and field computation for 3dMap data export\n");

  char MapInfoFile[256];
  FILE *fMapInfo;
  strcpy(MapInfoFile, BCOutDir);
  strcat(MapInfoFile, "/MapInfo.out");
  fMapInfo = fopen(MapInfoFile, "w");
  if (fMapInfo == NULL) {
    neBEMMessage("MapFPR - MapInfoFile");
    return -1;
  }
  if (dbgFn) {
    printf("MapInfoFile.out created ...\n");
    fflush(stdout);
  }

  // in certain versions, we may have only the version number in the header and
  // nothing more. In that case, it is unsafe to assume that OptMap or
  // OptStaggerMap will be at all present in the output file. This decision
  // may need to be taken immediately after reading the MapVersion value.
  fprintf(fMapInfo, "%s\n", MapVersion);
  fprintf(fMapInfo, "%d\n", OptMap);
  fprintf(fMapInfo, "%d\n", OptStaggerMap);
  fprintf(fMapInfo, "%d\n", Map.NbXCells + 1);
  fprintf(fMapInfo, "%d\n", Map.NbYCells + 1);
  fprintf(fMapInfo, "%d\n", Map.NbZCells + 1);
  fprintf(fMapInfo, "%le %le\n", Map.Xmin * 100.0, Map.Xmax * 100.0);
  fprintf(fMapInfo, "%le %le\n", Map.Ymin * 100.0, Map.Ymax * 100.0);
  fprintf(fMapInfo, "%le %le\n", Map.Zmin * 100.0, Map.Zmax * 100.0);
  fprintf(fMapInfo, "%le\n", Map.XStagger * 100.0);
  fprintf(fMapInfo, "%le\n", Map.YStagger * 100.0);
  fprintf(fMapInfo, "%le\n", Map.ZStagger * 100.0);
  fprintf(fMapInfo, "MapFPR.out\n");
  // if(OptStaggerMap) fprintf(fMapInfo, "StgrMapFPR.out\n"); not being read
  fclose(fMapInfo);

  char MapFile[256];
  FILE *fMap;
  strcpy(MapFile, BCOutDir);
  strcat(MapFile, "/MapFPR.out");
  fMap = fopen(MapFile, "w");
  if (fMap == NULL) {
    neBEMMessage("MapFPR - MapFile");
    return -1;
  }
  if (dbgFn) {
    printf("MapFPR.out created ...\n");
    fflush(stdout);
  }

  fprintf(
      fMap,
      "# X(cm)\tY(cm)\tZ(cm)\tFX(V/cm)\tFY(V/cm)\tFZ(V/cm)\tPot(V)\tRegion\n");

  nbXCells = Map.NbXCells;
  nbYCells = Map.NbYCells;
  nbZCells = Map.NbZCells;
  startX = Map.Xmin;
  startY = Map.Ymin;
  startZ = Map.Zmin;
  delX = (Map.Xmax - Map.Xmin) / nbXCells;
  delY = (Map.Ymax - Map.Ymin) / nbYCells;
  delZ = (Map.Zmax - Map.Zmin) / nbZCells;

  int ivol;  // relates XYZ position to volume number using Garfield subroutine
  double *MapFX, *MapFY, *MapFZ, *MapP;
  MapFX = dvector(0, nbZCells + 1);
  MapFY = dvector(0, nbZCells + 1);
  MapFZ = dvector(0, nbZCells + 1);
  MapP = dvector(0, nbZCells + 1);

  if (dbgFn) {
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    fflush(stdout);
  }

  for (int i = 1; i <= nbXCells + 1; ++i) {
    for (int j = 1; j <= nbYCells + 1; ++j) {
      if (dbgFn) {
        printf("MapFPR => i: %4d, j: %4d", i, j);
        fflush(stdout);
      }

      Point3D point;
#ifdef _OPENMP
      int nthreads = 1, tid = 0;
      #pragma omp parallel private(nthreads, tid)
#endif
      {
#ifdef _OPENMP
        if (dbgFn) {
          tid = omp_get_thread_num();
          if (tid == 0) {
            nthreads = omp_get_num_threads();
            printf("Starting voxel computation with %d threads\n", nthreads);
          }
        }
#endif
        int k;
        Vector3D field;
        double potential;
#ifdef _OPENMP
        #pragma omp for private(k, point, potential, field)
#endif
        for (k = 1; k <= nbZCells + 1; ++k) {
          point.X = startX + (i - 1) * delX;  // all 3 components need to be
          point.Y = startY + (j - 1) * delY;  // evaluated after pragma omp
          point.Z = startZ + (k - 1) * delZ;

          if (dbgFn) {
            printf("i, j, k: %d, %d, %d\n", i, j, k);
            printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                   point.X / LengthScale, point.Y / LengthScale,
                   point.Z / LengthScale);
            fflush(stdout);
          }

          if (OptReadFastPF) {
            fstatus = FastPFAtPoint(&point, &potential, &field);
            if (fstatus != 0) {
              neBEMMessage("wrong FastPFAtPoint return value in MapFPR\n");
              // return -1;
            }
          } else {
            neBEMMessage(
                "Suggestion: Use of FastVol can expedite generation of Map.\n");
            fstatus = PFAtPoint(&point, &potential, &field);
            if (fstatus != 0) {
              neBEMMessage("wrong PFAtPoint return value in MapFPR\n");
              // return -1;
            }
          }

          if (dbgFn) {
            printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                   point.X / LengthScale, point.Y / LengthScale,
                   point.Z / LengthScale, field.X, field.Y, field.Z,
                   potential / LengthScale);
            fflush(stdout);
          }

          MapFX[k] = field.X;
          MapFY[k] = field.Y;
          MapFZ[k] = field.Z;
          MapP[k] = potential;
        }  // loop k
      }    // pragma omp parallel

      for (int k = 1; k <= nbZCells + 1; ++k)  // file output
      {
        point.X = startX + (i - 1) * delX;
        point.Y = startY + (j - 1) * delY;
        point.Z = startZ + (k - 1) * delZ;

        ivol = neBEMVolumePoint(point.X, point.Y, point.Z);
        /*
        volMaterial[ivol];	// region linked to material
        neBEMVolumeDescription(ivol, &vshp, &vmat, &veps, &vpot, &vq, &vtype);
if(dbgFn)
{
printf("volref: %d\n", ivol);
printf("shape: %d,  material: %d\n", volShape[ivol], volMaterial[ivol]);
printf("eps: %d,  pot: %d\n", volEpsilon[ivol], volPotential[ivol]);
printf("q: %d,  type: %d\n", volCharge[ivol], volBoundaryType[ivol]);
printf("shape: %d,  material: %d\n", vshp, vmat);
printf("eps: %d,  pot: %d\n", veps, vpot);
printf("q: %d,  type: %d\n", vq, vtype);
}
        */

        fprintf(fMap, "%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%4d\n",
                100.0 * point.X / LengthScale, 100.0 * point.Y / LengthScale,
                100.0 * point.Z / LengthScale, MapFX[k] / 100.0,
                MapFY[k] / 100.0, MapFZ[k] / 100.0, MapP[k] / LengthScale,
                ivol + 1);
      }
      fflush(fMap);  // file output over

      // printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    }  // loop j
  }    // loop i

  fclose(fMap);

  free_dvector(MapFX, 0, nbZCells + 1);
  free_dvector(MapFY, 0, nbZCells + 1);
  free_dvector(MapFZ, 0, nbZCells + 1);
  free_dvector(MapP, 0, nbZCells + 1);

  // If staggered map
  if (OptStaggerMap) {
    char MapFile[256];
    FILE *fMap;
    strcpy(MapFile, BCOutDir);
    strcat(MapFile, "/StgrMapFPR.out");
    fMap = fopen(MapFile, "w");
    if (fMap == NULL) {
      neBEMMessage("StgrMapFPR - Staggered MapFile");
      return -1;
    }
    if (dbgFn) {
      printf("StgrMapFPR.out created ...\n");
      fflush(stdout);
    }

    fprintf(
        fMap,
        "# "
        "X(cm)\tY(cm)\tZ(cm)\tFX(V/cm)\tFY(V/cm)\tFZ(V/cm)\tPot(V)\tRegion\n");

    double LX = (Map.Xmax - Map.Xmin);
    Map.Xmin = Map.Xmax;  // very static stagger where X-shift is one map long
    Map.Xmax = Map.Xmin + LX;
    double LY = (Map.Ymax - Map.Ymin);
    Map.Ymin = Map.Ymin + Map.YStagger;
    Map.Ymax = Map.Ymin + LY;
    nbXCells = Map.NbXCells;
    nbYCells = Map.NbYCells;
    nbZCells = Map.NbZCells;
    startX = Map.Xmin;
    startY =
        Map.Ymin + Map.YStagger;  // and y-shift is of the presecribed amount
    startZ = Map.Zmin;
    delX = (Map.Xmax - Map.Xmin) / nbXCells;
    delY = (Map.Ymax - Map.Ymin) / nbYCells;
    delZ = (Map.Zmax - Map.Zmin) / nbZCells;

    int ivol;  // relates XYZ position to volume number using Garfield
               // subroutine
    double *MapFX, *MapFY, *MapFZ, *MapP;
    MapFX = dvector(0, nbZCells + 1);
    MapFY = dvector(0, nbZCells + 1);
    MapFZ = dvector(0, nbZCells + 1);
    MapP = dvector(0, nbZCells + 1);

    if (dbgFn) {
      printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
             nbZCells);
      printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
      printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
      fflush(stdout);
    }

    for (int i = 1; i <= nbXCells + 1; ++i) {
      for (int j = 1; j <= nbYCells + 1; ++j) {
        if (dbgFn) {
          printf("StgrMapFPR => i: %4d, j: %4d", i, j);
          fflush(stdout);
        }

        Point3D point;
#ifdef _OPENMP
        int nthreads = 1, tid = 0;
        #pragma omp parallel private(nthreads, tid)
#endif
        {
#ifdef _OPENMP
          if (dbgFn) {
            tid = omp_get_thread_num();
            if (tid == 0) {
              nthreads = omp_get_num_threads();
              printf("Starting voxel computation with %d threads\n", nthreads);
            }
          }  // if dbgFn
#endif
          int k;
          Vector3D field;
          double potential;
#ifdef _OPENMP
          #pragma omp for private(k, point, potential, field)
#endif
          for (k = 1; k <= nbZCells + 1; ++k) {
            point.X = startX + (i - 1) * delX;  // all 3 components need to be
            point.Y = startY + (j - 1) * delY;  // evaluated after pragma omp
            point.Z = startZ + (k - 1) * delZ;

            if (dbgFn) {
              printf("i, j, k: %d, %d, %d\n", i, j, k);
              printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale);
              fflush(stdout);
            }

            if (OptReadFastPF) {
              fstatus = FastPFAtPoint(&point, &potential, &field);
              if (fstatus != 0) {
                neBEMMessage("wrong FastPFAtPoint return value in MapFPR\n");
                // return -1;
              }
            } else {
              neBEMMessage(
                  "Suggestion: Use of FastVol can expedite generation of "
                  "Map.\n");
              fstatus = PFAtPoint(&point, &potential, &field);
              if (fstatus != 0) {
                neBEMMessage("wrong PFAtPoint return value in MapFPR\n");
                // return -1;
              }
            }

            if (dbgFn) {
              printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale, field.X, field.Y, field.Z,
                     potential / LengthScale);
              fflush(stdout);
            }

            MapFX[k] = field.X;
            MapFY[k] = field.Y;
            MapFZ[k] = field.Z;
            MapP[k] = potential;
          }  // loop k
        }    // pragma omp parallel

        for (int k = 1; k <= nbZCells + 1; ++k)  // file output
        {
          point.X = startX + (i - 1) * delX;
          point.Y = startY + (j - 1) * delY;
          point.Z = startZ + (k - 1) * delZ;

          ivol = neBEMVolumePoint(point.X, point.Y, point.Z);
          /*
          volMaterial[ivol];	// region linked to material
          neBEMVolumeDescription(ivol, &vshp, &vmat, &veps, &vpot, &vq, &vtype);
  if(dbgFn)
  {
  printf("volref: %d\n", ivol);
  printf("shape: %d,  material: %d\n", volShape[ivol], volMaterial[ivol]);
  printf("eps: %d,  pot: %d\n", volEpsilon[ivol], volPotential[ivol]);
  printf("q: %d,  type: %d\n", volCharge[ivol], volBoundaryType[ivol]);
  printf("shape: %d,  material: %d\n", vshp, vmat);
  printf("eps: %d,  pot: %d\n", veps, vpot);
  printf("q: %d,  type: %d\n", vq, vtype);
  }
          */

          fprintf(
              fMap, "%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%4d\n",
              100.0 * point.X / LengthScale, 100.0 * point.Y / LengthScale,
              100.0 * point.Z / LengthScale, MapFX[k] / 100.0, MapFY[k] / 100.0,
              MapFZ[k] / 100.0, MapP[k] / LengthScale, ivol + 1);
        }              // for k <= nbZCells
        fflush(fMap);  // file output over

        // printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
      }  // loop j
    }    // loop i

    fclose(fMap);

    free_dvector(MapFX, 0, nbZCells + 1);
    free_dvector(MapFY, 0, nbZCells + 1);
    free_dvector(MapFZ, 0, nbZCells + 1);
    free_dvector(MapP, 0, nbZCells + 1);
  }  // If staggered map

  return 0;
}  // end of MapFPR

// Compute potential and field in a mesh within the Fast Volume
int FastVolPF(void) {

  // The following may be necessary only during the first time step / iteration
  // At present, FastVolElePF() considers both element and KnCh effects
  // and create one combined fast volume.
  int fstatus = FastVolElePF();
  if (fstatus) {
    printf(
        "Problem in FastVolElePF being called from FastVolPF ... returning\n");
    return -1;
  }

  /*
  // The following is likely to change throughout the computation and necessary
  // at all time steps. However, in order to achieve computational economy, it
  // may be prudent to carry out the following only after several time steps.
  if (OptKnCh) {
    fstatus = FastVolKnChPF();
    if (fstatus) {
      printf("Problem in FastVolKnChPF being called from FastVolPF... returning\n"); 
      return -1;
    }
  }
  */

  return 0;
}  // FastVolPF ends

// Compute potential and field in a mesh within the Fast Volume
// Possible pitfall: evaluation of n-skips
// As the name implies, this function uses the ElePFAtPoint function only.
// As a result, KnChPFAtPoint is not included in the resulting values.
// The effects due to known charge distributions can be included separately
// as if they are perturbations on the background system.
int FastVolElePF(void) {
  int dbgFn = 0;
  int fstatus;

  int nbXCells;
  int nbYCells;
  int nbZCells;
  double startX;
  double startY;
  double startZ;
  double delX;
  double delY;
  double delZ;

  printf("\nPotential and field computation within basic fast volume\n");
  int bskip = 0, iskip = 0, jskip = 0, kskip = 0;

  // calculate n-skips based on NbPtSkip
  if (NbPtSkip) {
    int volptcnt = 0, endskip = 0;

    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      nbXCells = BlkNbXCells[block];
      nbYCells = BlkNbYCells[block];
      nbZCells = BlkNbZCells[block];
      for (int i = 1; i <= nbXCells + 1; ++i) {
        for (int j = 1; j <= nbYCells + 1; ++j) {
          for (int k = 1; k <= nbZCells + 1; ++k) {
            ++volptcnt;

            if (volptcnt == NbPtSkip) {
              bskip = block - 1;
              iskip = i - 1;
              jskip = j - 1;
              kskip = k;
              endskip = 1;
            }

            if (endskip) break;
          }
          if (endskip) break;
        }
        if (endskip) break;
      }
      if (endskip) break;
    }
    if (dbgFn) {
      printf(
          "Basic fast volume => bskip, iskip, jskip, kskip: %d, %d, %d, %d\n",
          bskip, iskip, jskip, kskip);
    }
  }  // NbPtSkip

  char FastVolPFFile[256];
  FILE *fFastVolPF;
  strcpy(FastVolPFFile, BCOutDir);
  strcat(FastVolPFFile, "/FastVolPF.out");
  fFastVolPF = fopen(FastVolPFFile, "w");
  if (fFastVolPF == NULL) {
    neBEMMessage("FastVolPF - FastVolPFFile");
    return -1;
  }
  fprintf(fFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

  if (dbgFn) {
    printf("FastVolPF.out created ...\n");
    fflush(stdout);
  }

  for (int block = 1 + bskip; block <= FastVol.NbBlocks; ++block) {
    nbXCells = BlkNbXCells[block];
    nbYCells = BlkNbYCells[block];
    nbZCells = BlkNbZCells[block];
    startX = FastVol.CrnrX;
    startY = FastVol.CrnrY;
    startZ = BlkCrnrZ[block];
    delX = FastVol.LX / nbXCells;
    delY = FastVol.LY / nbYCells;
    delZ = BlkLZ[block] / nbZCells;
    printf(
        "NbBlocks: %d, block: %d, nbXCells: %d, nbYCells: %d, nbZCells: %d\n",
        FastVol.NbBlocks, block, nbXCells, nbYCells, nbZCells);

    if (dbgFn) {
      printf("block: %d\n", block);
      printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
             nbZCells);
      printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
      printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
      printf("bskip, iskip, jskip, kskip: %d, %d, %d, %d\n", bskip, iskip,
             jskip, kskip);
      fflush(stdout);
    }
    // total number of points in a given block
    // int blktotpt = (nbXCells + 1) * (nbYCells + 1) * (nbZCells + 1);
    for (int i = 1 + iskip; i <= nbXCells + 1; ++i) {
      for (int j = 1 + jskip; j <= nbYCells + 1; ++j) {
        printf("Fast volume => block: %3d, i: %4d, j: %4d", block, i, j);
        fflush(stdout);

        Point3D point;
#ifdef _OPENMP
        int nthreads = 1, tid = 0;
        #pragma omp parallel private(nthreads, tid)
#endif
        {
#ifdef _OPENMP
          if (dbgFn) {
            tid = omp_get_thread_num();
            if (tid == 0) {
              nthreads = omp_get_num_threads();
              printf("Starting fast volume computation with %d threads\n",
                     nthreads);
            }
          }
#endif
          int k;
          int omitFlag;
          double potential;
          Vector3D field;
#ifdef _OPENMP
          #pragma omp for private(k, point, omitFlag, potential, field)
#endif
          for (k = 1 + kskip; k <= nbZCells + 1; ++k) {
            potential = 0.0;
            field.X = field.Y = field.Z = 0.0;

            point.X = startX + (i - 1) * delX;
            point.Y = startY + (j - 1) * delY;
            point.Z = startZ + (k - 1) * delZ;

            // Check whether the point falls within a volume that should be
            // ignored
            omitFlag = 0;
            for (int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
              if ((point.X > OmitVolCrnrX[omit]) &&
                  (point.X < OmitVolCrnrX[omit] + OmitVolLX[omit]) &&
                  (point.Y > OmitVolCrnrY[omit]) &&
                  (point.Y < OmitVolCrnrY[omit] + OmitVolLY[omit]) &&
                  (point.Z > OmitVolCrnrZ[omit]) &&
                  (point.Z < OmitVolCrnrZ[omit] + OmitVolLZ[omit])) {
                omitFlag = 1;
                break;
              }
            }  // loop over omitted volumes

            if (dbgFn) {
              printf("block, i, j, k: %d, %d, %d, %d\n", block, i, j, k);
              printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale);
              printf("omitFlag: %d\n", omitFlag);
              fflush(stdout);
            }

            if (omitFlag) {
              potential = field.X = field.Y = field.Z = 0.0;
            } else {
              // fstatus = ElePFAtPoint(&point, &potential, &field);
              fstatus =
                  PFAtPoint(&point, &potential, &field);  // both ele and KnCh
              if (fstatus != 0) {
                neBEMMessage(
                    "wrong ElePFAtPoint return value in FastVolElePF.\n");
                // return -1;
              }
            }
            if (dbgFn) {
              printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale, potential / LengthScale, field.X,
                     field.Y, field.Z);
              fflush(stdout);
            }

            FastPot[block][i][j][k] = potential;
            FastFX[block][i][j][k] = field.X;
            FastFY[block][i][j][k] = field.Y;
            FastFZ[block][i][j][k] = field.Z;
          }  // loop k
        }    // pragma omp parallel

        for (int k = 1 + kskip; k <= nbZCells + 1; ++k)  // file output
        {
          point.X = startX + (i - 1) * delX;
          point.Y = startY + (j - 1) * delY;
          point.Z = startZ + (k - 1) * delZ;

          fprintf(fFastVolPF,
                  "%4d\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                  block, point.X / LengthScale, point.Y / LengthScale,
                  point.Z / LengthScale, FastPot[block][i][j][k],
                  FastFX[block][i][j][k], FastFY[block][i][j][k],
                  FastFZ[block][i][j][k]);
        }
        fflush(fFastVolPF);  // file output over

        printf(
            "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
            "\b\b\b\b\b\b\b\b\b\b");
      }  // loop j
    }    // loop i
  }      // loop block

  fclose(fFastVolPF);

  if (OptStaggerFastVol) {
    printf("Potential and field computation within staggered fast volume\n");

    bskip = iskip = jskip = kskip = 0;

    // calculate n-skips based on NbStgPtSkip
    if (NbStgPtSkip) {
      int volptcnt = 0, endskip = 0;

      for (int block = 1; block <= FastVol.NbBlocks; ++block) {
        nbXCells = BlkNbXCells[block];
        nbYCells = BlkNbYCells[block];
        nbZCells = BlkNbZCells[block];
        for (int i = 1; i <= nbXCells + 1; ++i) {
          for (int j = 1; j <= nbYCells + 1; ++j) {
            for (int k = 1; k <= nbZCells + 1; ++k) {
              ++volptcnt;

              if (volptcnt == NbStgPtSkip) {
                bskip = block - 1;
                iskip = i - 1;
                jskip = j - 1;
                kskip = k;
                endskip = 1;
              }

              if (endskip) break;
            }
            if (endskip) break;
          }
          if (endskip) break;
        }
        if (endskip) break;
      }
      if (dbgFn) {
        printf(
            "Staggered volume => bskip, iskip, jskip, kskip: %d, %d, %d, %d\n",
            bskip, iskip, jskip, kskip);
      }
    }  // NbStgPtSkip

    char FastStgVolPFFile[256];
    FILE *fFastStgVolPF;
    strcpy(FastStgVolPFFile, BCOutDir);
    strcat(FastStgVolPFFile, "/FastStgVolPF.out");
    fFastStgVolPF = fopen(FastStgVolPFFile, "w");
    if (fFastStgVolPF == NULL) {
      neBEMMessage("FastVolPF - FastStgVolPFFile");
      return -1;
    }
    fprintf(fFastStgVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

    if (dbgFn) {
      printf("FastStgVolPF.out created ...\n");
      fflush(stdout);
    }

    for (int block = 1 + bskip; block <= FastVol.NbBlocks; ++block) {
      nbXCells = BlkNbXCells[block];
      nbYCells = BlkNbYCells[block];
      nbZCells = BlkNbZCells[block];
      startX = FastVol.CrnrX + FastVol.LX;
      startY = FastVol.CrnrY + FastVol.YStagger;
      startZ = BlkCrnrZ[block];
      delX = FastVol.LX / nbXCells;
      delY = FastVol.LY / nbYCells;
      delZ = BlkLZ[block] / nbZCells;
      printf(
          "NbBlocks: %d, block: %d, nbXCells: %d, nbYCells: %d, nbZCells: %d\n",
          FastVol.NbBlocks, block, nbXCells, nbYCells, nbZCells);

      if (dbgFn) {
        printf("block: %d\n", block);
        printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
               nbZCells);
        printf("startX, startY, startZ: %le, %le, %le\n", startX, startY,
               startZ);
        printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
        printf("bskip, iskip, jskip, kskip: %d, %d, %d, %d\n", bskip, iskip,
               jskip, kskip);
        fflush(stdout);
      }

      // int blktotpt = (nbXCells + 1) * (nbYCells + 1) * (nbZCells + 1);
      for (int i = 1 + iskip; i <= nbXCells + 1; ++i) {
        for (int j = 1 + jskip; j <= nbYCells + 1; ++j) {
          printf("Fast volume => block: %3d, i: %4d, j: %4d", block, i, j);
          fflush(stdout);

          Point3D point;
#ifdef _OPENMP
          int nthreads = 1, tid = 0;
          #pragma omp parallel private(nthreads, tid)
#endif
          {
#ifdef _OPENMP
            if (dbgFn) {
              tid = omp_get_thread_num();
              if (tid == 0) {
                nthreads = omp_get_num_threads();
                printf(
                    "Starting staggered fast volume computation with %d "
                    "threads\n",
                    nthreads);
              }
            }
#endif
            int k;
            int omitFlag;
            double potential;
            Vector3D field;
#ifdef _OPENMP
            #pragma omp for private(k, point, omitFlag, potential, field)
#endif
            for (k = 1 + kskip; k <= nbZCells + 1; ++k) {
              potential = 0.0;
              field.X = field.Y = field.Z = 0.0;

              point.X = startX + (i - 1) * delX;
              point.Y = startY + (j - 1) * delY;
              point.Z = startZ + (k - 1) * delZ;

              // Check whether point falls within a volume that should be
              // ignored
              omitFlag = 0;
              for (int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
                if ((point.X > OmitVolCrnrX[omit] + FastVol.LX) &&
                    (point.X <
                     OmitVolCrnrX[omit] + OmitVolLX[omit] + FastVol.LX) &&
                    (point.Y > OmitVolCrnrY[omit] + FastVol.YStagger) &&
                    (point.Y <
                     OmitVolCrnrY[omit] + OmitVolLY[omit] + FastVol.YStagger) &&
                    (point.Z > OmitVolCrnrZ[omit]) &&
                    (point.Z < OmitVolCrnrZ[omit] + OmitVolLZ[omit])) {
                  omitFlag = 1;
                  break;
                }
              }  // loop over omitted volumes

              if (dbgFn) {
                printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                       point.X / LengthScale, point.Y / LengthScale,
                       point.Z / LengthScale);
                printf("omitFlag: %d\n", omitFlag);
                fflush(stdout);
              }

              if (omitFlag) {
                potential = field.X = field.Y = field.Z = 0.0;
              } else {
                // fstatus = ElePFAtPoint(&point, &potential, &field);
                fstatus =
                    PFAtPoint(&point, &potential, &field);  // both ele & KnCh
                if (fstatus != 0) {
                  neBEMMessage(
                      "wrong PFAtPoint return value in FastVolElePF.\n");
                  // return -1;
                }
              }
              if (dbgFn) {
                printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                       point.X / LengthScale, point.Y / LengthScale,
                       point.Z / LengthScale, potential / LengthScale, field.X,
                       field.Y, field.Z);
                fflush(stdout);
              }

              FastStgPot[block][i][j][k] = potential;
              FastStgFX[block][i][j][k] = field.X;
              FastStgFY[block][i][j][k] = field.Y;
              FastStgFZ[block][i][j][k] = field.Z;
            }  // loop k
          }    // pragma omp

          for (int k = 1 + kskip; k <= nbZCells + 1; ++k)  // file output
          {
            point.X = startX + (i - 1) * delX;
            point.Y = startY + (j - 1) * delY;
            point.Z = startZ + (k - 1) * delZ;

            fprintf(fFastStgVolPF,
                    "%4d\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                    block, point.X / LengthScale, point.Y / LengthScale,
                    point.Z / LengthScale, FastPot[block][i][j][k],
                    FastFX[block][i][j][k], FastFY[block][i][j][k],
                    FastFZ[block][i][j][k]);
          }
          fflush(fFastStgVolPF);  // file output over

          printf(
              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
              "\b\b\b\b\b\b\b\b\b\b\b");
        }  // loop j
      }    // loop i
    }      // loop block

    fclose(fFastStgVolPF);
  }  // if OptStaggerFastVol

  return 0;
}  // FastVolElePF ends

/* This at present is unnecessary since FastVolElePF considers effects due to
 * elements, as well as known charges.
// Compute potential and field in a mesh within the Fast Volume
// Possible pitfall: evaluation of n-skips
// As the name implies, this function does NOT use the ElePFAtPoint function.
// Only KnChPFAtPoint is included in the resulting values.
// These effects due to known charge distributions are expected to be
// perturbations on the background system.
int FastVolKnChPF(void)
{
int dbgFn = 0;
int fstatus;

int nbXCells;
int nbYCells;
int nbZCells;
double startX;
double startY;
double startZ;
double delX;
double delY;
double delZ;

printf("\nPotential and field computation within basic fast volume\n");
int blktotpt;	// total number of points in a given block
int bskip = 0, iskip = 0, jskip = 0, kskip = 0;

// calculate n-skips based on NbPtSkip
if(NbPtSkip)
        {
        int volptcnt = 0, endskip = 0;

        for(int block = 1; block <= FastVol.NbBlocks; ++block)
                {
                nbXCells = BlkNbXCells[block];
                nbYCells = BlkNbYCells[block];
                nbZCells = BlkNbZCells[block];
                for(int i = 1; i <= nbXCells+1; ++i)
                        {
                for(int j = 1; j <= nbYCells+1; ++j)
                {
                for(int k = 1; k <= nbZCells+1; ++k)
                {
                                        ++volptcnt;

                                        if(volptcnt == NbPtSkip)
                                                {
                                                bskip = block-1; iskip = i-1;
jskip = j-1; kskip = k; endskip = 1;
                                                }

                                        if(endskip) break;
                                        }
                                if(endskip) break;
                                }
                        if(endskip) break;
                        }
                if(endskip) break;
                }
        if(dbgFn)
                {
                printf("Basic fast volume => bskip, iskip, jskip, kskip: %d, %d,
%d, %d\n", bskip, iskip, jskip, kskip);
                }
        }	// NbPtSkip

char FastVolKnChPFFile[256];
FILE *fFastVolKnChPF;
strcpy(FastVolKnChPFFile, BCOutDir);
strcat(FastVolKnChPFFile, "/FastVolKnChPF.out");
fFastVolKnChPF = fopen(FastVolKnChPFFile, "w");
if(fFastVolKnChPF == NULL)
        {
        neBEMMessage("FastVolKnChPF - FastVolKnChPFFile");
        return -1;
        }
fprintf(fFastVolKnChPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

if(dbgFn)
        {
        printf("FastVolKnChPF.out created ...\n"); fflush(stdout);
        }

for(int block = 1+bskip; block <= FastVol.NbBlocks; ++block)
        {
        nbXCells = BlkNbXCells[block];
        nbYCells = BlkNbYCells[block];
        nbZCells = BlkNbZCells[block];
        startX = FastVol.CrnrX;
        startY = FastVol.CrnrY;
        startZ = BlkCrnrZ[block];
        delX = FastVol.LX / nbXCells;
        delY = FastVol.LY / nbYCells;
        delZ = BlkLZ[block] / nbZCells;
        printf("NbBlocks: %d, block: %d, nbXCells: %d, nbYCells: %d, nbZCells:
%d\n", FastVol.NbBlocks, block, nbXCells, nbYCells, nbZCells);

        if(dbgFn)
                {
                printf("block: %d\n", block);
                printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n",
                                                nbXCells, nbYCells, nbZCells);
                printf("startX, startY, startZ: %le, %le, %le\n", startX,
startY, startZ); printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
                printf("bskip, iskip, jskip, kskip: %d, %d, %d, %d\n",
                                                bskip, iskip, jskip, kskip);
                fflush(stdout);
                }

        blktotpt = (nbXCells+1)*(nbYCells+1)*(nbZCells+1);
        for(int i = 1+iskip; i <= nbXCells+1; ++i)
                {
        for(int j = 1+jskip; j <= nbYCells+1; ++j)
        {

                        printf("Fast volume => block: %3d, i: %4d, j: %4d",
block, i, j); fflush(stdout);

                        Point3D point;
                        int nthreads = 1, tid = 0;
                        #pragma omp parallel private(nthreads, tid)
                        {
                        if(dbgFn)
                                {
                                tid = omp_get_thread_num();
                                if (tid == 0)
                                {
                                nthreads = omp_get_num_threads();
                                printf("Starting fast volume computation with %d
threads\n", nthreads);
                                }
                                }

                        int k;
                        int omitFlag;
                        double potential;
                        Vector3D field;
                        #pragma omp for private(k, point, omitFlag, potential,
field) for(k = 1+kskip; k <= nbZCells+1; ++k)
        {
                        point.X = startX + (i-1)*delX;
                point.Y = startY + (j-1)*delY;
        point.Z = startZ + (k-1)*delZ;

                                // Check whether the point falls within a volume
that should be ignored omitFlag = 0; for(int omit = 1; omit <=
FastVol.NbOmitVols; ++omit)
                                {
                                        if((point.X > OmitVolCrnrX[omit])
                                                        && (point.X <
OmitVolCrnrX[omit]+OmitVolLX[omit])
                                                        && (point.Y >
OmitVolCrnrY[omit])
                                                        && (point.Y <
OmitVolCrnrY[omit]+OmitVolLY[omit])
                                                        && (point.Z >
OmitVolCrnrZ[omit])
                                                        && (point.Z <
OmitVolCrnrZ[omit]+OmitVolLZ[omit]))
                                                {
                                                omitFlag = 1;
                                                break;
                                                }
                                        }	// loop over omitted volumes

                                if(dbgFn)
                                        {
                                        printf("block, i, j, k: %d, %d, %d,
%d\n", block, i, j, k); printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                        point.X/LengthScale, point.Y/LengthScale,
point.Z/LengthScale); printf("omitFlag: %d\n", omitFlag); fflush(stdout);
                                        }

                                if(omitFlag)
                                        {
                                        potential = field.X = field.Y = field.Z
= 0.0;
                                        }
                                else
                                        {
                fstatus = KnChPFAtPoint(&point, &potential, &field);
                if(fstatus != 0)
                                                {
                                                neBEMMessage("wrong
KnChPFAtPoint return value in FastVolKnChPF.\n");
                                                // return -1;
                                                }
                                        }
                                if(dbgFn)
                                        {
                printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                        point.X/LengthScale, point.Y/LengthScale,
point.Z/LengthScale, potential/LengthScale, field.X, field.Y, field.Z);
                                        fflush(stdout);
                                        }

                                FastPot[block][i][j][k] = potential;
                                FastFX[block][i][j][k] = field.X;
                                FastFY[block][i][j][k] = field.Y;
                                FastFZ[block][i][j][k] = field.Z;
                }	// loop k
                        }	// pragma omp parallel

                        for(int k = 1+kskip; k <= nbZCells+1; ++k)	// file
output
                        {
                        point.X = startX + (i-1)*delX;
                point.Y = startY + (j-1)*delY;
                        point.Z = startZ + (k-1)*delZ;

                fprintf(fFastVolKnChPF,
                        "%4d\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                        block,
                                                                point.X/LengthScale,
point.Y/LengthScale, point.Z/LengthScale,
                        FastPot[block][i][j][k],FastFX[block][i][j][k],
                                                                FastFY[block][i][j][k],
FastFZ[block][i][j][k]);
                                }
                        fflush(fFastVolKnChPF);	// file output over

                        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                }	// loop j
                }	// loop i
        }	// loop block

fclose(fFastVolKnChPF);

if(OptStaggerFastVol)
        {
        printf("Potential and field computation within staggered fast
volume\n");

        bskip = iskip = jskip = kskip = 0;

        // calculate n-skips based on NbStgPtSkip
        if(NbStgPtSkip)
                {
                int volptcnt = 0, endskip = 0;

                for(int block = 1; block <= FastVol.NbBlocks; ++block)
                        {
                        nbXCells = BlkNbXCells[block];
                        nbYCells = BlkNbYCells[block];
                        nbZCells = BlkNbZCells[block];
                        for(int i = 1; i <= nbXCells+1; ++i)
                                {
                        for(int j = 1; j <= nbYCells+1; ++j)
                        {
                        for(int k = 1; k <= nbZCells+1; ++k)
                        {
                                                ++volptcnt;

                                                if(volptcnt == NbStgPtSkip)
                                                        {
                                                        bskip = block-1; iskip =
i-1; jskip = j-1; kskip = k; endskip = 1;
                                                        }

                                                if(endskip) break;
                                                }
                                        if(endskip) break;
                                        }
                                if(endskip) break;
                                }
                        if(endskip) break;
                        }
                if(dbgFn)
                        {
                        printf("Staggered volume => bskip, iskip, jskip, kskip:
%d, %d, %d, %d\n", bskip, iskip, jskip, kskip);
                        }
                }	// NbStgPtSkip

        char FastStgVolKnChPFFile[256];
        FILE *fFastStgVolKnChPF;
        strcpy(FastStgVolKnChPFFile, BCOutDir);
        strcat(FastStgVolKnChPFFile, "/FastStgVolKnChPF.out");
        fFastStgVolKnChPF = fopen(FastStgVolKnChPFFile, "w");
        if(fFastStgVolKnChPF == NULL)
                {
                neBEMMessage("FastVolKnChPF - FastStgVolKnChPFFile");
                return -1;
                }
        fprintf(fFastStgVolKnChPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

        if(dbgFn)
                {
                printf("FastStgVolKnChPF.out created ...\n"); fflush(stdout);
                }

        for(int block = 1+bskip; block <= FastVol.NbBlocks; ++block)
                {
                nbXCells = BlkNbXCells[block];
                nbYCells = BlkNbYCells[block];
                nbZCells = BlkNbZCells[block];
                startX = FastVol.CrnrX + FastVol.LX;
                startY = FastVol.CrnrY + FastVol.YStagger;
                startZ = BlkCrnrZ[block];
                delX = FastVol.LX / nbXCells;
                delY = FastVol.LY / nbYCells;
                delZ = BlkLZ[block] / nbZCells;
                printf(
                                        "NbBlocks: %d, block: %d, nbXCells: %d,
nbYCells: %d, nbZCells: %d\n", FastVol.NbBlocks, block, nbXCells, nbYCells,
nbZCells);

                if(dbgFn)
                        {
                        printf("block: %d\n", block);
                        printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n",
                                                        nbXCells, nbYCells,
nbZCells); printf("startX, startY, startZ: %le, %le, %le\n", startX, startY,
startZ); printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
                        printf("bskip, iskip, jskip, kskip: %d, %d, %d, %d\n",
                                                        bskip, iskip, jskip,
kskip); fflush(stdout);
                        }

                blktotpt = (nbXCells+1)*(nbYCells+1)*(nbZCells+1);
                for(int i = 1+iskip; i <= nbXCells+1; ++i)
                        {
                for(int j = 1+jskip; j <= nbYCells+1; ++j)
                {

                                printf("Fast volume => block: %3d, i: %4d, j:
%4d", block, i, j); fflush(stdout);

                                Point3D point;
                                int nthreads, tid;
                                #pragma omp parallel private(nthreads, tid)
                                {
                                if(dbgFn)
                                        {
                                        tid = omp_get_thread_num();
                                        if (tid == 0)
                                        {
                                        nthreads = omp_get_num_threads();
                                        printf(
                                                                "Starting
staggered fast volume computation with %d threads\n", nthreads);
                                        }
                                }

                                int k;
                                int omitFlag;
                                double potential;
                                Vector3D field;
                                #pragma omp for private(k, point, omitFlag,
potential, field) for(k = 1+kskip; k <= nbZCells+1; ++k)
                {
                                point.X = startX + (i-1)*delX;
                        point.Y = startY + (j-1)*delY;
                point.Z = startZ + (k-1)*delZ;

                                        // Check whether point falls within a
volume that should be ignored omitFlag = 0; for(int omit = 1; omit <=
FastVol.NbOmitVols; ++omit)
                                        {
                                                if((point.X >
OmitVolCrnrX[omit]+FastVol.LX)
                                                                && (point.X <
OmitVolCrnrX[omit]+OmitVolLX[omit]+FastVol.LX)
                                                                && (point.Y >
OmitVolCrnrY[omit]+FastVol.YStagger)
                                                                && (point.Y <
OmitVolCrnrY[omit]+OmitVolLY[omit]+FastVol.YStagger)
                                                                && (point.Z >
OmitVolCrnrZ[omit])
                                                                && (point.Z <
OmitVolCrnrZ[omit]+OmitVolLZ[omit]))
                                                        {
                                                        omitFlag = 1;
                                                        break;
                                                        }
                                                }	// loop over omitted
volumes

                                        if(dbgFn)
                                                {
                        printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                        point.X/LengthScale, point.Y/LengthScale,
point.Z/LengthScale); printf("omitFlag: %d\n", omitFlag); fflush(stdout);
                                                }

                                        if(omitFlag)
                                                {
                                                potential = field.X = field.Y =
field.Z = 0.0;
                                                }
                                        else
                                                {
                        fstatus = KnChPFAtPoint(&point, &potential, &field);
                        if(fstatus != 0)
                                                        {
                                                        neBEMMessage("wrong
KnChPFAtPoint return value in FastVolKnChPF.\n");
                                                        // return -1;
                                                        }
                                                }
                                        if(dbgFn)
                                                {
                        printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                        point.X/LengthScale, point.Y/LengthScale,
point.Z/LengthScale, potential/LengthScale, field.X, field.Y, field.Z);
                                                fflush(stdout);
                                                }

                                        FastStgPot[block][i][j][k] = potential;
                                        FastStgFX[block][i][j][k] = field.X;
                                        FastStgFY[block][i][j][k] = field.Y;
                                        FastStgFZ[block][i][j][k] = field.Z;
                }	// loop k
                                }	// pragma omp

                                for(int k = 1+kskip; k <= nbZCells+1; ++k)
// file output
                                {
                                point.X = startX + (i-1)*delX;
                        point.Y = startY + (j-1)*delY;
                                point.Z = startZ + (k-1)*delZ;

                        fprintf(fFastStgVolKnChPF,
                                "%4d\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                                block,
                                                                        point.X/LengthScale,
point.Y/LengthScale, point.Z/LengthScale,
                                FastPot[block][i][j][k],FastFX[block][i][j][k],
                                                                        FastFY[block][i][j][k],
FastFZ[block][i][j][k]);
                                        }
                                fflush(fFastStgVolKnChPF);	// file output
over

                                printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                }	// loop j
                }	// loop i
                }	// loop block

        fclose(fFastStgVolKnChPF);
        }	// if OptStaggerFastVol

return 0;
}	// FastVolKnChPF ends
*/

// Gives three components of the total Potential and flux in the global
// coordinate system due to all the elements using the results stored in
// the FAST volume mesh.
int FastPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;
  double Xpt = globalP->X;
  double Ypt = globalP->Y;
  double Zpt = globalP->Z;
  double RptVolLX = FastVol.LX;
  double RptVolLY = FastVol.LY;
  double RptVolLZ = FastVol.LZ;
  double CornerX = FastVol.CrnrX;
  double CornerY = FastVol.CrnrY;
  double CornerZ = FastVol.CrnrZ;
  double TriLin(double xd, double yd, double zd, double c000, double c100,
                double c010, double c001, double c110, double c101, double c011,
                double c111);

  // First of all, check how the point in question should be treated ...

  // Check whether the point falls within a volume that is not regarded as
  // FastVol
  for (int ignore = 1; ignore <= FastVol.NbIgnoreVols; ++ignore) {
    if ((Xpt >= (IgnoreVolCrnrX[ignore])) &&
        (Xpt <= (IgnoreVolCrnrX[ignore] + IgnoreVolLX[ignore])) &&
        (Ypt >= (IgnoreVolCrnrY[ignore])) &&
        (Ypt <= (IgnoreVolCrnrY[ignore] + IgnoreVolLY[ignore])) &&
        (Zpt >= (IgnoreVolCrnrZ[ignore])) &&
        (Zpt <= (IgnoreVolCrnrZ[ignore] + IgnoreVolLZ[ignore]))) {
      if (dbgFn)
        neBEMMessage("In FastPFAtPoint: point in an ignored volume!\n");

      int fstatus = PFAtPoint(globalP, Potential, globalF);
      if (fstatus != 0) {
        neBEMMessage("wrong PFAtPoint return value in FastVolPF.\n");
        return -1;
      } else
        return 0;
    }
  }  // loop over ignored volumes

  // If not ignored, the point qualifies for FastVol evaluation ...

  // for a staggered fast volume, the volume repeated in X is larger
  if (OptStaggerFastVol) {
    RptVolLX += FastVol.LX;
  }
  if (dbgFn) {
    printf("\nin FastPFAtPoint\n");
    printf("x, y, z: %g, %g, %g\n", Xpt, Ypt, Zpt);
    printf("RptVolLX, RptVolLY, RptVolLZ: %g, %g, %g\n", RptVolLX, RptVolLY,
           RptVolLZ);
    printf("CornerX, CornerY, CornerZ: %g, %g, %g\n", CornerX, CornerY,
           CornerZ);
    printf("Nb of blocks: %d\n", FastVol.NbBlocks);
    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      printf("NbOfXCells: %d\n", BlkNbXCells[block]);
      printf("NbOfYCells: %d\n", BlkNbYCells[block]);
      printf("NbOfZCells: %d\n", BlkNbZCells[block]);
      printf("LZ: %le\n", BlkLZ[block]);
      printf("CornerZ: %le\n", BlkCrnrZ[block]);
    }
  }

  // Find equivalent position inside the basic / staggered volume.
  // real distance from volume corner
  double dx = Xpt - CornerX;
  double dy = Ypt - CornerY;
  double dz = Zpt - CornerZ;
  if (dbgFn)
    printf("real dx, dy, dz from volume corner: %g, %g, %g\n", dx, dy, dz);

  int NbFastVolX = (int)(dx / RptVolLX);
  if (dx < 0.0) --NbFastVolX;
  int NbFastVolY = (int)(dy / RptVolLY);
  if (dy < 0.0) --NbFastVolY;
  int NbFastVolZ = (int)(dz / RptVolLZ);
  if (dz < 0.0) --NbFastVolZ;
  if (dbgFn)
    printf("Volumes in x, y, z: %d, %d, %d\n", NbFastVolX, NbFastVolY,
           NbFastVolZ);

  // equivalent distances from fast volume corner
  dx -= NbFastVolX * RptVolLX;
  dy -= NbFastVolY * RptVolLY;
  dz -= NbFastVolZ * RptVolLZ;
  // The following conditions should never happen - generate an error message
  if (dx < 0.0) {
    dx = 0.0;
    neBEMMessage("equiv dx < 0.0 - not correct!\n");
  }
  if (dy < 0.0) {
    dy = 0.0;
    neBEMMessage("equiv dy < 0.0 - not correct!\n");
  }
  if (dz < 0.0) {
    dz = 0.0;
    neBEMMessage("equiv dz < 0.0 - not correct!\n");
  }
  if (dx > RptVolLX) {
    dx = RptVolLX;
    neBEMMessage("equiv dx > RptVolLX - not correct!\n");
  }
  if (dy > RptVolLY) {
    dy = RptVolLY;
    neBEMMessage("equiv dy > RptVolLY - not correct!\n");
  }
  if (dz > RptVolLZ) {
    dz = RptVolLZ;
    neBEMMessage("equiv dz > RptVolLZ - not correct!\n");
  }
  if (dbgFn)
    printf("equivalent dist from corner - dx, dy, dz: %g, %g, %g\n", dx, dy,
           dz);

  // Take care of possible trouble-makers
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= FastVol.LX) && (FastVol.LX - dx) < MINDIST)
    dx = FastVol.LX - MINDIST;
  // else if((dx > FastVol.LX) && (fabs(FastVol.LX-dx) < MINDIST))
  else if ((dx > FastVol.LX) && (dx - FastVol.LX) < MINDIST)  // more inttuitive
    dx = FastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted - dx, dy, dz: %g, %g, %g\n", dx, dy, dz);

  // If volume is staggered, we have a few more things to do before finalizing
  // the values of equivalent distance
  // sector identification
  // _................__________________
  // |    .     .     |    Sector 3    |
  // |     .   .      |                |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 2    |    .     .     |
  // |----------------|     .   .      |
  // |                |       .        |
  // |       .        |                |
  // |     .   .      |                |
  // |    .     .     |----------------|
  // |     .   .      |    Sector 4    |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 1    |    .     .     |
  // |----------------|................|

  int sector = 1;  // kept outside `if' since this is necessary further below
  if (OptStaggerFastVol) {
    if ((dx >= 0.0) && (dx <= FastVol.LX) && (dy >= 0.0) &&
        (dy <= FastVol.LY)) {
      // point lies in sector 1, everything remains unchanged
      sector = 1;
    } else if ((dx >= 0.0) && (dx <= FastVol.LX) && (dy > FastVol.LY) &&
               (dy <= FastVol.LY + FastVol.YStagger)) {
      // point lies in sector 2, move basic volume one step up
      sector = 2;
      ++NbFastVolY;
      CornerY += FastVol.LY;  // repeat length in Y is LY
      dy -= FastVol.LY;
    } else if ((dx > FastVol.LX) && (dx <= 2.0 * FastVol.LX) &&
               (dy >= FastVol.YStagger) &&
               (dy <= FastVol.LY + FastVol.YStagger)) {
      // point lies in sector 3, pt in staggered vol, change corner coords
      sector = 3;
      CornerX += FastVol.LX;
      CornerY += FastVol.YStagger;
      dx -= FastVol.LX;
      dy -= FastVol.YStagger;
    } else if ((dx > FastVol.LX) && (dx <= 2.0 * FastVol.LX) && (dy >= 0.0) &&
               (dy < FastVol.YStagger)) {
      // point lies in sector 4, move basic volume one step down and consider
      // staggered fast volume
      sector = 4;
      --NbFastVolY;
      CornerX += FastVol.LX;  // in the staggered part of the repeated volume
      CornerY -= (FastVol.LY - FastVol.YStagger);
      dx -= FastVol.LX;
      dy += (FastVol.LY - FastVol.YStagger);
    } else {
      neBEMMessage("FastPFAtPoint: point in none of the sectors!\n");
    }
    if (dbgFn) printf("stagger modified dx, dy, dz: %g, %g, %g\n", dx, dy, dz);
  }

  // Take care of possible trouble-makers - once more
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= FastVol.LX) && (FastVol.LX - dx) < MINDIST)
    dx = FastVol.LX - MINDIST;
  // else if((dx > FastVol.LX) && (fabs(FastVol.LX-dx) < MINDIST))
  else if ((dx > FastVol.LX) && (dx - FastVol.LX) < MINDIST)  // more intuitive
    dx = FastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted for staggered: %g, %g, %g\n", dx, dy, dz);

  /*
  // Check whether the point falls within a volume that is omitted
  for(int omit = 1; omit <= FastVol.NbOmitVols; ++omit)
    {
          if((dx >= (OmitVolCrnrX[omit]-FastVol.CrnrX))
                          && (dx <=
  (OmitVolCrnrX[omit]+OmitVolLX[omit]-FastVol.CrnrX))
                          && (dy >= (OmitVolCrnrY[omit]-FastVol.CrnrY))
                          && (dy <=
  (OmitVolCrnrY[omit]+OmitVolLY[omit]-FastVol.CrnrY))
                          && (dz >= (OmitVolCrnrZ[omit]-FastVol.CrnrZ))
                          && (dz <=
  (OmitVolCrnrZ[omit]+OmitVolLZ[omit]-FastVol.CrnrZ)))
                  {
                  neBEMMessage("In FastPFAtPoint: point in an omitted
  volume!\n"); *Potential = 0.0; globalF->X = 0.0; globalF->Y = 0.0; globalF->Z
  = 0.0;
                  }
          }	// loop over omitted volumes
  */

  // Find the block in which the point lies
  /*
  int thisBlock = 1;
  if(FastVol.NbBlocks > 1)
          {
          for(int block = 1; block <= FastVol.NbBlocks; ++block)
                  {
                  if(dbgFn)
                          {
                          printf("dz,(BlkCrnrZ-CornerZ),(BlkCrnrZ+BlkLZ-CornerZ):
  %lg, %lg, %lg\n", dz, (BlkCrnrZ[block]-CornerZ),
                                                          (BlkCrnrZ[block]+BlkLZ[block]-CornerZ));
                          }
                  if((dz >= (BlkCrnrZ[block]-CornerZ))
                                  && (dz <=
  (BlkCrnrZ[block]+BlkLZ[block]-CornerZ)))
                          {
                          thisBlock = block;
                          break;
                          }
                  }
          }	// if NbBlocks > 1
  */

  int thisBlock = 0;
  for (int block = 1; block <= FastVol.NbBlocks; ++block) {
    double blkBtmZ = BlkCrnrZ[block] - CornerZ;  // since CornerZ has been
    double blkTopZ = blkBtmZ + BlkLZ[block];     // subtracted from dz already
    if (dbgFn) {
      printf("block, dz, blkBtmZ, blkTopZ: %d, %lg, %lg, %lg\n", block, dz,
             blkBtmZ, blkTopZ);
    }

    // take care of difficult situations
    if ((dz <= blkBtmZ) && ((blkBtmZ - dz) < MINDIST)) dz = blkBtmZ - MINDIST;
    if ((dz >= blkBtmZ) && ((dz - blkBtmZ) < MINDIST)) dz = blkBtmZ + MINDIST;
    if ((dz <= blkTopZ) && ((blkTopZ - dz) < MINDIST)) dz = blkTopZ - MINDIST;
    if ((dz >= blkTopZ) && ((dz - blkTopZ) < MINDIST)) dz = blkTopZ + MINDIST;

    if ((dz >= blkBtmZ) && (dz <= blkTopZ)) {
      thisBlock = block;
      break;
    }
  }
  if (!thisBlock) {
    neBEMMessage("FastPFAtPoint: point in none of the blocks!\n");
  }

  int nbXCells = BlkNbXCells[thisBlock];
  int nbYCells = BlkNbYCells[thisBlock];
  int nbZCells = BlkNbZCells[thisBlock];
  double delX = FastVol.LX / nbXCells;
  double delY = FastVol.LY / nbYCells;
  double delZ = BlkLZ[thisBlock] / nbZCells;
  dz -= (BlkCrnrZ[thisBlock] - CornerZ);  // distance from the block corner

  if (dbgFn) {
    printf("thisBlock: %d\n", thisBlock);
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("BlkCrnrZ: %lg\n", BlkCrnrZ[thisBlock]);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    printf("dz: %lg\n", dz);
    fflush(stdout);
  }

  // Find cell in block of basic / staggered volume within which the point lies
  int celli = (int)(dx / delX) + 1;  // Find cell in which the point lies
  if (celli < 1) {
    celli = 1;
    dx = 0.5 * delX;
    neBEMMessage("FastPFAtPoint - celli < 1\n");
  }
  if (celli > nbXCells) {
    celli = nbXCells;
    dx = FastVol.LX - 0.5 * delX;
    neBEMMessage("FastPFAtPoint - celli > nbXCells\n");
  }
  int cellj = (int)(dy / delY) + 1;
  if (cellj < 1) {
    cellj = 1;
    dy = 0.5 * delY;
    neBEMMessage("FastPFAtPoint - cellj < 1\n");
  }
  if (cellj > nbYCells) {
    cellj = nbYCells;
    dy = FastVol.LY - 0.5 * delY;
    neBEMMessage("FastPFAtPoint - cellj > nbYCells\n");
  }
  int cellk = (int)(dz / delZ) + 1;
  if (cellk < 1) {
    cellk = 1;
    dz = 0.5 * delX;
    neBEMMessage("FastPFAtPoint - cellk < 1\n");
  }
  if (cellk > nbZCells) {
    cellk = nbZCells;
    dz = FastVol.LZ - 0.5 * delZ;
    neBEMMessage("FastPFAtPoint - cellk > nbZCells\n");
  }
  if (dbgFn) printf("Cells in x, y, z: %d, %d, %d\n", celli, cellj, cellk);

  // Interpolate potential and field at the point using the corner values of
  // of the cell and, if necessary, of the neighbouring cells
  // These gradients can also be calculated while computing the potential and
  // field at the cells and stored in memory, provided enough memory is
  // available

  // distances from cell corner
  double dxcellcrnr = dx - (double)(celli - 1) * delX;
  double dycellcrnr = dy - (double)(cellj - 1) * delY;
  double dzcellcrnr = dz - (double)(cellk - 1) * delZ;
  if (dbgFn)
    printf("cell crnr dx, dy, dz: %g, %g, %g\n", dxcellcrnr, dycellcrnr,
           dzcellcrnr);

  // normalized distances
  double xd = dxcellcrnr / delX;  // xd = (x-x0)/(x1-x0)
  double yd = dycellcrnr / delY;  // etc
  double zd = dzcellcrnr / delZ;
  if (xd <= 0.0) xd = 0.0;
  if (yd <= 0.0) yd = 0.0;
  if (zd <= 0.0) zd = 0.0;
  if (xd >= 1.0) xd = 1.0;
  if (yd >= 1.0) yd = 1.0;
  if (zd >= 1.0) zd = 1.0;

  // corner values of potential and field
  double P000 = FastPot[thisBlock][celli][cellj][cellk];  // lowest corner
  double FX000 = FastFX[thisBlock][celli][cellj][cellk];
  double FY000 = FastFY[thisBlock][celli][cellj][cellk];
  double FZ000 = FastFZ[thisBlock][celli][cellj][cellk];
  double P100 = FastPot[thisBlock][celli + 1][cellj][cellk];
  double FX100 = FastFX[thisBlock][celli + 1][cellj][cellk];
  double FY100 = FastFY[thisBlock][celli + 1][cellj][cellk];
  double FZ100 = FastFZ[thisBlock][celli + 1][cellj][cellk];
  double P010 = FastPot[thisBlock][celli][cellj + 1][cellk];
  double FX010 = FastFX[thisBlock][celli][cellj + 1][cellk];
  double FY010 = FastFY[thisBlock][celli][cellj + 1][cellk];
  double FZ010 = FastFZ[thisBlock][celli][cellj + 1][cellk];
  double P001 = FastPot[thisBlock][celli][cellj][cellk + 1];
  double FX001 = FastFX[thisBlock][celli][cellj][cellk + 1];
  double FY001 = FastFY[thisBlock][celli][cellj][cellk + 1];
  double FZ001 = FastFZ[thisBlock][celli][cellj][cellk + 1];
  double P110 = FastPot[thisBlock][celli + 1][cellj + 1][cellk];
  double FX110 = FastFX[thisBlock][celli + 1][cellj + 1][cellk];
  double FY110 = FastFY[thisBlock][celli + 1][cellj + 1][cellk];
  double FZ110 = FastFZ[thisBlock][celli + 1][cellj + 1][cellk];
  double P101 = FastPot[thisBlock][celli + 1][cellj][cellk + 1];
  double FX101 = FastFX[thisBlock][celli + 1][cellj][cellk + 1];
  double FY101 = FastFY[thisBlock][celli + 1][cellj][cellk + 1];
  double FZ101 = FastFZ[thisBlock][celli + 1][cellj][cellk + 1];
  double P011 = FastPot[thisBlock][celli][cellj + 1][cellk + 1];
  double FX011 = FastFX[thisBlock][celli][cellj + 1][cellk + 1];
  double FY011 = FastFY[thisBlock][celli][cellj + 1][cellk + 1];
  double FZ011 = FastFZ[thisBlock][celli][cellj + 1][cellk + 1];
  double P111 = FastPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FX111 = FastFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FY111 = FastFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FZ111 = FastFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
  if (OptStaggerFastVol) {
    if (sector == 1) {  // nothing to be done
    }
    if (sector == 2) {  // volume shifted up but point not in the staggered part
    }
    if (sector == 3) {  // staggered volume
      P000 = FastStgPot[thisBlock][celli][cellj][cellk];
      FX000 = FastStgFX[thisBlock][celli][cellj][cellk];
      FY000 = FastStgFY[thisBlock][celli][cellj][cellk];
      FZ000 = FastStgFZ[thisBlock][celli][cellj][cellk];
      P100 = FastStgPot[thisBlock][celli + 1][cellj][cellk];
      FX100 = FastStgFX[thisBlock][celli + 1][cellj][cellk];
      FY100 = FastStgFY[thisBlock][celli + 1][cellj][cellk];
      FZ100 = FastStgFZ[thisBlock][celli + 1][cellj][cellk];
      P010 = FastStgPot[thisBlock][celli][cellj + 1][cellk];
      FX010 = FastStgFX[thisBlock][celli][cellj + 1][cellk];
      FY010 = FastStgFY[thisBlock][celli][cellj + 1][cellk];
      FZ010 = FastStgFZ[thisBlock][celli][cellj + 1][cellk];
      P001 = FastStgPot[thisBlock][celli][cellj][cellk + 1];
      FX001 = FastStgFX[thisBlock][celli][cellj][cellk + 1];
      FY001 = FastStgFY[thisBlock][celli][cellj][cellk + 1];
      FZ001 = FastStgFZ[thisBlock][celli][cellj][cellk + 1];
      P110 = FastStgPot[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = FastStgFX[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = FastStgFY[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = FastStgFZ[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = FastStgPot[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = FastStgFX[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = FastStgFY[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = FastStgFZ[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = FastStgPot[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = FastStgFX[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = FastStgFY[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = FastStgFZ[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = FastStgPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = FastStgFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = FastStgFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = FastStgFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
    if (sector == 4) {  // volume shifted down and point in the staggered part
      P000 = FastStgPot[thisBlock][celli][cellj][cellk];
      FX000 = FastStgFX[thisBlock][celli][cellj][cellk];
      FY000 = FastStgFY[thisBlock][celli][cellj][cellk];
      FZ000 = FastStgFZ[thisBlock][celli][cellj][cellk];
      P100 = FastStgPot[thisBlock][celli + 1][cellj][cellk];
      FX100 = FastStgFX[thisBlock][celli + 1][cellj][cellk];
      FY100 = FastStgFY[thisBlock][celli + 1][cellj][cellk];
      FZ100 = FastStgFZ[thisBlock][celli + 1][cellj][cellk];
      P010 = FastStgPot[thisBlock][celli][cellj + 1][cellk];
      FX010 = FastStgFX[thisBlock][celli][cellj + 1][cellk];
      FY010 = FastStgFY[thisBlock][celli][cellj + 1][cellk];
      FZ010 = FastStgFZ[thisBlock][celli][cellj + 1][cellk];
      P001 = FastStgPot[thisBlock][celli][cellj][cellk + 1];
      FX001 = FastStgFX[thisBlock][celli][cellj][cellk + 1];
      FY001 = FastStgFY[thisBlock][celli][cellj][cellk + 1];
      FZ001 = FastStgFZ[thisBlock][celli][cellj][cellk + 1];
      P110 = FastStgPot[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = FastStgFX[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = FastStgFY[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = FastStgFZ[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = FastStgPot[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = FastStgFX[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = FastStgFY[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = FastStgFZ[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = FastStgPot[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = FastStgFX[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = FastStgFY[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = FastStgFZ[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = FastStgPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = FastStgFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = FastStgFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = FastStgFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
  }  // if OptStaggerFastVol

  double intP =
      TriLin(xd, yd, zd, P000, P100, P010, P001, P110, P101, P011, P111);
  double intFX = TriLin(xd, yd, zd, FX000, FX100, FX010, FX001, FX110, FX101,
                        FX011, FX111);
  double intFY = TriLin(xd, yd, zd, FY000, FY100, FY010, FY001, FY110, FY101,
                        FY011, FY111);
  double intFZ = TriLin(xd, yd, zd, FZ000, FZ100, FZ010, FZ001, FZ110, FZ101,
                        FZ011, FZ111);

  *Potential = intP;
  globalF->X = intFX;
  globalF->Y = intFY;
  globalF->Z = intFZ;

  if (dbgFn) {
    printf("Cell corner values:\n");
    printf("Potential: %g, %g, %g, %g\n", P000, P100, P010, P001);
    printf("Potential: %g, %g, %g, %g\n", P110, P101, P011, P111);
    printf("FastFX: %g, %g, %g, %g\n", FX000, FX100, FX010, FX001);
    printf("FastFX: %g, %g, %g, %g\n", FX110, FX101, FX011, FX111);
    printf("FastFY: %g, %g, %g, %g\n", FY000, FY100, FY010, FY001);
    printf("FastFY: %g, %g, %g, %g\n", FY110, FY101, FY011, FY111);
    printf("FastFZ: %g, %g, %g, %g\n", FZ000, FZ100, FZ010, FZ001);
    printf("FastFZ: %g, %g, %g, %g\n", FZ110, FZ101, FZ011, FZ111);
    printf("Pot, FX, FY, FZ: %g, %g, %g, %g\n", *Potential, globalF->X,
           globalF->Y, globalF->Z);
  }

  if (dbgFn) {
    printf("out FastPFAtPoint\n");
    fflush(stdout);
  }

  return 0;
}  // FastPFAtPoint ends

// There could be a function
int FastElePFAtPoint(Point3D* /*globalP*/, double* /*Potential*/, 
                     Vector3D* /*globalF*/) {
  return 0;
}  // FastElePFAtPoint ends

// Gives three components of the total Potential and flux in the global
// coordinate system due to all the known charges using the results stored in
// the FAST volume KnCh mesh.
int FastKnChPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;
  double Xpt = globalP->X;
  double Ypt = globalP->Y;
  double Zpt = globalP->Z;
  double RptVolLX = FastVol.LX;
  double RptVolLY = FastVol.LY;
  double RptVolLZ = FastVol.LZ;
  double CornerX = FastVol.CrnrX;
  double CornerY = FastVol.CrnrY;
  double CornerZ = FastVol.CrnrZ;
  double TriLin(double xd, double yd, double zd, double c000, double c100,
                double c010, double c001, double c110, double c101, double c011,
                double c111);

  // First of all, check how the point in question should be treated ...

  // Check whether the point falls within a volume that is not regarded as
  // FastVol
  for (int ignore = 1; ignore <= FastVol.NbIgnoreVols; ++ignore) {
    if ((Xpt >= (IgnoreVolCrnrX[ignore])) &&
        (Xpt <= (IgnoreVolCrnrX[ignore] + IgnoreVolLX[ignore])) &&
        (Ypt >= (IgnoreVolCrnrY[ignore])) &&
        (Ypt <= (IgnoreVolCrnrY[ignore] + IgnoreVolLY[ignore])) &&
        (Zpt >= (IgnoreVolCrnrZ[ignore])) &&
        (Zpt <= (IgnoreVolCrnrZ[ignore] + IgnoreVolLZ[ignore]))) {
      if (dbgFn)
        neBEMMessage("In FastKnChPFAtPoint: point in an ignored volume!\n");

      int fstatus = KnChPFAtPoint(globalP, Potential, globalF);
      if (fstatus != 0) {
        neBEMMessage("wrong KnChPFAtPoint return value in FastVolKnChPF.\n");
        return -1;
      } else
        return 0;
    }
  }  // loop over ignored volumes

  // If not ignored, the point qualifies for FastVol evaluation ...

  // for a staggered fast volume, the volume repeated in X is larger
  if (OptStaggerFastVol) {
    RptVolLX += FastVol.LX;
  }
  if (dbgFn) {
    printf("\nin FastKnChPFAtPoint\n");
    printf("x, y, z: %g, %g, %g\n", Xpt, Ypt, Zpt);
    printf("RptVolLX, RptVolLY, RptVolLZ: %g, %g, %g\n", RptVolLX, RptVolLY,
           RptVolLZ);
    printf("CornerX, CornerY, CornerZ: %g, %g, %g\n", CornerX, CornerY,
           CornerZ);
    printf("Nb of blocks: %d\n", FastVol.NbBlocks);
    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      printf("NbOfXCells: %d\n", BlkNbXCells[block]);
      printf("NbOfYCells: %d\n", BlkNbYCells[block]);
      printf("NbOfZCells: %d\n", BlkNbZCells[block]);
      printf("LZ: %le\n", BlkLZ[block]);
      printf("CornerZ: %le\n", BlkCrnrZ[block]);
    }
  }

  // Find equivalent position inside the basic / staggered volume.
  // real distance from volume corner
  double dx = Xpt - CornerX;
  double dy = Ypt - CornerY;
  double dz = Zpt - CornerZ;
  if (dbgFn)
    printf("real dx, dy, dz from volume corner: %g, %g, %g\n", dx, dy, dz);

  int NbFastVolX = (int)(dx / RptVolLX);
  if (dx < 0.0) --NbFastVolX;
  int NbFastVolY = (int)(dy / RptVolLY);
  if (dy < 0.0) --NbFastVolY;
  int NbFastVolZ = (int)(dz / RptVolLZ);
  if (dz < 0.0) --NbFastVolZ;
  if (dbgFn)
    printf("Volumes in x, y, z: %d, %d, %d\n", NbFastVolX, NbFastVolY,
           NbFastVolZ);

  // equivalent distances from fast volume corner
  dx -= NbFastVolX * RptVolLX;
  dy -= NbFastVolY * RptVolLY;
  dz -= NbFastVolZ * RptVolLZ;
  // The following conditions should never happen - generate an error message
  if (dx < 0.0) {
    dx = 0.0;
    neBEMMessage("equiv dx < 0.0 - not correct!\n");
  }
  if (dy < 0.0) {
    dy = 0.0;
    neBEMMessage("equiv dy < 0.0 - not correct!\n");
  }
  if (dz < 0.0) {
    dz = 0.0;
    neBEMMessage("equiv dz < 0.0 - not correct!\n");
  }
  if (dx > RptVolLX) {
    dx = RptVolLX;
    neBEMMessage("equiv dx > RptVolLX - not correct!\n");
  }
  if (dy > RptVolLY) {
    dy = RptVolLY;
    neBEMMessage("equiv dy > RptVolLY - not correct!\n");
  }
  if (dz > RptVolLZ) {
    dz = RptVolLZ;
    neBEMMessage("equiv dz > RptVolLZ - not correct!\n");
  }
  if (dbgFn)
    printf("equivalent dist from corner - dx, dy, dz: %g, %g, %g\n", dx, dy,
           dz);

  // Take care of possible trouble-makers
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= FastVol.LX) && (FastVol.LX - dx) < MINDIST)
    dx = FastVol.LX - MINDIST;
  // else if((dx > FastVol.LX) && (fabs(FastVol.LX-dx) < MINDIST))
  else if ((dx > FastVol.LX) && (dx - FastVol.LX) < MINDIST)  // more inttuitive
    dx = FastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted - dx, dy, dz: %g, %g, %g\n", dx, dy, dz);

  // If volume is staggered, we have a few more things to do before finalizing
  // the values of equivalent distance
  // sector identification
  // _................__________________
  // |    .     .     |    Sector 3    |
  // |     .   .      |                |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 2    |    .     .     |
  // |----------------|     .   .      |
  // |                |       .        |
  // |       .        |                |
  // |     .   .      |                |
  // |    .     .     |----------------|
  // |     .   .      |    Sector 4    |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 1    |    .     .     |
  // |----------------|................|

  int sector = 1;  // kept outside `if' since this is necessary further below
  if (OptStaggerFastVol) {
    if ((dx >= 0.0) && (dx <= FastVol.LX) && (dy >= 0.0) &&
        (dy <= FastVol.LY)) {
      // point lies in sector 1, everything remains unchanged
      sector = 1;
    } else if ((dx >= 0.0) && (dx <= FastVol.LX) && (dy > FastVol.LY) &&
               (dy <= FastVol.LY + FastVol.YStagger)) {
      // point lies in sector 2, move basic volume one step up
      sector = 2;
      ++NbFastVolY;
      CornerY += FastVol.LY;  // repeat length in Y is LY
      dy -= FastVol.LY;
    } else if ((dx > FastVol.LX) && (dx <= 2.0 * FastVol.LX) &&
               (dy >= FastVol.YStagger) &&
               (dy <= FastVol.LY + FastVol.YStagger)) {
      // point lies in sector 3, pt in staggered vol, change corner coords
      sector = 3;
      CornerX += FastVol.LX;
      CornerY += FastVol.YStagger;
      dx -= FastVol.LX;
      dy -= FastVol.YStagger;
    } else if ((dx > FastVol.LX) && (dx <= 2.0 * FastVol.LX) && (dy >= 0.0) &&
               (dy < FastVol.YStagger)) {
      // point lies in sector 4, move basic volume one step down and consider
      // staggered fast volume
      sector = 4;
      --NbFastVolY;
      CornerX += FastVol.LX;  // in the staggered part of the repeated volume
      CornerY -= (FastVol.LY - FastVol.YStagger);
      dx -= FastVol.LX;
      dy += (FastVol.LY - FastVol.YStagger);
    } else {
      neBEMMessage("FastKnChPFAtPoint: point in none of the sectors!\n");
    }
    if (dbgFn) printf("stagger modified dx, dy, dz: %g, %g, %g\n", dx, dy, dz);
  }

  // Take care of possible trouble-makers - once more
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= FastVol.LX) && (FastVol.LX - dx) < MINDIST)
    dx = FastVol.LX - MINDIST;
  // else if((dx > FastVol.LX) && (fabs(FastVol.LX-dx) < MINDIST))
  else if ((dx > FastVol.LX) && (dx - FastVol.LX) < MINDIST)  // more intuitive
    dx = FastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted for staggered: %g, %g, %g\n", dx, dy, dz);

  /*
  // Check whether the point falls within a volume that is omitted
  for(int omit = 1; omit <= FastVol.NbOmitVols; ++omit)
    {
          if((dx >= (OmitVolCrnrX[omit]-FastVol.CrnrX))
                          && (dx <=
  (OmitVolCrnrX[omit]+OmitVolLX[omit]-FastVol.CrnrX))
                          && (dy >= (OmitVolCrnrY[omit]-FastVol.CrnrY))
                          && (dy <=
  (OmitVolCrnrY[omit]+OmitVolLY[omit]-FastVol.CrnrY))
                          && (dz >= (OmitVolCrnrZ[omit]-FastVol.CrnrZ))
                          && (dz <=
  (OmitVolCrnrZ[omit]+OmitVolLZ[omit]-FastVol.CrnrZ)))
                  {
                  neBEMMessage("In FastKnChPFAtPoint: point in an omitted
  volume!\n"); *Potential = 0.0; globalF->X = 0.0; globalF->Y = 0.0; globalF->Z
  = 0.0;
                  }
          }	// loop over omitted volumes

  // Find the block in which the point lies
  int thisBlock = 1;
  if(FastVol.NbBlocks > 1)
          {
          for(int block = 1; block <= FastVol.NbBlocks; ++block)
                  {
                  if(dbgFn)
                          {
                          printf("dz,(BlkCrnrZ-CornerZ),(BlkCrnrZ+BlkLZ-CornerZ):
  %lg, %lg, %lg\n", dz, (BlkCrnrZ[block]-CornerZ),
                                                          (BlkCrnrZ[block]+BlkLZ[block]-CornerZ));
                          }
                  if((dz >= (BlkCrnrZ[block]-CornerZ))
                                  && (dz <=
  (BlkCrnrZ[block]+BlkLZ[block]-CornerZ)))
                          {
                          thisBlock = block;
                          break;
                          }
                  }
          }	// if NbBlocks > 1
  */

  int thisBlock = 0;
  for (int block = 1; block <= FastVol.NbBlocks; ++block) {
    double blkBtmZ = BlkCrnrZ[block] - CornerZ;  // since CornerZ has been
    double blkTopZ = blkBtmZ + BlkLZ[block];     // subtracted from dz already
    if (dbgFn) {
      printf("block, dz, blkBtmZ, blkTopZ: %d, %lg, %lg, %lg\n", block, dz,
             blkBtmZ, blkTopZ);
    }

    // take care of difficult situations
    if ((dz <= blkBtmZ) && ((blkBtmZ - dz) < MINDIST)) dz = blkBtmZ - MINDIST;
    if ((dz >= blkBtmZ) && ((dz - blkBtmZ) < MINDIST)) dz = blkBtmZ + MINDIST;
    if ((dz <= blkTopZ) && ((blkTopZ - dz) < MINDIST)) dz = blkTopZ - MINDIST;
    if ((dz >= blkTopZ) && ((dz - blkTopZ) < MINDIST)) dz = blkTopZ + MINDIST;

    if ((dz >= blkBtmZ) && (dz <= blkTopZ)) {
      thisBlock = block;
      break;
    }
  }
  if (!thisBlock) {
    neBEMMessage("FastKnChPFAtPoint: point in none of the blocks!\n");
  }

  int nbXCells = BlkNbXCells[thisBlock];
  int nbYCells = BlkNbYCells[thisBlock];
  int nbZCells = BlkNbZCells[thisBlock];
  double delX = FastVol.LX / nbXCells;
  double delY = FastVol.LY / nbYCells;
  double delZ = BlkLZ[thisBlock] / nbZCells;
  dz -= (BlkCrnrZ[thisBlock] - CornerZ);  // distance from the block corner

  if (dbgFn) {
    printf("thisBlock: %d\n", thisBlock);
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("BlkCrnrZ: %lg\n", BlkCrnrZ[thisBlock]);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    printf("dz: %lg\n", dz);
    fflush(stdout);
  }

  // Find cell in block of basic / staggered volume within which the point lies
  int celli = (int)(dx / delX) + 1;  // Find cell in which the point lies
  if (celli < 1) {
    celli = 1;
    dx = 0.5 * delX;
    neBEMMessage("FastKnChPFAtPoint - celli < 1\n");
  }
  if (celli > nbXCells) {
    celli = nbXCells;
    dx = FastVol.LX - 0.5 * delX;
    neBEMMessage("FastKnChPFAtPoint - celli > nbXCells\n");
  }
  int cellj = (int)(dy / delY) + 1;
  if (cellj < 1) {
    cellj = 1;
    dy = 0.5 * delY;
    neBEMMessage("FastKnChPFAtPoint - cellj < 1\n");
  }
  if (cellj > nbYCells) {
    cellj = nbYCells;
    dy = FastVol.LY - 0.5 * delY;
    neBEMMessage("FastKnChPFAtPoint - cellj > nbYCells\n");
  }
  int cellk = (int)(dz / delZ) + 1;
  if (cellk < 1) {
    cellk = 1;
    dz = 0.5 * delX;
    neBEMMessage("FastKnChPFAtPoint - cellk < 1\n");
  }
  if (cellk > nbZCells) {
    cellk = nbZCells;
    dz = FastVol.LZ - 0.5 * delZ;
    neBEMMessage("FastKnChPFAtPoint - cellk > nbZCells\n");
  }
  if (dbgFn) printf("Cells in x, y, z: %d, %d, %d\n", celli, cellj, cellk);

  // Interpolate potential and field at the point using the corner values of
  // of the cell and, if necessary, of the neighbouring cells
  // These gradients can also be calculated while computing the potential and
  // field at the cells and stored in memory, provided enough memory is
  // available

  // distances from cell corner
  double dxcellcrnr = dx - (double)(celli - 1) * delX;
  double dycellcrnr = dy - (double)(cellj - 1) * delY;
  double dzcellcrnr = dz - (double)(cellk - 1) * delZ;
  if (dbgFn)
    printf("cell crnr dx, dy, dz: %g, %g, %g\n", dxcellcrnr, dycellcrnr,
           dzcellcrnr);

  // normalized distances
  double xd = dxcellcrnr / delX;  // xd = (x-x0)/(x1-x0)
  double yd = dycellcrnr / delY;  // etc
  double zd = dzcellcrnr / delZ;
  if (xd <= 0.0) xd = 0.0;
  if (yd <= 0.0) yd = 0.0;
  if (zd <= 0.0) zd = 0.0;
  if (xd >= 1.0) xd = 1.0;
  if (yd >= 1.0) yd = 1.0;
  if (zd >= 1.0) zd = 1.0;

  // corner values of potential and field
  double P000 = FastPotKnCh[thisBlock][celli][cellj][cellk];  // lowest corner
  double FX000 = FastFXKnCh[thisBlock][celli][cellj][cellk];
  double FY000 = FastFYKnCh[thisBlock][celli][cellj][cellk];
  double FZ000 = FastFZKnCh[thisBlock][celli][cellj][cellk];
  double P100 = FastPotKnCh[thisBlock][celli + 1][cellj][cellk];
  double FX100 = FastFXKnCh[thisBlock][celli + 1][cellj][cellk];
  double FY100 = FastFYKnCh[thisBlock][celli + 1][cellj][cellk];
  double FZ100 = FastFZKnCh[thisBlock][celli + 1][cellj][cellk];
  double P010 = FastPotKnCh[thisBlock][celli][cellj + 1][cellk];
  double FX010 = FastFXKnCh[thisBlock][celli][cellj + 1][cellk];
  double FY010 = FastFYKnCh[thisBlock][celli][cellj + 1][cellk];
  double FZ010 = FastFZKnCh[thisBlock][celli][cellj + 1][cellk];
  double P001 = FastPotKnCh[thisBlock][celli][cellj][cellk + 1];
  double FX001 = FastFXKnCh[thisBlock][celli][cellj][cellk + 1];
  double FY001 = FastFYKnCh[thisBlock][celli][cellj][cellk + 1];
  double FZ001 = FastFZKnCh[thisBlock][celli][cellj][cellk + 1];
  double P110 = FastPotKnCh[thisBlock][celli + 1][cellj + 1][cellk];
  double FX110 = FastFXKnCh[thisBlock][celli + 1][cellj + 1][cellk];
  double FY110 = FastFYKnCh[thisBlock][celli + 1][cellj + 1][cellk];
  double FZ110 = FastFZKnCh[thisBlock][celli + 1][cellj + 1][cellk];
  double P101 = FastPotKnCh[thisBlock][celli + 1][cellj][cellk + 1];
  double FX101 = FastFXKnCh[thisBlock][celli + 1][cellj][cellk + 1];
  double FY101 = FastFYKnCh[thisBlock][celli + 1][cellj][cellk + 1];
  double FZ101 = FastFZKnCh[thisBlock][celli + 1][cellj][cellk + 1];
  double P011 = FastPotKnCh[thisBlock][celli][cellj + 1][cellk + 1];
  double FX011 = FastFXKnCh[thisBlock][celli][cellj + 1][cellk + 1];
  double FY011 = FastFYKnCh[thisBlock][celli][cellj + 1][cellk + 1];
  double FZ011 = FastFZKnCh[thisBlock][celli][cellj + 1][cellk + 1];
  double P111 = FastPotKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FX111 = FastFXKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FY111 = FastFYKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FZ111 = FastFZKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
  if (OptStaggerFastVol) {
    if (sector == 1) {  // nothing to be done
    }
    if (sector == 2) {  // volume shifted up but point not in the staggered part
    }
    if (sector == 3) {  // staggered volume
      P000 = FastStgPotKnCh[thisBlock][celli][cellj][cellk];
      FX000 = FastStgFXKnCh[thisBlock][celli][cellj][cellk];
      FY000 = FastStgFYKnCh[thisBlock][celli][cellj][cellk];
      FZ000 = FastStgFZKnCh[thisBlock][celli][cellj][cellk];
      P100 = FastStgPotKnCh[thisBlock][celli + 1][cellj][cellk];
      FX100 = FastStgFXKnCh[thisBlock][celli + 1][cellj][cellk];
      FY100 = FastStgFYKnCh[thisBlock][celli + 1][cellj][cellk];
      FZ100 = FastStgFZKnCh[thisBlock][celli + 1][cellj][cellk];
      P010 = FastStgPotKnCh[thisBlock][celli][cellj + 1][cellk];
      FX010 = FastStgFXKnCh[thisBlock][celli][cellj + 1][cellk];
      FY010 = FastStgFYKnCh[thisBlock][celli][cellj + 1][cellk];
      FZ010 = FastStgFZKnCh[thisBlock][celli][cellj + 1][cellk];
      P001 = FastStgPotKnCh[thisBlock][celli][cellj][cellk + 1];
      FX001 = FastStgFXKnCh[thisBlock][celli][cellj][cellk + 1];
      FY001 = FastStgFYKnCh[thisBlock][celli][cellj][cellk + 1];
      FZ001 = FastStgFZKnCh[thisBlock][celli][cellj][cellk + 1];
      P110 = FastStgPotKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = FastStgFXKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = FastStgFYKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = FastStgFZKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = FastStgPotKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = FastStgFXKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = FastStgFYKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = FastStgFZKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = FastStgPotKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = FastStgFXKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = FastStgFYKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = FastStgFZKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = FastStgPotKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = FastStgFXKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = FastStgFYKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = FastStgFZKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
    if (sector == 4) {  // volume shifted down and point in the staggered part
      P000 = FastStgPotKnCh[thisBlock][celli][cellj][cellk];
      FX000 = FastStgFXKnCh[thisBlock][celli][cellj][cellk];
      FY000 = FastStgFYKnCh[thisBlock][celli][cellj][cellk];
      FZ000 = FastStgFZKnCh[thisBlock][celli][cellj][cellk];
      P100 = FastStgPotKnCh[thisBlock][celli + 1][cellj][cellk];
      FX100 = FastStgFXKnCh[thisBlock][celli + 1][cellj][cellk];
      FY100 = FastStgFYKnCh[thisBlock][celli + 1][cellj][cellk];
      FZ100 = FastStgFZKnCh[thisBlock][celli + 1][cellj][cellk];
      P010 = FastStgPotKnCh[thisBlock][celli][cellj + 1][cellk];
      FX010 = FastStgFXKnCh[thisBlock][celli][cellj + 1][cellk];
      FY010 = FastStgFYKnCh[thisBlock][celli][cellj + 1][cellk];
      FZ010 = FastStgFZKnCh[thisBlock][celli][cellj + 1][cellk];
      P001 = FastStgPotKnCh[thisBlock][celli][cellj][cellk + 1];
      FX001 = FastStgFXKnCh[thisBlock][celli][cellj][cellk + 1];
      FY001 = FastStgFYKnCh[thisBlock][celli][cellj][cellk + 1];
      FZ001 = FastStgFZKnCh[thisBlock][celli][cellj][cellk + 1];
      P110 = FastStgPotKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = FastStgFXKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = FastStgFYKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = FastStgFZKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = FastStgPotKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = FastStgFXKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = FastStgFYKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = FastStgFZKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = FastStgPotKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = FastStgFXKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = FastStgFYKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = FastStgFZKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = FastStgPotKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = FastStgFXKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = FastStgFYKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = FastStgFZKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
  }  // if OptStaggerFastVol

  double intP =
      TriLin(xd, yd, zd, P000, P100, P010, P001, P110, P101, P011, P111);
  double intFX = TriLin(xd, yd, zd, FX000, FX100, FX010, FX001, FX110, FX101,
                        FX011, FX111);
  double intFY = TriLin(xd, yd, zd, FY000, FY100, FY010, FY001, FY110, FY101,
                        FY011, FY111);
  double intFZ = TriLin(xd, yd, zd, FZ000, FZ100, FZ010, FZ001, FZ110, FZ101,
                        FZ011, FZ111);

  *Potential = intP;
  globalF->X = intFX;
  globalF->Y = intFY;
  globalF->Z = intFZ;

  if (dbgFn) {
    printf("Cell corner values:\n");
    printf("Potential: %g, %g, %g, %g\n", P000, P100, P010, P001);
    printf("Potential: %g, %g, %g, %g\n", P110, P101, P011, P111);
    printf("FastFX: %g, %g, %g, %g\n", FX000, FX100, FX010, FX001);
    printf("FastFX: %g, %g, %g, %g\n", FX110, FX101, FX011, FX111);
    printf("FastFY: %g, %g, %g, %g\n", FY000, FY100, FY010, FY001);
    printf("FastFY: %g, %g, %g, %g\n", FY110, FY101, FY011, FY111);
    printf("FastFZ: %g, %g, %g, %g\n", FZ000, FZ100, FZ010, FZ001);
    printf("FastFZ: %g, %g, %g, %g\n", FZ110, FZ101, FZ011, FZ111);
    printf("Pot, FX, FY, FZ: %g, %g, %g, %g\n", *Potential, globalF->X,
           globalF->Y, globalF->Z);
  }

  if (dbgFn) {
    printf("out FastKnChPFAtPoint\n");
    fflush(stdout);
  }

  return 0;
}  // FastKnChPFAtPoint ends

// Gives three components of the total Potential and flux in the global
// coordinate system due to all the elements using the results stored in
// the FAST volume mesh. The Fast volume is generated in the normal manner
// but by making necessary changes in the boundary conditions. This Fast
// volume is then renamed. The same is true for the data in staggered volume.
// These names are provided to the code by the neBEMWtFldFastVol.inp
int WtFldFastPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;
  double Xpt = globalP->X;
  double Ypt = globalP->Y;
  double Zpt = globalP->Z;
  double RptVolLX = WtFldFastVol.LX;
  double RptVolLY = WtFldFastVol.LY;
  double RptVolLZ = WtFldFastVol.LZ;
  double CornerX = WtFldFastVol.CrnrX;
  double CornerY = WtFldFastVol.CrnrY;
  double CornerZ = WtFldFastVol.CrnrZ;
  double TriLin(double xd, double yd, double zd, double c000, double c100,
                double c010, double c001, double c110, double c101, double c011,
                double c111);

  // First of all, check how the point in question should be treated ...

  // Check whether the point falls within a volume that is not regarded as
  // FastVol
  for (int ignore = 1; ignore <= WtFldFastVol.NbIgnoreVols; ++ignore) {
    if ((Xpt >= (WtFldIgnoreVolCrnrX[ignore])) &&
        (Xpt <= (WtFldIgnoreVolCrnrX[ignore] + WtFldIgnoreVolLX[ignore])) &&
        (Ypt >= (WtFldIgnoreVolCrnrY[ignore])) &&
        (Ypt <= (WtFldIgnoreVolCrnrY[ignore] + WtFldIgnoreVolLY[ignore])) &&
        (Zpt >= (WtFldIgnoreVolCrnrZ[ignore])) &&
        (Zpt <= (WtFldIgnoreVolCrnrZ[ignore] + WtFldIgnoreVolLZ[ignore]))) {
      if (dbgFn)
        neBEMMessage("In WtFldFastPFAtPoint: point in an ignored volume!\n");

      // KnCh does not have any effect
      int fstatus = ElePFAtPoint(globalP, Potential, globalF);
      if (fstatus != 0) {
        neBEMMessage("wrong WtFldPFAtPoint return value in FastVolPF.\n");
        return -1;
      } else
        return 0;
    }
  }  // loop over ignored volumes

  // If not ignored, the point qualifies for FastVol evaluation ...

  // for a staggered fast volume, the volume repeated in X is larger
  if (OptWtFldStaggerFastVol) {
    RptVolLX += WtFldFastVol.LX;
  }
  if (dbgFn) {
    printf("\nin WtFldFastPFAtPoint\n");
    printf("x, y, z: %g, %g, %g\n", Xpt, Ypt, Zpt);
    printf("RptVolLX, RptVolLY, RptVolLZ: %g, %g, %g\n", RptVolLX, RptVolLY,
           RptVolLZ);
    printf("CornerX, CornerY, CornerZ: %g, %g, %g\n", CornerX, CornerY,
           CornerZ);
    printf("Nb of blocks: %d\n", WtFldFastVol.NbBlocks);
    for (int block = 1; block <= WtFldFastVol.NbBlocks; ++block) {
      printf("NbOfXCells: %d\n", WtFldBlkNbXCells[block]);
      printf("NbOfYCells: %d\n", WtFldBlkNbYCells[block]);
      printf("NbOfZCells: %d\n", WtFldBlkNbZCells[block]);
      printf("LZ: %le\n", WtFldBlkLZ[block]);
      printf("CornerZ: %le\n", WtFldBlkCrnrZ[block]);
    }
  }

  // Find equivalent position inside the basic / staggered volume.
  // real distance from volume corner
  double dx = Xpt - CornerX;
  double dy = Ypt - CornerY;
  double dz = Zpt - CornerZ;
  if (dbgFn)
    printf("real dx, dy, dz from volume corner: %g, %g, %g\n", dx, dy, dz);

  int NbFastVolX = (int)(dx / RptVolLX);
  if (dx < 0.0) --NbFastVolX;
  int NbFastVolY = (int)(dy / RptVolLY);
  if (dy < 0.0) --NbFastVolY;
  int NbFastVolZ = (int)(dz / RptVolLZ);
  if (dz < 0.0) --NbFastVolZ;
  if (dbgFn)
    printf("Volumes in x, y, z: %d, %d, %d\n", NbFastVolX, NbFastVolY,
           NbFastVolZ);

  // equivalent distances from fast volume corner
  dx -= NbFastVolX * RptVolLX;
  dy -= NbFastVolY * RptVolLY;
  dz -= NbFastVolZ * RptVolLZ;
  // The following conditions should never happen - generate an error message
  if (dx < 0.0) {
    dx = 0.0;
    neBEMMessage("equiv dx < 0.0 - not correct!\n");
  }
  if (dy < 0.0) {
    dy = 0.0;
    neBEMMessage("equiv dy < 0.0 - not correct!\n");
  }
  if (dz < 0.0) {
    dz = 0.0;
    neBEMMessage("equiv dz < 0.0 - not correct!\n");
  }
  if (dx > RptVolLX) {
    dx = RptVolLX;
    neBEMMessage("equiv dx > RptVolLX - not correct!\n");
  }
  if (dy > RptVolLY) {
    dy = RptVolLY;
    neBEMMessage("equiv dy > RptVolLY - not correct!\n");
  }
  if (dz > RptVolLZ) {
    dz = RptVolLZ;
    neBEMMessage("equiv dz > RptVolLZ - not correct!\n");
  }
  if (dbgFn)
    printf("equivalent dist from corner - dx, dy, dz: %g, %g, %g\n", dx, dy,
           dz);

  // Take care of possible trouble-makers
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= WtFldFastVol.LX) && (WtFldFastVol.LX - dx) < MINDIST)
    dx = WtFldFastVol.LX - MINDIST;
  else if ((dx > WtFldFastVol.LX) && (fabs(WtFldFastVol.LX - dx) < MINDIST))
    dx = WtFldFastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted - dx, dy, dz: %g, %g, %g\n", dx, dy, dz);

  // If volume is staggered, we have a few more things to do before finalizing
  // the values of equivalent distance
  // sector identification
  // _................__________________
  // |    .     .     |    Sector 3    |
  // |     .   .      |                |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 2    |    .     .     |
  // |----------------|     .   .      |
  // |                |       .        |
  // |       .        |                |
  // |     .   .      |                |
  // |    .     .     |----------------|
  // |     .   .      |    Sector 4    |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 1    |    .     .     |
  // |----------------|................|

  int sector = 1;  // kept outside `if' since this is necessary further below
  if (OptWtFldStaggerFastVol) {
    if ((dx >= 0.0) && (dx <= WtFldFastVol.LX) && (dy >= 0.0) &&
        (dy <= WtFldFastVol.LY)) {
      // point lies in sector 1, everything remains unchanged
      sector = 1;
    } else if ((dx >= 0.0) && (dx <= WtFldFastVol.LX) &&
               (dy > WtFldFastVol.LY) &&
               (dy <= WtFldFastVol.LY + WtFldFastVol.YStagger)) {
      // point lies in sector 2, move basic volume one step up
      sector = 2;
      ++NbFastVolY;
      CornerY += WtFldFastVol.LY;  // repeat length in Y is LY
      dy -= WtFldFastVol.LY;
    } else if ((dx > WtFldFastVol.LX) && (dx <= 2.0 * WtFldFastVol.LX) &&
               (dy >= WtFldFastVol.YStagger) &&
               (dy <= WtFldFastVol.LY + WtFldFastVol.YStagger)) {
      // point lies in sector 3, pt in staggered vol, change corner coords
      sector = 3;
      CornerX += WtFldFastVol.LX;
      CornerY += WtFldFastVol.YStagger;
      dx -= WtFldFastVol.LX;
      dy -= WtFldFastVol.YStagger;
    } else if ((dx > WtFldFastVol.LX) && (dx <= 2.0 * WtFldFastVol.LX) &&
               (dy >= 0.0) && (dy < WtFldFastVol.YStagger)) {
      // point lies in sector 4, move basic volume one step down and consider
      // staggered fast volume
      sector = 4;
      --NbFastVolY;
      CornerX +=
          WtFldFastVol.LX;  // in the staggered part of the repeated volume
      CornerY -= (WtFldFastVol.LY - WtFldFastVol.YStagger);
      dx -= WtFldFastVol.LX;
      dy += (WtFldFastVol.LY - WtFldFastVol.YStagger);
    } else {
      neBEMMessage("WtFldFastPFAtPoint: point in none of the sectors!\n");
    }
    if (dbgFn) printf("stagger modified dx, dy, dz: %g, %g, %g\n", dx, dy, dz);
  }

  // Take care of possible trouble-makers - once more
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= WtFldFastVol.LX) && (WtFldFastVol.LX - dx) < MINDIST)
    dx = WtFldFastVol.LX - MINDIST;
  else if ((dx > WtFldFastVol.LX) && (fabs(WtFldFastVol.LX - dx) < MINDIST))
    dx = WtFldFastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted for staggered: %g, %g, %g\n", dx, dy, dz);

  /*
  // Check whether the point falls within a volume that is omitted
  for(int omit = 1; omit <= WtFldFastVol.NbOmitVols; ++omit)
    {
          if((dx >= (WtFldOmitVolCrnrX[omit]-WtFldFastVol.CrnrX))
                          && (dx <=
  (WtFldOmitVolCrnrX[omit]+WtFldOmitVolLX[omit]-WtFldFastVol.CrnrX))
                          && (dy >=
  (WtFldOmitVolCrnrY[omit]-WtFldFastVol.CrnrY))
                          && (dy <=
  (WtFldOmitVolCrnrY[omit]+WtFldOmitVolLY[omit]-WtFldFastVol.CrnrY))
                          && (dz >=
  (WtFldOmitVolCrnrZ[omit]-WtFldFastVol.CrnrZ))
                          && (dz <=
  (WtFldOmitVolCrnrZ[omit]+WtFldOmitVolLZ[omit]-WtFldFastVol.CrnrZ)))
                  {
                  neBEMMessage("In FastPFAtPoint: point in an omitted
  volume!\n"); *Potential = 0.0; globalF->X = 0.0; globalF->Y = 0.0; globalF->Z
  = 0.0;
                  }
          }	// loop over omitted volumes
  */

  // Find the block in which the point lies
  /*
  int thisBlock = 1;
  if(WtFldFastVol.NbBlocks > 1)
          {
          for(int block = 1; block <= WtFldFastVol.NbBlocks; ++block)
                  {
                  if(dbgFn)
                          {
                          printf("dz,(WtFldBlkCrnrZ-CornerZ),(WtFldBlkCrnrZ+WtFldBlkLZ-CornerZ):
  %lg, %lg, %lg\n", dz, (WtFldBlkCrnrZ[block]-CornerZ),
                                                          (WtFldBlkCrnrZ[block]+WtFldBlkLZ[block]-CornerZ));
                          }
                  if((dz >= (WtFldBlkCrnrZ[block]-CornerZ))
                                  && (dz <=
  (WtFldBlkCrnrZ[block]+WtFldBlkLZ[block]-CornerZ)))
                          {
                          thisBlock = block;
                          break;
                          }
                  }
          }	// if NbBlocks > 1
  */

  int thisBlock = 0;
  for (int block = 1; block <= WtFldFastVol.NbBlocks; ++block) {
    double blkBtmZ = WtFldBlkCrnrZ[block] - CornerZ;  // since CornerZ has been
    double blkTopZ = blkBtmZ + WtFldBlkLZ[block];  // subtracted from dz already
    if (dbgFn) {
      printf("block, dz, blkBtmZ, blkTopZ: %d, %lg, %lg, %lg\n", block, dz,
             blkBtmZ, blkTopZ);
    }

    // take care of difficult situations
    if ((dz <= blkBtmZ) && ((blkBtmZ - dz) < MINDIST)) dz = blkBtmZ - MINDIST;
    if ((dz >= blkBtmZ) && ((dz - blkBtmZ) < MINDIST)) dz = blkBtmZ + MINDIST;
    if ((dz <= blkTopZ) && ((blkTopZ - dz) < MINDIST)) dz = blkTopZ - MINDIST;
    if ((dz >= blkTopZ) && ((dz - blkTopZ) < MINDIST)) dz = blkTopZ + MINDIST;

    if ((dz >= blkBtmZ) && (dz <= blkTopZ)) {
      thisBlock = block;
      break;
    }
  }
  if (!thisBlock) {
    neBEMMessage("WtFldFastPFAtPoint: point in none of the blocks!\n");
  }

  int nbXCells = WtFldBlkNbXCells[thisBlock];
  int nbYCells = WtFldBlkNbYCells[thisBlock];
  int nbZCells = WtFldBlkNbZCells[thisBlock];
  double delX = WtFldFastVol.LX / nbXCells;
  double delY = WtFldFastVol.LY / nbYCells;
  double delZ = WtFldBlkLZ[thisBlock] / nbZCells;
  dz -= (WtFldBlkCrnrZ[thisBlock] - CornerZ);  // distance from the block corner

  if (dbgFn) {
    printf("thisBlock: %d\n", thisBlock);
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("WtFldBlkCrnrZ: %lg\n", WtFldBlkCrnrZ[thisBlock]);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    printf("dz: %lg\n", dz);
    fflush(stdout);
  }

  // Find cell in block of basic / staggered volume within which the point lies
  int celli = (int)(dx / delX) + 1;  // Find cell in which the point lies
  if (celli < 1) {
    celli = 1;
    dx = 0.5 * delX;
    neBEMMessage("WtFldFastPFAtPoint - celli < 1\n");
  }
  if (celli > nbXCells) {
    celli = nbXCells;
    dx = WtFldFastVol.LX - 0.5 * delX;
    neBEMMessage("WtFldFastPFAtPoint - celli > nbXCells\n");
  }
  int cellj = (int)(dy / delY) + 1;
  if (cellj < 1) {
    cellj = 1;
    dy = 0.5 * delY;
    neBEMMessage("WtFldFastPFAtPoint - cellj < 1\n");
  }
  if (cellj > nbYCells) {
    cellj = nbYCells;
    dy = WtFldFastVol.LY - 0.5 * delY;
    neBEMMessage("WtFldFastPFAtPoint - cellj > nbYCells\n");
  }
  int cellk = (int)(dz / delZ) + 1;
  if (cellk < 1) {
    cellk = 1;
    dz = 0.5 * delX;
    neBEMMessage("WtFldFastPFAtPoint - cellk < 1\n");
  }
  if (cellk > nbZCells) {
    cellk = nbZCells;
    dz = WtFldFastVol.LZ - 0.5 * delZ;
    neBEMMessage("WtFldFastPFAtPoint - cellk > nbZCells\n");
  }
  if (dbgFn) printf("Cells in x, y, z: %d, %d, %d\n", celli, cellj, cellk);

  // Interpolate potential and field at the point using the corner values of
  // of the cell and, if necessary, of the neighbouring cells
  // These gradients can also be calculated while computing the potential and
  // field at the cells and stored in memory, provided enough memory is
  // available

  // distances from cell corner
  double dxcellcrnr = dx - (double)(celli - 1) * delX;
  double dycellcrnr = dy - (double)(cellj - 1) * delY;
  double dzcellcrnr = dz - (double)(cellk - 1) * delZ;
  if (dbgFn)
    printf("cell crnr dx, dy, dz: %g, %g, %g\n", dxcellcrnr, dycellcrnr,
           dzcellcrnr);

  // normalized distances
  double xd = dxcellcrnr / delX;  // xd = (x-x0)/(x1-x0)
  double yd = dycellcrnr / delY;  // etc
  double zd = dzcellcrnr / delZ;
  if (xd <= 0.0) xd = 0.0;
  if (yd <= 0.0) yd = 0.0;
  if (zd <= 0.0) zd = 0.0;
  if (xd >= 1.0) xd = 1.0;
  if (yd >= 1.0) yd = 1.0;
  if (zd >= 1.0) zd = 1.0;

  // corner values of potential and field
  double P000 = WtFldFastPot[thisBlock][celli][cellj][cellk];  // lowest corner
  double FX000 = WtFldFastFX[thisBlock][celli][cellj][cellk];
  double FY000 = WtFldFastFY[thisBlock][celli][cellj][cellk];
  double FZ000 = WtFldFastFZ[thisBlock][celli][cellj][cellk];
  double P100 = WtFldFastPot[thisBlock][celli + 1][cellj][cellk];
  double FX100 = WtFldFastFX[thisBlock][celli + 1][cellj][cellk];
  double FY100 = WtFldFastFY[thisBlock][celli + 1][cellj][cellk];
  double FZ100 = WtFldFastFZ[thisBlock][celli + 1][cellj][cellk];
  double P010 = WtFldFastPot[thisBlock][celli][cellj + 1][cellk];
  double FX010 = WtFldFastFX[thisBlock][celli][cellj + 1][cellk];
  double FY010 = WtFldFastFY[thisBlock][celli][cellj + 1][cellk];
  double FZ010 = WtFldFastFZ[thisBlock][celli][cellj + 1][cellk];
  double P001 = WtFldFastPot[thisBlock][celli][cellj][cellk + 1];
  double FX001 = WtFldFastFX[thisBlock][celli][cellj][cellk + 1];
  double FY001 = WtFldFastFY[thisBlock][celli][cellj][cellk + 1];
  double FZ001 = WtFldFastFZ[thisBlock][celli][cellj][cellk + 1];
  double P110 = WtFldFastPot[thisBlock][celli + 1][cellj + 1][cellk];
  double FX110 = WtFldFastFX[thisBlock][celli + 1][cellj + 1][cellk];
  double FY110 = WtFldFastFY[thisBlock][celli + 1][cellj + 1][cellk];
  double FZ110 = WtFldFastFZ[thisBlock][celli + 1][cellj + 1][cellk];
  double P101 = WtFldFastPot[thisBlock][celli + 1][cellj][cellk + 1];
  double FX101 = WtFldFastFX[thisBlock][celli + 1][cellj][cellk + 1];
  double FY101 = WtFldFastFY[thisBlock][celli + 1][cellj][cellk + 1];
  double FZ101 = WtFldFastFZ[thisBlock][celli + 1][cellj][cellk + 1];
  double P011 = WtFldFastPot[thisBlock][celli][cellj + 1][cellk + 1];
  double FX011 = WtFldFastFX[thisBlock][celli][cellj + 1][cellk + 1];
  double FY011 = WtFldFastFY[thisBlock][celli][cellj + 1][cellk + 1];
  double FZ011 = WtFldFastFZ[thisBlock][celli][cellj + 1][cellk + 1];
  double P111 = WtFldFastPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FX111 = WtFldFastFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FY111 = WtFldFastFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FZ111 = WtFldFastFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
  if (OptWtFldStaggerFastVol) {
    if (sector == 1) {  // nothing to be done
    }
    if (sector == 2) {  // volume shifted up but point not in the staggered part
    }
    if (sector == 3) {  // staggered volume
      P000 = WtFldFastStgPot[thisBlock][celli][cellj][cellk];
      FX000 = WtFldFastStgFX[thisBlock][celli][cellj][cellk];
      FY000 = WtFldFastStgFY[thisBlock][celli][cellj][cellk];
      FZ000 = WtFldFastStgFZ[thisBlock][celli][cellj][cellk];
      P100 = WtFldFastStgPot[thisBlock][celli + 1][cellj][cellk];
      FX100 = WtFldFastStgFX[thisBlock][celli + 1][cellj][cellk];
      FY100 = WtFldFastStgFY[thisBlock][celli + 1][cellj][cellk];
      FZ100 = WtFldFastStgFZ[thisBlock][celli + 1][cellj][cellk];
      P010 = WtFldFastStgPot[thisBlock][celli][cellj + 1][cellk];
      FX010 = WtFldFastStgFX[thisBlock][celli][cellj + 1][cellk];
      FY010 = WtFldFastStgFY[thisBlock][celli][cellj + 1][cellk];
      FZ010 = WtFldFastStgFZ[thisBlock][celli][cellj + 1][cellk];
      P001 = WtFldFastStgPot[thisBlock][celli][cellj][cellk + 1];
      FX001 = WtFldFastStgFX[thisBlock][celli][cellj][cellk + 1];
      FY001 = WtFldFastStgFY[thisBlock][celli][cellj][cellk + 1];
      FZ001 = WtFldFastStgFZ[thisBlock][celli][cellj][cellk + 1];
      P110 = WtFldFastStgPot[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = WtFldFastStgFX[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = WtFldFastStgFY[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = WtFldFastStgFZ[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = WtFldFastStgPot[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = WtFldFastStgFX[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = WtFldFastStgFY[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = WtFldFastStgFZ[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = WtFldFastStgPot[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = WtFldFastStgFX[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = WtFldFastStgFY[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = WtFldFastStgFZ[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = WtFldFastStgPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = WtFldFastStgFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = WtFldFastStgFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = WtFldFastStgFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
    if (sector == 4) {  // volume shifted down and point in the staggered part
      P000 = WtFldFastStgPot[thisBlock][celli][cellj][cellk];
      FX000 = WtFldFastStgFX[thisBlock][celli][cellj][cellk];
      FY000 = WtFldFastStgFY[thisBlock][celli][cellj][cellk];
      FZ000 = WtFldFastStgFZ[thisBlock][celli][cellj][cellk];
      P100 = WtFldFastStgPot[thisBlock][celli + 1][cellj][cellk];
      FX100 = WtFldFastStgFX[thisBlock][celli + 1][cellj][cellk];
      FY100 = WtFldFastStgFY[thisBlock][celli + 1][cellj][cellk];
      FZ100 = WtFldFastStgFZ[thisBlock][celli + 1][cellj][cellk];
      P010 = WtFldFastStgPot[thisBlock][celli][cellj + 1][cellk];
      FX010 = WtFldFastStgFX[thisBlock][celli][cellj + 1][cellk];
      FY010 = WtFldFastStgFY[thisBlock][celli][cellj + 1][cellk];
      FZ010 = WtFldFastStgFZ[thisBlock][celli][cellj + 1][cellk];
      P001 = WtFldFastStgPot[thisBlock][celli][cellj][cellk + 1];
      FX001 = WtFldFastStgFX[thisBlock][celli][cellj][cellk + 1];
      FY001 = WtFldFastStgFY[thisBlock][celli][cellj][cellk + 1];
      FZ001 = WtFldFastStgFZ[thisBlock][celli][cellj][cellk + 1];
      P110 = WtFldFastStgPot[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = WtFldFastStgFX[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = WtFldFastStgFY[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = WtFldFastStgFZ[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = WtFldFastStgPot[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = WtFldFastStgFX[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = WtFldFastStgFY[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = WtFldFastStgFZ[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = WtFldFastStgPot[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = WtFldFastStgFX[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = WtFldFastStgFY[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = WtFldFastStgFZ[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = WtFldFastStgPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = WtFldFastStgFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = WtFldFastStgFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = WtFldFastStgFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
  }

  double intP =
      TriLin(xd, yd, zd, P000, P100, P010, P001, P110, P101, P011, P111);
  double intFX = TriLin(xd, yd, zd, FX000, FX100, FX010, FX001, FX110, FX101,
                        FX011, FX111);
  double intFY = TriLin(xd, yd, zd, FY000, FY100, FY010, FY001, FY110, FY101,
                        FY011, FY111);
  double intFZ = TriLin(xd, yd, zd, FZ000, FZ100, FZ010, FZ001, FZ110, FZ101,
                        FZ011, FZ111);

  *Potential = intP;
  globalF->X = intFX;
  globalF->Y = intFY;
  globalF->Z = intFZ;

  if (dbgFn) {
    printf("WtFldCell corner values:\n");
    printf("Potential: %g, %g, %g, %g\n", P000, P100, P010, P001);
    printf("Potential: %g, %g, %g, %g\n", P110, P101, P011, P111);
    printf("FastFX: %g, %g, %g, %g\n", FX000, FX100, FX010, FX001);
    printf("FastFX: %g, %g, %g, %g\n", FX110, FX101, FX011, FX111);
    printf("FastFY: %g, %g, %g, %g\n", FY000, FY100, FY010, FY001);
    printf("FastFY: %g, %g, %g, %g\n", FY110, FY101, FY011, FY111);
    printf("FastFZ: %g, %g, %g, %g\n", FZ000, FZ100, FZ010, FZ001);
    printf("FastFZ: %g, %g, %g, %g\n", FZ110, FZ101, FZ011, FZ111);
    printf("Pot, FX, FY, FZ: %g, %g, %g, %g\n", *Potential, globalF->X,
           globalF->Y, globalF->Z);
  }

  if (dbgFn) {
    printf("out WtFldFastPFAtPoint\n");
    fflush(stdout);
  }

  return 0;
}  // WtFldFastPFAtPoint ends

// Potential and flux per unit charge density on an element returned as
// Pot, Fx, Fy, Fz components
// in the global coordinate system
void GetPFGCS(int ele, Point3D *localP, double *Potential, Vector3D *globalF,
              DirnCosn3D *DirCos) {
  Vector3D localF;

  switch ((EleArr + ele - 1)->G.Type) {
    case 4:  // rectangular element
      RecPF(ele, localP, Potential, &localF);
      break;
    case 3:  // triangular element
      TriPF(ele, localP, Potential, &localF);
      break;
    case 2:  // linear (wire) element
      WirePF(ele, localP, Potential, &localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends

  (*globalF) = RotateVector3D(&localF, DirCos, local2global);
}  // end of GetPFGCS

// Potential and flux per unit charge density on an element returned as
// Pot, Fx, Fy, Fz components
// in the local coordiante system
void GetPF(int ele, Point3D *localP, double *Potential, Vector3D *localF) {
  switch ((EleArr + ele - 1)->G.Type) {
    case 4:  // rectangular element
      RecPF(ele, localP, Potential, localF);
      break;
    case 3:  // triangular element
      TriPF(ele, localP, Potential, localF);
      break;
    case 2:  // linear (wire) element
      WirePF(ele, localP, Potential, localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends
}  // end of GetPF

// Flux per unit charge density on a rectangular element
// Are X and Z directions the same as obtained using the direction cosines?
void RecPF(int ele, Point3D *localP, double *Potential, Vector3D *localF) {
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  double diag = sqrt(a * a + b * b);  // diagonal of the element

  if (dist >= FarField * diag) {
    double dA = a * b;  // area of the rectangular element
    (*Potential) = dA / dist;
    double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    int fstatus =
        ExactRecSurf(xpt / a, ypt / a, zpt / a, -0.5, -(b / a) / 2.0,
                     0.5, (b / a) / 2.0, Potential, localF);
    if (fstatus) { // non-zero
      printf("problem in RecPF ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
    (*Potential) *= a;  // rescale - cannot be done outside because of the `if'
  }

  // Potential, Ex, Ey, Ez
#ifdef __cplusplus
  (*Potential) *= InvFourPiEps0;
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  (*Potential) /= MyFACTOR;
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
}  // end of RecPF

// Flux per unit charge density on a triangluar element
void TriPF(int ele, Point3D *localP, double *Potential, Vector3D *localF) {
  // printf("In TriPF\n");
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  // longest side (hypotenuse) of the element
  double diag = sqrt(a * a + b * b);  

  const double xm = xpt - a / 3.;
  const double zm = zpt - b / 3.;
  double dist = sqrt(xm * xm + ypt * ypt + zm * zm);

  if (dist >= FarField * diag) {
    double dA = 0.5 * a * b;  // area of the triangular element
    (*Potential) = dA / dist;
    double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    int fstatus =
        ExactTriSurf(b / a, xpt / a, ypt / a, zpt / a, Potential, localF);
    // fstatus = ApproxTriSurf(b/a, X/a, Y/a, Z/a, 5000, 5000, &Pot, &Flux);
    if (fstatus) { // non-zero
      printf("problem in TriPF ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
    (*Potential) *= a;  // rescale - cannot be done outside because of the `if'
  }

#ifdef __cplusplus
  (*Potential) *= InvFourPiEps0;
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  (*Potential) /= MyFACTOR;
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
  // printf("Out of TriPF\n");
}  // end of TriPF

// Flux per unit charge density on a wire element
void WirePF(int ele, Point3D *localP, double *Potential, Vector3D *localF) {
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double rW = (EleArr + ele - 1)->G.LX;
  double lW = (EleArr + ele - 1)->G.LZ;
  double dist =
      sqrt(xpt * xpt + ypt * ypt + zpt * zpt);  // fld pt from ele cntrd

  if (dist >= FarField * lW) { // all are distances and, hence, +ve
    double dA = 2.0 * ST_PI * rW * lW;
    (*Potential) = dA / dist;
    double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST)) {
      if (fabs(zpt) < MINDIST)
        (*Potential) = ExactCentroidalP_W(rW, lW);
      else
        (*Potential) = ExactAxialP_W(rW, lW, zpt);

      localF->X = localF->Y = 0.0;
      localF->Z = ExactThinFZ_W(rW, lW, xpt, ypt, zpt);
    } else {
      ExactThinWire(rW, lW, xpt, ypt, zpt, Potential, localF);
    }
  }

#ifdef __cplusplus
  (*Potential) *= InvFourPiEps0;
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  (*Potential) /= MyFACTOR;
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
}  // end of WirePF

// Potential and flux per unit charge density on an element returned as
// Pot, Fx, Fy, Fz components
// in the global coordiante system
void GetPrimPFGCS(int prim, Point3D *localP, double *Potential,
                  Vector3D *globalF, DirnCosn3D *DirCos) {
  Vector3D localF;

  switch (PrimType[prim]) {
    case 4:  // rectangular primitive
      RecPrimPF(prim, localP, Potential, &localF);
      break;
    case 3:  // triangular primitive
      TriPrimPF(prim, localP, Potential, &localF);
      break;
    case 2:  // linear (wire) primitive
      WirePrimPF(prim, localP, Potential, &localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends

  (*globalF) = RotateVector3D(&localF, DirCos, local2global);
}  // end of GetPrimPFGCS

// Potential and flux per unit charge density on an element returned as
// Pot, Fx, Fy, Fz components
// in the local coordiante system
void GetPrimPF(int prim, Point3D *localP, double *Potential, Vector3D *localF) {
  switch (PrimType[prim]) {
    case 4:  // rectangular primitive
      RecPrimPF(prim, localP, Potential, localF);
      break;
    case 3:  // triangular primitive
      TriPrimPF(prim, localP, Potential, localF);
      break;
    case 2:  // linear (wire) primitive
      WirePrimPF(prim, localP, Potential, localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends
}  // end of GetPrimPF

// Flux per unit charge density on a rectangular primitive
void RecPrimPF(int prim, Point3D *localP, double *Potential, Vector3D *localF) {
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  double a = PrimLX[prim];
  double b = PrimLZ[prim];
  double diag = sqrt(a * a + b * b);  // diagonal

  if (dist >= FarField * diag)  // all are distances and, hence, +ve
  {
    double dA = a * b;  // area
    (*Potential) = dA / dist;
    const double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    int fstatus =
        ExactRecSurf(xpt / a, ypt / a, zpt / a, -0.5, -(b / a) / 2.0,
                     0.5, (b / a) / 2.0, Potential, localF);
    if (fstatus) { // non-zero
      printf("problem in RecPrimPF ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
    (*Potential) *= a;  // rescale - cannot be done outside because of the `if'
  }

  // Potential, Ex, Ey, Ez
#ifdef __cplusplus
  (*Potential) *= InvFourPiEps0;
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  (*Potential) /= MyFACTOR;
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
}  // end of RecPrimPF

// Flux per unit charge density on a triangluar primitive
// Note that vertex[1] is the right angle corner
void TriPrimPF(int prim, Point3D *localP, double *Potential, Vector3D *localF) {
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double a = PrimLX[prim];
  double b = PrimLZ[prim];
  // longest side (hypotenuse)
  double diag = sqrt(a * a + b * b);  
  const double xm = xpt - a / 3.;
  const double zm = zpt - b / 3.;
  double dist = sqrt(xm * xm + ypt * ypt + zm * zm);

  if (dist >= FarField * diag) {
    double dA = 0.5 * a * b;  // area
    (*Potential) = dA / dist;
    double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    int fstatus =
        ExactTriSurf(b / a, xpt / a, ypt / a, zpt / a, Potential, localF);
    if (fstatus)  // non-zero
    {
      printf("problem in TriPrimPF ... \n");
    }
    (*Potential) *= a;  // rescale - cannot be done outside because of the `if'
  }

#ifdef __cplusplus
  (*Potential) *= InvFourPiEps0;
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  (*Potential) /= MyFACTOR;
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
  // printf("Out of TriPrimPF\n");
}  // end of TriPrimPF

// Flux per unit charge density on a wire primitive
void WirePrimPF(int prim, Point3D *localP, double *Potential,
                Vector3D *localF) {
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double rW = Radius[prim];
  double lW = PrimLZ[prim];
  double dist =
      sqrt(xpt * xpt + ypt * ypt + zpt * zpt);  // fld pt from ele cntrd

  if (dist >= FarField * lW)  // all are distances and, hence, +ve
  {
    double dA = 2.0 * ST_PI * rW * lW;
    (*Potential) = dA / dist;
    double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST)) {
      if (fabs(zpt) < MINDIST)
        (*Potential) = ExactCentroidalP_W(rW, lW);
      else
        (*Potential) = ExactAxialP_W(rW, lW, zpt);

      localF->X = localF->Y = 0.0;
      localF->Z = ExactThinFZ_W(rW, lW, xpt, ypt, zpt);
    } else {
      ExactThinWire(rW, lW, xpt, ypt, zpt, Potential, localF);
    }
  }

#ifdef __cplusplus
  (*Potential) *= InvFourPiEps0;
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  (*Potential) /= MyFACTOR;
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
}  // end of WirePrimPF

/*
// Gives three components of weighting field in the global coordinate system
// due to all the elements
// Note that local evaluation of influence and additional influences have not
// been incorporated here. Iff local evaluation show a substantial advantage
// over the cleaner function call, we'll implement the former in this function.
// This function should be merged with PFAtPoint since the only change is in
// the use of weighting field charge density instead of the physical charge
// denstiy. However, care should be taken to check the last to points mentioned
// in this function - VSystemChargeZero and effects of known charge densities.
int WtPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF,
                                                                int IdWtField)
{
int primsrc;
Point3D localP;
Point3D fldpt;
double tmpPot;
Vector3D localF, tmpF;
double xsrc, ysrc, zsrc;

// Compute field at different locations
(*Potential) = globalF->X = globalF->Y = globalF->Z = 0.0;
double xfld = globalP->X;
double yfld = globalP->Y;
double zfld = globalP->Z;
fldpt.X = xfld; fldpt.Y = yfld; fldpt.Z = zfld;

for(unsigned int ele = 1; ele <= NbElements; ++ele)
        {
        localF.X = localF.Y = localF.Z = 0.0;

  primsrc = (EleArr+ele-1)->PrimitiveNb; // not just convenience
  xsrc = (EleArr+ele-1)->G.Origin.X;	// allows faster computation
  ysrc = (EleArr+ele-1)->G.Origin.Y;
  zsrc = (EleArr+ele-1)->G.Origin.Z;

        {	// Rotate point3D from global to local system
        double InitialVector[4];
        double TransformationMatrix[4][4] = {{0.0, 0.0, 0.0, 0.0},
                                                                                        {0.0, 0.0, 0.0, 0.0},
                                                                                                {0.0, 0.0, 0.0, 0.0},
                                                                                                {0.0, 0.0, 0.0, 1.0}};
        DirnCosn3D *DirCos = &(EleArr+ele-1)->G.DC;
        double FinalVector[4];

        InitialVector[0] = xfld - xsrc; InitialVector[1] = yfld - ysrc;
        InitialVector[2] = zfld - zsrc; InitialVector[3] = 1.0;

  TransformationMatrix[0][0] = DirCos->XUnit.X;
  TransformationMatrix[0][1] = DirCos->XUnit.Y;
  TransformationMatrix[0][2] = DirCos->XUnit.Z;
  TransformationMatrix[1][0] = DirCos->YUnit.X;
  TransformationMatrix[1][1] = DirCos->YUnit.Y;
  TransformationMatrix[1][2] = DirCos->YUnit.Z;
  TransformationMatrix[2][0] = DirCos->ZUnit.X;
  TransformationMatrix[2][1] = DirCos->ZUnit.Y;
  TransformationMatrix[2][2] = DirCos->ZUnit.Z;

        for(int i = 0; i < 4; ++i)
                {
                FinalVector[i] = 0.0;
                for(int j = 0 ; j < 4; ++j)
                {
                FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
                }
                }

        localP.X = FinalVector[0];
        localP.Y = FinalVector[1];
        localP.Z = FinalVector[2];
        }	// Point3D rotated

        // Potential and flux (local system) due to base primitive
        GetPF(ele, &localP, &tmpPot, &tmpF);
  (*Potential) += WtFieldChDen[IdWtField][ele] * tmpPot;
  localF.X += WtFieldChDen[IdWtField][ele] * tmpF.X;
  localF.Y += WtFieldChDen[IdWtField][ele] * tmpF.Y;
  localF.Z += WtFieldChDen[IdWtField][ele] * tmpF.Z;

        {	// Mirror effect of base primitive
        Point3D srcpt;
        DirnCosn3D DirCos;

  srcpt.X = xsrc; srcpt.Y = ysrc; srcpt.Z = zsrc;

  if(MirrorTypeX[primsrc])
    { MirrorTypeY[primsrc] = 0; MirrorTypeZ[primsrc] = 0; }
  if(MirrorTypeY[primsrc]) MirrorTypeZ [primsrc]= 0;

  if(MirrorTypeX[primsrc])
    {
    localP = ReflectOnMirror('X', ele, srcpt, fldpt,
                              MirrorDistXFromOrigin[primsrc], &DirCos);
                GetPFGCS(ele, &localP, &tmpPot, &tmpF, &DirCos);

    if(MirrorTypeX[primsrc] == 1) // opposite charge density
                        {
                        (*Potential) -= WtFieldChDen[IdWtField][ele] * tmpPot;
                        globalF->X -= WtFieldChDen[IdWtField][ele] * tmpF.X;
                        globalF->Y -= WtFieldChDen[IdWtField][ele] * tmpF.Y;
                        globalF->Z -= WtFieldChDen[IdWtField][ele] * tmpF.Z;
                        }
    if(MirrorTypeX[primsrc] == 2) // same charge density
                        {
                        (*Potential) += WtFieldChDen[IdWtField][ele] * tmpPot;
                        globalF->X += WtFieldChDen[IdWtField][ele] * tmpF.X;
                        globalF->Y += WtFieldChDen[IdWtField][ele] * tmpF.Y;
                        globalF->Z += WtFieldChDen[IdWtField][ele] * tmpF.Z;
                        }
     }

  if(MirrorTypeY[primsrc])
    {
    localP = ReflectOnMirror('Y', ele, srcpt, fldpt,
                              MirrorDistYFromOrigin[primsrc], &DirCos);
                GetPFGCS(ele, &localP, &tmpPot, &tmpF, &DirCos);

    if(MirrorTypeY[primsrc] == 1) // opposite charge density
                        {
                        (*Potential) -= WtFieldChDen[IdWtField][ele] * tmpPot;
                        globalF->X -= WtFieldChDen[IdWtField][ele] * tmpF.X;
                        globalF->Y -= WtFieldChDen[IdWtField][ele] * tmpF.Y;
                        globalF->Z -= WtFieldChDen[IdWtField][ele] * tmpF.Z;
                        }
    if(MirrorTypeY[primsrc] == 2) // same charge density
                        {
                        (*Potential) += WtFieldChDen[IdWtField][ele] * tmpPot;
                        globalF->X += WtFieldChDen[IdWtField][ele] * tmpF.X;
                        globalF->Y += WtFieldChDen[IdWtField][ele] * tmpF.Y;
                        globalF->Z += WtFieldChDen[IdWtField][ele] * tmpF.Z;
                        }
    }

  if(MirrorTypeZ[primsrc])
    {
    localP = ReflectOnMirror('Z', ele, srcpt, fldpt,
                              MirrorDistZFromOrigin[primsrc], &DirCos);
                GetPFGCS(ele, &localP, &tmpPot, &tmpF, &DirCos);

    if(MirrorTypeZ[primsrc] == 1) // opposite charge density
                        {
                        (*Potential) -= WtFieldChDen[IdWtField][ele] * tmpPot;
                        globalF->X -= WtFieldChDen[IdWtField][ele] * tmpF.X;
                        globalF->Y -= WtFieldChDen[IdWtField][ele] * tmpF.Y;
                        globalF->Z -= WtFieldChDen[IdWtField][ele] * tmpF.Z;
                        }
    if(MirrorTypeZ[primsrc] == 2) // same charge density
                        {
                        (*Potential) += WtFieldChDen[IdWtField][ele] * tmpPot;
                        globalF->X += WtFieldChDen[IdWtField][ele] * tmpF.X;
                        globalF->Y += WtFieldChDen[IdWtField][ele] * tmpF.Y;
                        globalF->Z += WtFieldChDen[IdWtField][ele] * tmpF.Z;
                        }
                }
        }	// Mirror effect ends

        // Flux due to repeated primitives
        int prim = (EleArr+ele-1)->PrimitiveNb;
  if((PeriodicTypeX[prim] == 1) || (PeriodicTypeY[prim] == 1)
      || (PeriodicTypeZ[prim] == 1))
    {
    if(PeriodicInX[prim] || PeriodicInY[prim] || PeriodicInZ[prim])
      {
      double XOfRpt, YOfRpt, ZOfRpt;

      for(int xrpt = -PeriodicInX[prim];
            xrpt <= PeriodicInX[prim]; ++xrpt)
        {
        XOfRpt = xsrc + XPeriod[prim] * (double)xrpt;
        for(int yrpt = -PeriodicInY[prim];
              yrpt <= PeriodicInY[prim]; ++yrpt)
          {
          YOfRpt = ysrc + YPeriod[prim] * (double)yrpt;

          for(int zrpt = -PeriodicInZ[prim];
                zrpt <= PeriodicInZ[prim]; ++zrpt)
            {
            ZOfRpt = zsrc + ZPeriod[prim] * (double)zrpt;

            if( (xrpt == 0) && (yrpt == 0) && (zrpt == 0) )
                continue; // this is the base device

                                                {	// Rotate point3D from
global to local system double InitialVector[4]; double
TransformationMatrix[4][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0,
0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 1.0}}; DirnCosn3D *DirCos =
&(EleArr+ele-1)->G.DC; double FinalVector[4];

                                                InitialVector[0] = xfld -
XOfRpt; InitialVector[1] = yfld - YOfRpt; InitialVector[2] = zfld - ZOfRpt;
InitialVector[3] = 1.0;

                                        TransformationMatrix[0][0] =
DirCos->XUnit.X; TransformationMatrix[0][1] = DirCos->XUnit.Y;
                                        TransformationMatrix[0][2] =
DirCos->XUnit.Z; TransformationMatrix[1][0] = DirCos->YUnit.X;
                                        TransformationMatrix[1][1] =
DirCos->YUnit.Y; TransformationMatrix[1][2] = DirCos->YUnit.Z;
                                        TransformationMatrix[2][0] =
DirCos->ZUnit.X; TransformationMatrix[2][1] = DirCos->ZUnit.Y;
                                        TransformationMatrix[2][2] =
DirCos->ZUnit.Z;

                                                for(int i = 0; i < 4;
++i)
                                                        {
                                                        FinalVector[i] = 0.0;
                                                        for(int j = 0 ;
j < 4; ++j)
                                                        {
                                                        FinalVector[i] +=
TransformationMatrix[i][j] * InitialVector[j];
                                                        }
                                                        }

                                                localP.X = FinalVector[0];
                                                localP.Y = FinalVector[1];
                                                localP.Z = FinalVector[2];
                                                }	// Point3D rotated

                GetPF(ele, &localP, &tmpPot, &tmpF);
            (*Potential) += WtFieldChDen[IdWtField][ele] * tmpPot;
            localF.X += WtFieldChDen[IdWtField][ele] * tmpF.X;
            localF.Y += WtFieldChDen[IdWtField][ele] * tmpF.Y;
            localF.Z += WtFieldChDen[IdWtField][ele] * tmpF.Z;

                                                {	// Mirror effect of
repetition Point3D srcpt; DirnCosn3D DirCos;

                                        srcpt.X = XOfRpt; srcpt.Y = YOfRpt;
srcpt.Z = ZOfRpt;

                                        if(MirrorTypeX[primsrc])
                                        { MirrorTypeY[primsrc] = 0;
MirrorTypeZ[primsrc] = 0; } if(MirrorTypeY[primsrc]) MirrorTypeZ [primsrc]= 0;

                                        if(MirrorTypeX[primsrc])
                                        {
                                        localP = ReflectOnMirror('X', ele,
srcpt, fldpt, MirrorDistXFromOrigin[primsrc], &DirCos); GetPFGCS(ele, &localP,
&tmpPot, &tmpF, &DirCos);

                                        if(MirrorTypeX[primsrc] == 1) //
opposite charge density
                                                                {
                                                                (*Potential) -=
WtFieldChDen[IdWtField][ele] * tmpPot; globalF->X -=
WtFieldChDen[IdWtField][ele] * tmpF.X; globalF->Y -=
WtFieldChDen[IdWtField][ele] * tmpF.Y; globalF->Z -=
WtFieldChDen[IdWtField][ele] * tmpF.Z;
                                                                }
                                        if(MirrorTypeX[primsrc] == 2) // same
charge density
                                                                {
                                                                (*Potential) +=
WtFieldChDen[IdWtField][ele] * tmpPot; globalF->X +=
WtFieldChDen[IdWtField][ele] * tmpF.X; globalF->Y +=
WtFieldChDen[IdWtField][ele] * tmpF.Y; globalF->Z +=
WtFieldChDen[IdWtField][ele] * tmpF.Z;
                                                                }
                                        }

                                        if(MirrorTypeY[primsrc])
                                        {
                                        localP = ReflectOnMirror('Y', ele,
srcpt, fldpt, MirrorDistYFromOrigin[primsrc], &DirCos); GetPFGCS(ele, &localP,
&tmpPot, &tmpF, &DirCos);

                                        if(MirrorTypeY[primsrc] == 1) //
opposite charge density
                                                                {
                                                                (*Potential) -=
WtFieldChDen[IdWtField][ele] * tmpPot; globalF->X -=
WtFieldChDen[IdWtField][ele] * tmpF.X; globalF->Y -=
WtFieldChDen[IdWtField][ele] * tmpF.Y; globalF->Z -=
WtFieldChDen[IdWtField][ele] * tmpF.Z;
                                                                }
                                        if(MirrorTypeY[primsrc] == 2) // same
charge density
                                                                {
                                                                (*Potential) +=
WtFieldChDen[IdWtField][ele] * tmpPot; globalF->X +=
WtFieldChDen[IdWtField][ele] * tmpF.X; globalF->Y +=
WtFieldChDen[IdWtField][ele] * tmpF.Y; globalF->Z +=
WtFieldChDen[IdWtField][ele] * tmpF.Z;
                                                                }
                                        }

                                        if(MirrorTypeZ[primsrc])
                                        {
                                        localP = ReflectOnMirror('Z', ele,
srcpt, fldpt, MirrorDistZFromOrigin[primsrc], &DirCos); GetPFGCS(ele, &localP,
&tmpPot, &tmpF, &DirCos);

                                        if(MirrorTypeZ[primsrc] == 1) //
opposite charge density
                                                                {
                                                                (*Potential) -=
WtFieldChDen[IdWtField][ele] * tmpPot; globalF->X -=
WtFieldChDen[IdWtField][ele] * tmpF.X; globalF->Y -=
WtFieldChDen[IdWtField][ele] * tmpF.Y; globalF->Z -=
WtFieldChDen[IdWtField][ele] * tmpF.Z;
                                                                }
                                        if(MirrorTypeZ[primsrc] == 2) // same
charge density
                                                                {
                                                                (*Potential) +=
WtFieldChDen[IdWtField][ele] * tmpPot; globalF->X +=
WtFieldChDen[IdWtField][ele] * tmpF.X; globalF->Y +=
WtFieldChDen[IdWtField][ele] * tmpF.Y; globalF->Z +=
WtFieldChDen[IdWtField][ele] * tmpF.Z;
                                                                }
                                                        }
                                                }	// Mirror effect ends

            } // for zrpt
          } // for yrpt
        } // for xrpt
      } // PeriodicInX || PeriodicInY || PeriodicInZ
    } // PeriodicType == 1

        tmpF = RotateVector3D(&localF, &(EleArr+ele-1)->G.DC, local2global);
        globalF->X += tmpF.X;
        globalF->Y += tmpF.Y;
        globalF->Z += tmpF.Z;
        }	// for ele

// Check the following addition - is it necessary for weighting Potential?
(*Potential) += VSystemChargeZero;	// respect total system charge
constraint

// effect of known charges?

return(0);
}	// end of WtPFAtPoint
*/

// Gives three components of weighting field in the global coordinate system
// due to all the elements
// Note that local evaluation of influence and additional influences have not
// been incorporated here. Iff local evaluation show a substantial advantage
// over the cleaner function call, we'll implement the former in this function.
// This function should be merged with PFAtPoint since the only change is in
// the use of weighting field charge density instead of the physical charge
// denstiy. However, care should be taken to check the last to points mentioned
// in this function - VSystemChargeZero and effects of known charge densities.
// Multi-threading implemented in the following routine.
// Gives three components of the total Potential and flux in the global
// coordinate system due to all the elements
int WtPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF,
                int IdWtField) {
  int dbgFn = 0;

  const double xfld = globalP->X;
  const double yfld = globalP->Y;
  const double zfld = globalP->Z;

  // Compute Potential and field at different locations
  *Potential = globalF->X = globalF->Y = globalF->Z = 0.0;

  // Effects due to base primitives and their repetitions are considered in the
  // local coordinate system of the primitive (or element), while effects due to
  // mirror elements and their repetitions are considered in the global
  // coordinate system (GCS). This works because the direction cosines of a
  // primitive (and its elements) and those of its repetitions are the same.
  // As a result, we can do just one transformation from local to global at the
  // end of calculations related to a primitive. This can save substantial
  // computation if a discretized version of the primitive is being used since
  // we avoid one unnecessary transformation for each element that comprises a
  // primitive.
  // Begin with primitive description of the device

  // Scope in OpenMP: Variables in the global data space are accessible to all
  // threads, while variables in a thread's private space is accessible to the
  // thread only (there are several variations - copying outside region etc)
  // Field point remains the same - kept outside private
  // source point changes with change in primitive - private
  // TransformationMatrix changes - kept within private (Note: matrices with
  // fixed dimensions can be maintained, but those with dynamic allocation
  // can not).
  double *pPot = dvector(1, NbPrimitives);
  double *plFx = dvector(1, NbPrimitives);  // field components in LCS
  double *plFy = dvector(1, NbPrimitives);  // for a primitive
  double *plFz = dvector(1, NbPrimitives);  // and its other incarnations

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    pPot[prim] = plFx[prim] = plFy[prim] = plFz[prim] = 0.0;
  }

#ifdef _OPENMP
  int tid = 0, nthreads = 1;
  #pragma omp parallel private(tid, nthreads)
#endif
  {

#ifdef _OPENMP
    if (dbgFn) {
      tid = omp_get_thread_num();
      if (tid == 0) {
        nthreads = omp_get_num_threads();
        printf("PFAtPoint computation with %d threads\n", nthreads);
      }
    }
#endif
// by default, nested parallelization is off in C
#ifdef _OPENMP
#pragma omp for
#endif
    for (int primsrc = 1; primsrc <= NbPrimitives; ++primsrc) {
      if (dbgFn) {
        printf("Evaluating effect of primsrc %d using on %lg, %lg, %lg\n",
               primsrc, xfld, yfld, zfld);
        fflush(stdout);
      }

      const double xpsrc = PrimOriginX[primsrc];
      const double ypsrc = PrimOriginY[primsrc];
      const double zpsrc = PrimOriginZ[primsrc];

      // Field in the local frame.
      double lFx = 0.;
      double lFy = 0.;
      double lFz = 0.;

      // Set up transform matrix for this primitive, which is also the same
      // for all the elements belonging to this primitive.
      double TransformationMatrix[3][3];
      TransformationMatrix[0][0] = PrimDC[primsrc].XUnit.X;
      TransformationMatrix[0][1] = PrimDC[primsrc].XUnit.Y;
      TransformationMatrix[0][2] = PrimDC[primsrc].XUnit.Z;
      TransformationMatrix[1][0] = PrimDC[primsrc].YUnit.X;
      TransformationMatrix[1][1] = PrimDC[primsrc].YUnit.Y;
      TransformationMatrix[1][2] = PrimDC[primsrc].YUnit.Z;
      TransformationMatrix[2][0] = PrimDC[primsrc].ZUnit.X;
      TransformationMatrix[2][1] = PrimDC[primsrc].ZUnit.Y;
      TransformationMatrix[2][2] = PrimDC[primsrc].ZUnit.Z;

      // The total influence is due to primitives on the basic device and due to
      // virtual primitives arising out of repetition, reflection etc and not
      // residing on the basic device

      {  // basic primitive
        // point translated to the ECS origin, but axes direction global
        Point3D localPP;
        {  // Rotate point3D from global to local system
          double InitialVector[3] = {xfld - xpsrc, yfld - ypsrc, zfld - zpsrc};
          double FinalVector[3] = {0., 0., 0.};
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
            }
          }
          localPP.X = FinalVector[0];
          localPP.Y = FinalVector[1];
          localPP.Z = FinalVector[2];
        }  // Point3D rotated

        // evaluate possibility whether primitive influence is accurate enough
        // This could be based on localPP and the subtended solid angle
        // If 1, then only primitive influence will be considered
        int PrimOK = 0;
        if (PrimOK) {
          // Potential and flux (local system) due to base primitive
          double tmpPot;
          Vector3D tmpF;
          GetPrimPF(primsrc, &localPP, &tmpPot, &tmpF);
          const double qpr = AvWtChDen[IdWtField][primsrc];
          pPot[primsrc] += qpr * tmpPot;
          lFx += qpr * tmpF.X;
          lFy += qpr * tmpF.Y;
          lFz += qpr * tmpF.Z;
          // if(DebugLevel == 301)
          if (dbgFn) {
            printf("PFAtPoint base primitive =>\n");
            printf("primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n",
                   primsrc, localPP.X, localPP.Y, localPP.Z);
            printf("primsrc: %d, Pot: %lg, Fx: %lg, Fx: %lg, Fz: %lg\n",
                   primsrc, tmpPot, tmpF.X, tmpF.Y, tmpF.Z);
            printf("primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
                   primsrc, pPot[primsrc], lFx, lFy, lFz);
            fflush(stdout);
            // exit(-1);
          }
        } else {
          // element influence
          double tPot;
          Vector3D tF;
          double ePot = 0.;
          Vector3D eF;
          eF.X = 0.0;
          eF.Y = 0.0;
          eF.Z = 0.0;

          const int eleMin = ElementBgn[primsrc];
          const int eleMax = ElementEnd[primsrc];
          for (int ele = eleMin; ele <= eleMax; ++ele) {
            const double xsrc = (EleArr + ele - 1)->G.Origin.X;
            const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
            const double zsrc = (EleArr + ele - 1)->G.Origin.Z;
            // Rotate point3D from global to local system; matrix as for
            // primitive
            double InitialVector[3] = {xfld - xsrc, yfld - ysrc, zfld - zsrc};
            double FinalVector[3] = {0., 0., 0.};
            for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                FinalVector[i] +=
                    TransformationMatrix[i][j] * InitialVector[j];
              }
            }
            Point3D localPE;
            localPE.X = FinalVector[0];
            localPE.Y = FinalVector[1];
            localPE.Z = FinalVector[2];

            // Potential and flux (local system) due to base primitive
            GetPF(ele, &localPE, &tPot, &tF);
            const double qel = WtFieldChDen[IdWtField][ele];
            ePot += qel * tPot;
            eF.X += qel * tF.X;
            eF.Y += qel * tF.Y;
            eF.Z += qel * tF.Z;
            // if(DebugLevel == 301)
            if (dbgFn) {
              printf("PFAtPoint base primitive:%d\n", primsrc);
              printf("ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n", ele,
                     localPE.X, localPE.Y, localPE.Z);
              printf(
                  "ele: %d, tPot: %lg, tFx: %lg, tFy: %lg, tFz: %lg, Solution: "
                  "%g\n",
                  ele, tPot, tF.X, tF.Y, tF.Z, qel);
              printf("ele: %d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n", ele,
                     ePot, eF.X, eF.Y, eF.Z);
              fflush(stdout);
            }
          }  // for all the elements on this primsrc primitive

          pPot[primsrc] += ePot;
          lFx += eF.X;
          lFy += eF.Y;
          lFz += eF.Z;
          if (dbgFn) {
            printf(
                "prim%d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n",
                primsrc, ePot, eF.X, eF.Y, eF.Z);
            printf("prim%d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n", primsrc,
                   pPot[primsrc], lFx, lFy, lFz);
            fflush(stdout);
          }
        }  // else elements influence

        // if(DebugLevel == 301)
        if (dbgFn) {
          printf("basic primtive\n");
          printf("primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
                 primsrc, pPot[primsrc], lFx, lFy, lFz);
          fflush(stdout);
        }
      }  // basic primitive ends

      if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
          MirrorTypeZ[primsrc]) {  // Mirror effect of base primitives
        printf("Mirror may not be correctly implemented ...\n");
        exit(0);
      }        // Mirror effect ends

      // Flux due to repeated primitives
      if ((PeriodicTypeX[primsrc] == 1) || (PeriodicTypeY[primsrc] == 1) ||
          (PeriodicTypeZ[primsrc] == 1)) {
        const int perx = PeriodicInX[primsrc];
        const int pery = PeriodicInY[primsrc];
        const int perz = PeriodicInZ[primsrc];
        if (perx || pery || perz) {
          for (int xrpt = -perx; xrpt <= perx; ++xrpt) {
            const double xShift = XPeriod[primsrc] * (double)xrpt;
            double XPOfRpt = xpsrc + xShift;
            for (int yrpt = -pery; yrpt <= pery; ++yrpt) {
              const double yShift = YPeriod[primsrc] * (double)yrpt;
              double YPOfRpt = ypsrc + yShift;
              for (int zrpt = -perz; zrpt <= perz; ++zrpt) {
                const double zShift = ZPeriod[primsrc] * (double)zrpt;
                double ZPOfRpt = zpsrc + zShift;
                // Skip the base device.
                if ((xrpt == 0) && (yrpt == 0) && (zrpt == 0)) continue;
                {  // basic primitive repeated
                  Point3D localPPR;
                  {  // Rotate point3D from global to local system
                    double InitialVector[3] = {xfld - XPOfRpt, yfld - YPOfRpt, zfld - ZPOfRpt};
                    double FinalVector[3] = {0., 0., 0.};
                    for (int i = 0; i < 3; ++i) {
                      for (int j = 0; j < 3; ++j) {
                        FinalVector[i] +=
                            TransformationMatrix[i][j] * InitialVector[j];
                      }
                    }
                    localPPR.X = FinalVector[0];
                    localPPR.Y = FinalVector[1];
                    localPPR.Z = FinalVector[2];
                  }  // Point3D rotated

                  int PrimOK = 0;

                  // consider primitive representation accurate enough if it is
                  // repeated and beyond PrimAfter repetitions.
                  if (PrimAfter == 0) {
                    // If PrimAfter is zero, PrimOK is always zero
                    PrimOK = 0;
                  } else if ((abs(xrpt) > PrimAfter) && (abs(yrpt) > PrimAfter)) {
                    PrimOK = 1;
                  }
                  if (PrimOK) {  // use primitive representation
                    // Potential and flux (local system) due to repeated
                    // primitive
                    double tmpPot;
                    Vector3D tmpF;
                    GetPrimPF(primsrc, &localPPR, &tmpPot, &tmpF);
                    const double qpr = AvWtChDen[IdWtField][primsrc];
                    pPot[primsrc] += qpr * tmpPot;
                    lFx += qpr * tmpF.X;
                    lFy += qpr * tmpF.Y;
                    lFz += qpr * tmpF.Z;
                    // if(DebugLevel == 301)
                    if (dbgFn) {
                      printf(
                          "primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal: "
                          "%lg\n",
                          primsrc, localPPR.X, localPPR.Y, localPPR.Z);
                      printf(
                          "primsrc: %d, Pot: %lg, Fx: %lg, Fy: %lg, Fz: %lg\n",
                          primsrc, tmpPot * qpr, tmpF.X * qpr, tmpF.Y * qpr,
                          tmpF.Z * qpr);
                      printf(
                          "primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: "
                          "%lg\n",
                          primsrc, pPot[primsrc], lFx, lFy, lFz);
                      fflush(stdout);
                    }
                  } else {
                    // use discretized representation of a repeated primitive
                    double tPot;
                    Vector3D tF;
                    double erPot = 0.;
                    Vector3D erF;
                    erF.X = 0.0;
                    erF.Y = 0.0;
                    erF.Z = 0.0;

                    const int eleMin = ElementBgn[primsrc];
                    const int eleMax = ElementEnd[primsrc];
                    for (int ele = eleMin; ele <= eleMax; ++ele) {
                      const double xrsrc = (EleArr + ele - 1)->G.Origin.X;
                      const double yrsrc = (EleArr + ele - 1)->G.Origin.Y;
                      const double zrsrc = (EleArr + ele - 1)->G.Origin.Z;

                      const double XEOfRpt = xrsrc + xShift;
                      const double YEOfRpt = yrsrc + yShift;
                      const double ZEOfRpt = zrsrc + zShift;

                      // Rotate point from global to local system.
                      double InitialVector[3] = {xfld - XEOfRpt, yfld - YEOfRpt, zfld - ZEOfRpt};
                      double FinalVector[3] = {0., 0., 0.};
                      for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                          FinalVector[i] +=
                              TransformationMatrix[i][j] * InitialVector[j];
                        }
                      }
                      Point3D localPER;
                      localPER.X = FinalVector[0];
                      localPER.Y = FinalVector[1];
                      localPER.Z = FinalVector[2];

                      // Allowed, because all the local coordinates have the
                      // same orientations. Only the origins are mutually
                      // displaced along a line.
                      GetPF(ele, &localPER, &tPot, &tF);
                      const double qel = WtFieldChDen[IdWtField][ele];
                      erPot += qel * tPot;
                      erF.X += qel * tF.X;
                      erF.Y += qel * tF.Y;
                      erF.Z += qel * tF.Z;
                      // if(DebugLevel == 301)
                      if (dbgFn) {
                        printf("PFAtPoint base primitive:%d\n", primsrc);
                        printf(
                            "ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n",
                            ele, localPER.X, localPER.Y, localPER.Z);
                        printf(
                            "ele: %d, tPot: %lg, tFx: %lg, tFy: %lg, tFz: %lg, "
                            "Solution: %g\n",
                            ele, tPot, tF.X, tF.Y, tF.Z, qel);
                        printf(
                            "ele: %d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: "
                            "%lg\n",
                            ele, erPot, erF.X, erF.Y, erF.Z);
                        fflush(stdout);
                      }
                    }  // for all the elements on this primsrc repeated
                       // primitive

                    pPot[primsrc] += erPot;
                    lFx += erF.X;
                    lFy += erF.Y;
                    lFz += erF.Z;
                  }  // else discretized representation of this primitive

                  // if(DebugLevel == 301)
                  if (dbgFn) {
                    printf("basic repeated xrpt: %d. yrpt: %d, zrpt: %d\n",
                           xrpt, yrpt, zrpt);
                    printf(
                        "primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: "
                        "%lg\n",
                        primsrc, pPot[primsrc], lFx, lFy, lFz);
                    fflush(stdout);
                  }
                }  // repitition of basic primitive

                if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
                    MirrorTypeZ[primsrc]) {  
                  // Mirror effect of repeated primitives - not parallelized
                  printf(
                      "Mirror not correctly implemented in this version of "
                      "neBEM ...\n");
                  exit(0);
                }        // Mirror effect for repeated primitives ends

              }  // for zrpt
            }    // for yrpt
          }      // for xrpt
        }        // PeriodicInX || PeriodicInY || PeriodicInZ
      }          // PeriodicType == 1
      Vector3D localF;
      localF.X = lFx;
      localF.Y = lFy;
      localF.Z = lFz;
      Vector3D tmpF = RotateVector3D(&localF, &PrimDC[primsrc], local2global);
      plFx[primsrc] = tmpF.X;  // local fluxes lFx, lFy, lFz in GCS
      plFy[primsrc] = tmpF.Y;
      plFz[primsrc] = tmpF.Z;
    }  // for all primitives: basic device, mirror reflections and repetitions
  }    // pragma omp parallel

  double totPot = 0.0;
  Vector3D totF;
  totF.X = totF.Y = totF.Z = 0.0;
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    totPot += pPot[prim];
    totF.X += plFx[prim];
    totF.Y += plFy[prim];
    totF.Z += plFz[prim];
  }

  // This should be done at the end of the function - before freeing memory
  *Potential = totPot;
  globalF->X = totF.X;
  globalF->Y = totF.Y;
  globalF->Z = totF.Z;

  /* for weighting field, effect of KnCh is possibly zero.
  double tmpPot; Vector3D tmpF;
  // ExactPointP and ExactPointF should also have an ExactPointPF
  // Similarly for area and volume element related functions
  // since there is no intermediate function that interfaces ExactPointP etc
  // division by MyFACTOR is necessary
  // Do parallelize before using these known charges - points or distributions
  for (int point = 1; point <= NbPtsKnCh; ++point) {
    tmpPot = ExactPointP(&(PtKnChArr+point-1)->P, globalP);
    (*Potential) += (PtKnChArr+point-1)->Assigned * tmpPot / MyFACTOR;
    ExactPointF(&(PtKnChArr+point-1)->P, globalP, &tmpF);
    globalF->X += (PtKnChArr+point-1)->Assigned * tmpF.X / MyFACTOR;
    globalF->Y += (PtKnChArr+point-1)->Assigned * tmpF.Y / MyFACTOR;
    globalF->Z += (PtKnChArr+point-1)->Assigned * tmpF.Z / MyFACTOR;
  } // for all points

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    (*Potential) += 0.0;
    globalF->X += 0.0;
    globalF->Y += 0.0;
    globalF->Z += 0.0;
  } // for all lines

  for (int area = 1; area <= NbAreasKnCh; ++area) {
    (*Potential) += 0.0;
    globalF->X += 0.0;
    globalF->Y += 0.0;
    globalF->Z += 0.0;
  } // for all areas

  for (int vol = 1; vol <= NbVolsKnCh; ++vol) {
    (*Potential) += 0.0;
    globalF->X += 0.0;
    globalF->Y += 0.0;
    globalF->Z += 0.0;
  } // for all volumes

  // This should be the final position
  // *Potential = totPot;
  // globalF->X = totF.X;
  // globalF->Y = totF.Y;
  // globalF->Z = totF.Z;
  // effect of KnCh is possibly zero on weighting field 
  */

  // (*Potential) += VSystemChargeZero;  // respect total system charge constraint

  if (dbgFn) {
    printf("Final values due to all primitives and other influences: ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n");	// refer, do not
    // uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", xfld, yfld, zfld,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  free_dvector(pPot, 1, NbPrimitives);
  free_dvector(plFx, 1, NbPrimitives);
  free_dvector(plFy, 1, NbPrimitives);
  free_dvector(plFz, 1, NbPrimitives);

  return (0);
}  // end of WtPFAtPoint

double TriLin(double xd, double yd, double zd, double c000, double c100,
              double c010, double c001, double c110, double c101, double c011,
              double c111) {
  double c00 = c000 * (1.0 - xd) + c100 * xd;
  double c10 = c010 * (1.0 - xd) + c110 * xd;
  double c01 = c001 * (1.0 - xd) + c101 * xd;
  double c11 = c011 * (1.0 - xd) + c111 * xd;
  double c0 = c00 * (1.0 - yd) + c10 * yd;
  double c1 = c01 * (1.0 - yd) + c11 * yd;
  return (c0 * (1.0 - zd) + c1 * zd);
}

#ifdef __cplusplus
} // namespace
#endif
