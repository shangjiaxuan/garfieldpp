/*
(c) 2005, Supratik Mukhopadhayay, Nayana Majumdar
*/

#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "neBEMInterface.h"
#include "Isles.h"
#include "NR.h"
#include "Vector.h"
#include "neBEM.h"

#define MyPI 3.14159265358979323846

#ifdef __cplusplus
namespace neBEM {
#endif

// GCS: global coordinate system
// PCS: primitive coordinate system
// ECS: element coordinate system

// How do we decide on the number of elements, in each direction, appropriate
// for a given surface?
// Since no linear element is being considered, no assumption is being made
// regarding a linear element for the time being.
// shall we return the pointer to the element array here? That seems to be
// more intuitive!
// Note that for the right triangle, the second vertex (vertex[1], since the
// vector begins from 0) is the 90 degree corner.
int SurfaceElements(int prim, int nvertex, double xvert[], double yvert[],
                    double zvert[], double xnorm, double ynorm, double znorm,
                    int volref1, int volref2, int inttype, double potential,
                    double charge, double lambda, int NbSegX, int NbSegZ) {
  // Decide the geometry of this primitive - if it is a rectangle, our job is
  // greatly simplified. To begin with, we check the number of vertices to
  // take the decision automatically.
  // Note that a triangle is the next best bet. All other primitive will have to
  // be reduced to rectangles (as many as possible) and triangles.
  // Incidentally, a PrimType (determined in neBEMReadGeom (neBEMInterface.c))
  // determines whether the primitive is a surface (2D) or a wire (1D),
  // for a given surface, SurfShape determines whether it is a triangle,
  // rectangle, or any other polygon besides these two (square is, of course, a
  // special case of rectangle)
  int fstatus;
  switch (nvertex) {
    case 3:  // triangle
      fstatus = DiscretizeTriangle(prim, nvertex, xvert, yvert, zvert, xnorm,
                                   ynorm, znorm, volref1, volref2, inttype,
                                   potential, charge, lambda, NbSegX, NbSegZ);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("SurfaceElements - DiscretizeTriangle");
        return -1;
      }
      break;

    case 4:  // rectangle
      fstatus = DiscretizeRectangle(prim, nvertex, xvert, yvert, zvert, xnorm,
                                    ynorm, znorm, volref1, volref2, inttype,
                                    potential, charge, lambda, NbSegX, NbSegZ);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("SurfaceElements - DiscretizeRectangle");
        return -1;
      }
      break;

    default:
      printf("nvertex out of bounds in SurfaceElements ... exiting ...\n");
      exit(-1);
  }

  return (0);
}  // end of SurfaceElements

// Analyze wire and break up into smaller wire elements
// How do we decide on the number of elements? Currently, it is based on user
// inputs that need to be made automatic and adaptive
int WireElements(int prim, int nvertex, double xvert[], double yvert[],
                 double zvert[], double radius, int volref1, int volref2,
                 int inttype, double potential, double charge, double lambda,
                 int NbSegs) {
  int fstatus;

  switch (nvertex) {
    case 2:  // wire
      fstatus =
          DiscretizeWire(prim, nvertex, xvert, yvert, zvert, radius, volref1,
                         volref2, inttype, potential, charge, lambda, NbSegs);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("WireElements - DiscretizeWire");
        return -1;
      }
      break;

    default:
      printf("nvertex out of bounds in WireElements ... exiting ...\n");
      exit(-1);
  }

  return (0);
}  // end of WireElement

// Try to set up elements on this primitive such that the average element area
// is close to what has been requested.
// If the number of elements is too small (<5), over-ride and have 5*5 elements
// If the size of the element is too small (<MINDIST2), or a side of the element
// is too small (<MINDIST), over-ride 5*5 and create elements (3*3) and (1*1)
// such that they are of acceptable size and shape (the latter, if possible)
// If the primitive is smaller than MINDIST * MINDIST, report and quit
// If the number of elements is more than a limit, make the element size
// larger than that requested. Report the incidences.
// Since, for a surface primitive, the larger number between Coord1Seg and
// Coord2Seg is used for the longer arm (look for OverSmart in
// DiscretizeTriangle and DiscretizeRectagnle), the aspect ratio problem is
// expected to be taken care of to a large extent.
int AnalyzePrimitive(int prim, int *NbSegCoord1, int *NbSegCoord2) {
  int fstatus;

  switch (NbVertices[prim]) {
    case 2:  // wire
      fstatus = AnalyzeWire(prim, NbSegCoord1);
      *NbSegCoord2 = 0;
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("AnalyzePrimitive - AnalyzeWire");
        return -1;
      }
      return (2);
      break;
    case 3:  // triangle
      fstatus = AnalyzeSurface(prim, NbSegCoord1, NbSegCoord2);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("AnalyzePrimitive - AnalyzeSurface");
        return -1;
      }
      return (3);
      break;
    case 4:  // rectangle
      fstatus = AnalyzeSurface(prim, NbSegCoord1, NbSegCoord2);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("AnalyzePrimitive - AnalyzeSurface");
        return -1;
      }
      return (4);
      break;
    default:  // no shape!
      return (0);
  }
}  // end of AnalyzePrimitive

int AnalyzeWire(int prim, int *NbSeg) {
  int nb = *NbSeg;

  if (nb < 1)  // absurd! use the trio: target, min, max
  {
    double lWire = (XVertex[prim][1] - XVertex[prim][0]) *
                       (XVertex[prim][1] - XVertex[prim][0]) +
                   (YVertex[prim][1] - YVertex[prim][0]) *
                       (YVertex[prim][1] - YVertex[prim][0]) +
                   (ZVertex[prim][1] - ZVertex[prim][0]) *
                       (ZVertex[prim][1] - ZVertex[prim][0]);
    lWire = sqrt(lWire);

    nb = (int)(lWire / ElementLengthRqstd);

    if ((nb > MinNbElementsOnLength) &&
        (nb < MaxNbElementsOnLength)) {  // nothing to be done
    }
    // Check whether the length of the wire primitive is long enough
    else if (lWire < MINDIST) {
      nb = 1;
      fprintf(fMeshLog, "Wire element too small on primitive %d!\n", prim);
    }  // if lWire < MINDIST
    // Need to have at least MinNbElementsOnLength elements per wire primitive
    else if (nb < MinNbElementsOnLength) {
      nb = MinNbElementsOnLength;
      double ellength = lWire / (double)nb;
      if (ellength <
          MINDIST)  // which may not be possible if the length is small
      {
        nb = (int)(lWire / MINDIST);
        if (nb < 1)  // However, it is necessary to have at least one element!
        {
          nb = 1;
          fprintf(fMeshLog, "Wire element very small on primitive %d!\n", prim);
        }
      }
    }  // if nb < MinNbElementsOnLength
    else if (nb > MaxNbElementsOnLength) {
      nb = MaxNbElementsOnLength;
      fprintf(fMeshLog, "Too many elements on wire primitive %d!\n", prim);
      fprintf(fMeshLog, "Number of elements reduced to maximum allowed %d\n",
              MaxNbElementsOnLength);
    }  // if nb > MaxNbElementsOnLength

    *NbSeg = nb;

    fprintf(fMeshLog, "Number of elements on wire primitive %d is %d.\n\n",
            prim, *NbSeg);
  } else {  // number of dicretization specified by user
    double lWire = (XVertex[prim][1] - XVertex[prim][0]) *
                       (XVertex[prim][1] - XVertex[prim][0]) +
                   (YVertex[prim][1] - YVertex[prim][0]) *
                       (YVertex[prim][1] - YVertex[prim][0]) +
                   (ZVertex[prim][1] - ZVertex[prim][0]) *
                       (ZVertex[prim][1] - ZVertex[prim][0]);
    lWire = sqrt(lWire);

    double ellength = lWire / (double)nb;

    if (lWire < MINDIST)  // at least one element is a necessity
    {
      nb = 1;
      fprintf(fMeshLog, "Fatal: Wire element too small on primitive %d!\n",
              prim);
    }                             // if lWire < MINDIST
    else if (ellength < MINDIST)  // element length more than twice MINDIST
    {
      nb = (int)(lWire / (2.0 * MINDIST));
      if (nb < 1) {
        nb = 1;
        fprintf(fMeshLog, "Fatal: Wire element too small on primitive %d!\n",
                prim);
      }
    }  // if ellength < MINDIST

    *NbSeg = nb;

    fprintf(fMeshLog, "Number of elements on wire primitive %d is %d.\n\n",
            prim, *NbSeg);
  }

  if (nb)
    return 0;
  else
    return -1;
}  // AnalyzeWire ends

int AnalyzeSurface(int prim, int *NbSegCoord1, int *NbSegCoord2) {
  int nb1 = *NbSegCoord1, nb2 = *NbSegCoord2;

  if ((nb1 < 1) || (nb2 < 1))  // absurd! use the trio: target, min, max
  {
    // Triangle primitives have their right angle on vertex 1
    double l1 = (XVertex[prim][0] - XVertex[prim][1]) *
                    (XVertex[prim][0] - XVertex[prim][1]) +
                (YVertex[prim][0] - YVertex[prim][1]) *
                    (YVertex[prim][0] - YVertex[prim][1]) +
                (ZVertex[prim][0] - ZVertex[prim][1]) *
                    (ZVertex[prim][0] - ZVertex[prim][1]);
    l1 = sqrt(l1);
    double l2 = (XVertex[prim][2] - XVertex[prim][1]) *
                    (XVertex[prim][2] - XVertex[prim][1]) +
                (YVertex[prim][2] - YVertex[prim][1]) *
                    (YVertex[prim][2] - YVertex[prim][1]) +
                (ZVertex[prim][2] - ZVertex[prim][1]) *
                    (ZVertex[prim][2] - ZVertex[prim][1]);
    l2 = sqrt(l2);

    // We can use the lengths independently and forget about area
    // for the time being

    // double area = l1 * l2;

    nb1 = (int)(l1 / ElementLengthRqstd);
    if ((nb1 > MinNbElementsOnLength) &&
        (nb1 < MaxNbElementsOnLength)) {  
      // nothing to be done
    } else if (l1 < MINDIST) {
      fprintf(fMeshLog, "Length1 too small on primitive %d!\n", prim);
      nb1 = 1;
    } else if (nb1 < MinNbElementsOnLength) {
      // Need to have at least MinNbElementsOnLength elements per wire primitive
      nb1 = MinNbElementsOnLength;
      double ellength = l1 / (double)nb1;
      // which may not be possible if the length is small
      if (ellength < MINDIST) {
        nb1 = (int)(l1 / MINDIST);
        // However, it is necessary to have at least one element!
        if (nb1 < 1) {
          fprintf(fMeshLog, "Length1 very small on primitive %d!\n", prim);
          nb1 = 1;
        }
      }
    }  // else if nb1 < MinNbElementsOnLength

    if (nb1 > MaxNbElementsOnLength) {
      fprintf(fMeshLog, "Too many elements on Length1 for primitive %d!\n",
              prim);
      fprintf(fMeshLog, "Number of elements reduced to maximum allowed %d\n",
              MaxNbElementsOnLength);
      nb1 = MaxNbElementsOnLength;
    }

    nb2 = (int)(l2 / ElementLengthRqstd);
    if ((nb2 > MinNbElementsOnLength) &&
        (nb2 < MaxNbElementsOnLength)) {
      // nothing to be done
    } else if (l2 < MINDIST) {
      fprintf(fMeshLog, "Length2 element too small on primitive %d!\n", prim);
      nb2 = 1;
    } else if (nb2 < MinNbElementsOnLength) {
      // Need to have at least MinNbElementsOnLength elements per wire primitive
      nb2 = MinNbElementsOnLength;
      double ellength = l2 / (double)nb2;
      // which may not be possible if the length is small
      if (ellength < MINDIST) {
        // However, it is necessary to have at least one element!
        nb2 = (int)(l2 / MINDIST);
        if (nb2 < 1) {
          fprintf(fMeshLog, "Length2 element very small on primitive %d!\n",
                  prim);
          nb2 = 1;
        }
      }
    }  // else if nb2 < MinNbElementsOnLength

    if (nb2 > MaxNbElementsOnLength) {
      fprintf(fMeshLog, "Too many elements on Length2 of primitive %d!\n",
              prim);
      fprintf(fMeshLog, "Number of elements reduced to maximum allowed %d\n",
              MaxNbElementsOnLength);
      nb2 = MaxNbElementsOnLength;
    }

    *NbSegCoord1 = nb1;
    *NbSegCoord2 = nb2;

    fprintf(fMeshLog,
            "Number of elements on surface primitive %d is %d X %d.\n\n", prim,
            *NbSegCoord1, *NbSegCoord2);
  } else {  // number of discretization specified by the user
    // Triangle primitives have their right angle on the vertex 1
    double l1 = (XVertex[prim][0] - XVertex[prim][1]) *
                    (XVertex[prim][0] - XVertex[prim][1]) +
                (YVertex[prim][0] - YVertex[prim][1]) *
                    (YVertex[prim][0] - YVertex[prim][1]) +
                (ZVertex[prim][0] - ZVertex[prim][1]) *
                    (ZVertex[prim][0] - ZVertex[prim][1]);
    l1 = sqrt(l1);
    double l2 = (XVertex[prim][2] - XVertex[prim][1]) *
                    (XVertex[prim][2] - XVertex[prim][1]) +
                (YVertex[prim][2] - YVertex[prim][1]) *
                    (YVertex[prim][2] - YVertex[prim][1]) +
                (ZVertex[prim][2] - ZVertex[prim][1]) *
                    (ZVertex[prim][2] - ZVertex[prim][1]);
    l2 = sqrt(l2);

    if (l1 > l2) {
      if (nb2 > nb1) { 
        // swap numbers to allow larger number to larger side
        int tmpnb = nb1;
        nb1 = nb2;
        nb2 = tmpnb;
      }
    }  // if l1 > l2

    double ellength1 = l1 / (double)nb1;
    double ellength2 = l2 / (double)nb2;

    if (l1 < MINDIST) {
      nb1 = 1;
      fprintf(fMeshLog, "Fatal: Side length l1 too small! prim: %d\n", prim);
    } else if (ellength1 < MINDIST)  // element length more than twice MINDIST
    {
      nb1 = (int)(l1 / (2.0 * MINDIST));
      if (nb1 < 1) {
        nb1 = 1;
        fprintf(fMeshLog, "Fatal: Side length l1 too small on primitive %d!\n",
                prim);
      }
    }  // if ellength1 < MINDIST

    if (l2 < MINDIST) {
      nb2 = 1;
      fprintf(fMeshLog, "Fatal: Side length l2 too small! prim: %d\n", prim);
    } else if (ellength2 < MINDIST)  // element length more than twice MINDIST
    {
      nb2 = (int)(l2 / (2.0 * MINDIST));
      if (nb2 < 1) {
        nb2 = 1;
        fprintf(fMeshLog, "Fatal: Side length l2 too small on primitive %d!\n",
                prim);
      }
    }  // if ellength2 < MINDIST

    *NbSegCoord1 = nb1;
    *NbSegCoord2 = nb2;

    fprintf(fMeshLog,
            "Number of elements on surface primitive %d is %d X %d.\n\n", prim,
            *NbSegCoord1, *NbSegCoord2);
  }

  if ((nb1 > 0) && (nb2 > 0))
    return 0;
  else
    return -1;
}  // AnalyzeSurface ends

// Discretize wire into linear wire elements
int DiscretizeWire(int prim, int nvertex, double xvert[], double yvert[],
                   double zvert[], double radius, int volref1, int volref2,
                   int inttype, double potential, double charge, double lambda,
                   int NbSegs) {
  int WireParentObj, WireEType;
  double WireR, WireL;
  double WireLambda, WireV;
  double WireElX, WireElY, WireElZ, WireElL;
  DirnCosn3D PrimDirnCosn;  // direction cosine of the current primitive
  char primstr[10];
  char gpElem[256];
  FILE *fPrim, *fElem, *fgpPrim, *fgpElem;

  // Check inputs
  if (PrimType[prim] != 2) {
    neBEMMessage("DiscretizeWire - PrimType in DiscretizeWire");
    return -1;
  }
  if (nvertex != 2) {
    neBEMMessage("DiscretizeWire - nvertex in DiscretizeWire");
    return -1;
  }
  if (radius < MINDIST) {
    neBEMMessage("DiscretizeWire - radius in DiscretizeWire");
    return -1;
  }
  if (NbSegs <= 0) {
    neBEMMessage("DiscretizeWire - NbSegs in DiscretizeWire");
    return -1;
  }

  if (OptPrintVertexAndNormal) {
    printf("nvertex: %d\n", nvertex);
    for (int vert = 0; vert < nvertex; ++vert) {
      printf("vert: %d, x: %lg, y: %lg, z: %lg\n", vert, xvert[vert],
             yvert[vert], zvert[vert]);
    }
    printf("radius: %lg\n", radius);
  }  // if OptPrintVertexAndNormal

  // necessary for separating filenames
  sprintf(primstr, "%d", prim);

  // in order to avoid warning messages
  fPrim = NULL;
  fElem = NULL;
  fgpElem = NULL;

  WireParentObj = 1;  // ParentObj not being used now

  WireL = sqrt((xvert[1] - xvert[0]) * (xvert[1] - xvert[0])  // length of wire
               + (yvert[1] - yvert[0]) * (yvert[1] - yvert[0]) +
               (zvert[1] - zvert[0]) * (zvert[1] - zvert[0]));
  WireR = radius;
  WireElL = WireL / NbSegs;  // length of each wire element

  // Direction cosines along the wire - note difference from surface primitives!
  // The direction along the wire is considered to be the z axis of the LCS
  // So, let us fix that axial vector first
  PrimDirnCosn.ZUnit.X = (xvert[1] - xvert[0]) / WireL;  // useful
  PrimDirnCosn.ZUnit.Y = (yvert[1] - yvert[0]) / WireL;
  PrimDirnCosn.ZUnit.Z = (zvert[1] - zvert[0]) / WireL;  // useful
  // Next, let us find out the coefficients of a plane that passes through the
  // wire centroid and is normal to the axis of the wire. This is basically the
  // mid-plane of the cylindrical wire
  // Any vector OR on the plane normal to the axial direction satisfies
  // \vec{OR} . \vec{OA} = 0
  // where O is the wire centroid, A is a point on the axis and R is a point on
  // the cylindrical surface of the wire. \vec{OA} can be easily replaced by the
  // axial vector that is equivalent to the vector PrimDirnCosn.ZUnit
  // The equation of the plane can be shown to be:
  // XCoef * X + YCoef * Y + ZCoef * Z = Const
  // double XCoef, YCoef, ZCoef, Const; - not needed any more
  // double rnorm;
  // Point3D O, R;
  // Vector3D OR;
  // O = CreatePoint3D(WireX, WireY, WireZ);
  {
    Vector3D XUnit, YUnit, ZUnit;
    ZUnit.X = PrimDirnCosn.ZUnit.X;
    ZUnit.Y = PrimDirnCosn.ZUnit.Y;
    ZUnit.Z = PrimDirnCosn.ZUnit.Z;

    /* old code	- abs instead of fabs??!!
    XCoef = ZUnit.X;
    YCoef = ZUnit.Y;
    ZCoef = ZUnit.Z;
    double WireX = 0.5 * (xvert[1] + xvert[0]);
    double WireY = 0.5 * (yvert[1] + yvert[0]);
    double WireZ = 0.5 * (zvert[1] + zvert[0]);
    Const = WireX * ZUnit.X + WireY * ZUnit.Y + WireZ * ZUnit.Z;
    if(abs(XCoef) < 1.0e-12)	// X can be anything!
            {
            XUnit.X = 1.0;
            XUnit.Y = 0.0;
            XUnit.Z = 0.0;
            YUnit = Vector3DCrossProduct(ZUnit, XUnit);
            }
    else
            {
            // For a point on the above surface where both Y and Z are zero
            O = CreatePoint3D(WireX, WireY, WireZ);
            R = CreatePoint3D(Const, 0, 0);
            // Create the vector joining O and R; find X and Y unit vectors
            OR = CreateDistanceVector3D(O,R);
            XUnit = UnitVector3D(OR);
            YUnit = Vector3DCrossProduct(ZUnit, XUnit);
            }
    old code */

    // replaced following Rob's suggestions (used functions instead of direct
    // evaluation, although the latter is probably faster)
    // x-Axis: orthogonal in the 2 largest components.
    if (fabs(ZUnit.X) >= fabs(ZUnit.Z) && fabs(ZUnit.Y) >= fabs(ZUnit.Z)) {
      XUnit.X = -ZUnit.Y;
      XUnit.Y = ZUnit.X;
      XUnit.Z = 0.0;
    } else if (fabs(ZUnit.X) >= fabs(ZUnit.Y) &&
               fabs(ZUnit.Z) >= fabs(ZUnit.Y)) {
      XUnit.X = -ZUnit.Z;
      XUnit.Y = 0.0;
      XUnit.Z = ZUnit.X;
    } else {
      XUnit.X = 0.0;
      XUnit.Y = ZUnit.Z;
      XUnit.Z = -ZUnit.Y;
    }
    XUnit = UnitVector3D(&XUnit);

    // y-Axis: vectorial product of axes 1 and 2.
    YUnit = Vector3DCrossProduct(&ZUnit, &XUnit);
    YUnit = UnitVector3D(&YUnit);
    // end of replacement

    PrimDirnCosn.XUnit.X = XUnit.X;
    PrimDirnCosn.XUnit.Y = XUnit.Y;
    PrimDirnCosn.XUnit.Z = XUnit.Z;
    PrimDirnCosn.YUnit.X = YUnit.X;
    PrimDirnCosn.YUnit.Y = YUnit.Y;
    PrimDirnCosn.YUnit.Z = YUnit.Z;
  }  // X and Y direction cosines computed

  // primitive direction cosine assignments
  PrimDC[prim].XUnit.X = PrimDirnCosn.XUnit.X;
  PrimDC[prim].XUnit.Y = PrimDirnCosn.XUnit.Y;
  PrimDC[prim].XUnit.Z = PrimDirnCosn.XUnit.Z;
  PrimDC[prim].YUnit.X = PrimDirnCosn.YUnit.X;
  PrimDC[prim].YUnit.Y = PrimDirnCosn.YUnit.Y;
  PrimDC[prim].YUnit.Z = PrimDirnCosn.YUnit.Z;
  PrimDC[prim].ZUnit.X = PrimDirnCosn.ZUnit.X;
  PrimDC[prim].ZUnit.Y = PrimDirnCosn.ZUnit.Y;
  PrimDC[prim].ZUnit.Z = PrimDirnCosn.ZUnit.Z;

  // primitive origin: also the barycenter for a wire element
  PrimOriginX[prim] = 0.5 * (xvert[0] + xvert[1]);
  PrimOriginY[prim] = 0.5 * (yvert[0] + yvert[1]);
  PrimOriginZ[prim] = 0.5 * (zvert[0] + zvert[1]);
  PrimLX[prim] = WireR;  // radius for wire
  PrimLZ[prim] = WireL;  // length of wire

  WireEType = inttype;
  WireV = potential;
  WireLambda = lambda;

  // file output for a primitive
  if (OptPrimitiveFiles) {
    char OutPrim[256];
    strcpy(OutPrim, ModelOutDir);
    strcat(OutPrim, "/Primitives/Primitive");
    strcat(OutPrim, primstr);
    strcat(OutPrim, ".out");
    fPrim = fopen(OutPrim, "w");
    if (fPrim == NULL) {
      neBEMMessage("DiscretizeWire - OutPrim");
      return -1;
    }
    fprintf(fPrim, "#prim: %d, nvertex: %d\n", prim, nvertex);
    fprintf(fPrim, "Node1: %lg\t%lg\t%lg\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fPrim, "Node2: %lg\t%lg\t%lg\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fPrim, "PrimOrigin: %lg\t%lg\t%lg\n", PrimOriginX[prim],
            PrimOriginY[prim], PrimOriginZ[prim]);
    fprintf(fPrim, "Primitive lengths: %lg\t%lg\n", PrimLX[prim], PrimLZ[prim]);
    fprintf(fPrim, "#DirnCosn: \n");
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.XUnit.X,
            PrimDirnCosn.XUnit.Y, PrimDirnCosn.XUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.YUnit.X,
            PrimDirnCosn.YUnit.Y, PrimDirnCosn.YUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.ZUnit.X,
            PrimDirnCosn.ZUnit.Y, PrimDirnCosn.ZUnit.Z);
    fprintf(fPrim, "#volref1: %d, volref2: %d\n", volref1, volref2);
    fprintf(fPrim, "#NbSegs: %d\n", NbSegs);
    fprintf(fPrim, "#ParentObj: %d\tEType: %d\n", WireParentObj, WireEType);
    fprintf(fPrim, "#WireR: %lg\tWireL: %lg\n", WireR, WireL);
    fprintf(fPrim, "#SurfLambda: %lg\tSurfV: %lg\n", WireLambda, WireV);
  }

  // necessary for gnuplot
  if (OptGnuplot && OptGnuplotPrimitives) {
    char gpPrim[256];
    strcpy(gpPrim, MeshOutDir);
    strcat(gpPrim, "/GViewDir/gpPrim");
    strcat(gpPrim, primstr);
    strcat(gpPrim, ".out");
    fgpPrim = fopen(gpPrim, "w");
    if (fgpPrim == NULL) {
      neBEMMessage("DiscretizeWire - OutgpPrim");
      return -1;
    }
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[1], yvert[1], zvert[1]);
    fclose(fgpPrim);

    if (prim == 1)
      fprintf(fgnuPrim, " '%s\' w l", gpPrim);
    else
      fprintf(fgnuPrim, ", \\\n \'%s\' w l", gpPrim);
  }

  // file outputs for elements on primitive
  if (OptElementFiles) {
    char OutElem[256];
    strcpy(OutElem, MeshOutDir);
    strcat(OutElem, "/Elements/ElemOnPrim");
    strcat(OutElem, primstr);
    strcat(OutElem, ".out");
    fElem = fopen(OutElem, "w");
    if (fElem == NULL) {
      neBEMMessage("DiscretizeWire - OutElem");
      return -1;
    }
  }
  // gnuplot friendly file outputs for elements on primitive
  if (OptGnuplot && OptGnuplotElements) {
    strcpy(gpElem, MeshOutDir);
    strcat(gpElem, "/GViewDir/gpElemOnPrim");
    strcat(gpElem, primstr);
    strcat(gpElem, ".out");
    fgpElem = fopen(gpElem, "w");
    if (fgpElem == NULL) {
      neBEMMessage("DiscretizeWire - OutgpElem");
      if (fElem) fclose(fElem);
      return -1;
    }
  }

  double xincr = (xvert[1] - xvert[0]) / (double)NbSegs;
  double yincr = (yvert[1] - yvert[0]) / (double)NbSegs;
  double zincr = (zvert[1] - zvert[0]) / (double)NbSegs;

  ElementBgn[prim] = EleCntr + 1;
  double xv0, yv0, zv0, xv1, yv1, zv1;
  for (int seg = 1; seg <= NbSegs; ++seg) {
    xv0 = xvert[0] + ((double)seg - 1.0) * xincr;
    yv0 = yvert[0] + ((double)seg - 1.0) * yincr;
    zv0 = zvert[0] + ((double)seg - 1.0) * zincr;
    xv1 = xvert[0] + ((double)seg) * xincr;
    yv1 = yvert[0] + ((double)seg) * yincr;
    zv1 = zvert[0] + ((double)seg) * zincr;
    WireElX = xvert[0] + ((double)seg - 1.0) * xincr + 0.5 * xincr;
    WireElY = yvert[0] + ((double)seg - 1.0) * yincr + 0.5 * yincr;
    WireElZ = zvert[0] + ((double)seg - 1.0) * zincr + 0.5 * zincr;

    // Assign element values and write in the file
    // If element counter exceeds the maximum allowed number of elements, warn!
    ++EleCntr;
    if (EleCntr > NbElements) {
      neBEMMessage("DiscretizeWire - EleCntr more than NbElements!");
      if (fgpElem) fclose(fgpElem);
      if (fElem) fclose(fElem);
      return -1;
    }

    (EleArr + EleCntr - 1)->DeviceNb =
        1;  // At present, there is only one device
    (EleArr + EleCntr - 1)->ComponentNb = WireParentObj;
    (EleArr + EleCntr - 1)->PrimitiveNb = prim;
    (EleArr + EleCntr - 1)->Id = EleCntr;
    (EleArr + EleCntr - 1)->G.Type = 2;  // linear (wire) here
    (EleArr + EleCntr - 1)->G.Origin.X = WireElX;
    (EleArr + EleCntr - 1)->G.Origin.Y = WireElY;
    (EleArr + EleCntr - 1)->G.Origin.Z = WireElZ;
    (EleArr + EleCntr - 1)->G.Vertex[0].X = xv0;
    (EleArr + EleCntr - 1)->G.Vertex[0].Y = yv0;
    (EleArr + EleCntr - 1)->G.Vertex[0].Z = zv0;
    (EleArr + EleCntr - 1)->G.Vertex[1].X = xv1;
    (EleArr + EleCntr - 1)->G.Vertex[1].Y = yv1;
    (EleArr + EleCntr - 1)->G.Vertex[1].Z = zv1;
    (EleArr + EleCntr - 1)->G.LX = WireR;    // radius of the wire element
    (EleArr + EleCntr - 1)->G.LZ = WireElL;  // wire element length
    (EleArr + EleCntr - 1)->G.dA = 2.0 * MyPI * (EleArr + EleCntr - 1)->G.LX *
                                   (EleArr + EleCntr - 1)->G.LZ;
    (EleArr + EleCntr - 1)->G.DC.XUnit.X = PrimDirnCosn.XUnit.X;
    (EleArr + EleCntr - 1)->G.DC.XUnit.Y = PrimDirnCosn.XUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.XUnit.Z = PrimDirnCosn.XUnit.Z;
    (EleArr + EleCntr - 1)->G.DC.YUnit.X = PrimDirnCosn.YUnit.X;
    (EleArr + EleCntr - 1)->G.DC.YUnit.Y = PrimDirnCosn.YUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.YUnit.Z = PrimDirnCosn.YUnit.Z;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.X = PrimDirnCosn.ZUnit.X;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.Y = PrimDirnCosn.ZUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.Z = PrimDirnCosn.ZUnit.Z;
    (EleArr + EleCntr - 1)->E.Type = WireEType;
    (EleArr + EleCntr - 1)->E.Lambda = WireLambda;
    (EleArr + EleCntr - 1)->Solution = 0.0;
    (EleArr + EleCntr - 1)->Assigned = charge;
    (EleArr + EleCntr - 1)->BC.NbOfBCs = 1;  // assume one BC per element
    (EleArr + EleCntr - 1)->BC.CollPt.X =
        (EleArr + EleCntr - 1)->G.Origin.X;  // modify
    (EleArr + EleCntr - 1)->BC.CollPt.Y =
        (EleArr + EleCntr - 1)->G.Origin.Y;  // to be on
    (EleArr + EleCntr - 1)->BC.CollPt.Z =
        (EleArr + EleCntr - 1)->G.Origin.Z;  // surface?

    // File operations begin
    // rfw = fwrite(&Ele, sizeof(Element), 1, fpEle);
    // printf("Return of fwrite is %d\n", rfw);

    if (OptElementFiles) {
      fprintf(fElem, "##Element Counter: %d\n", EleCntr);
      fprintf(fElem, "#DevNb\tCompNb\tPrimNb\tId\n");
      fprintf(fElem, "%d\t%d\t%d\t%d\n", (EleArr + EleCntr - 1)->DeviceNb,
              (EleArr + EleCntr - 1)->ComponentNb,
              (EleArr + EleCntr - 1)->PrimitiveNb, (EleArr + EleCntr - 1)->Id);
      fprintf(fElem, "#GType\tX\tY\tZ\tLX\tLZ\tdA\n");
      fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\n",
              (EleArr + EleCntr - 1)->G.Type,
              (EleArr + EleCntr - 1)->G.Origin.X,
              (EleArr + EleCntr - 1)->G.Origin.Y,
              (EleArr + EleCntr - 1)->G.Origin.Z, (EleArr + EleCntr - 1)->G.LX,
              (EleArr + EleCntr - 1)->G.LZ, (EleArr + EleCntr - 1)->G.dA);
      fprintf(fElem, "#DirnCosn: \n");
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.XUnit.X,
              (EleArr + EleCntr - 1)->G.DC.XUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.XUnit.Z);
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.YUnit.X,
              (EleArr + EleCntr - 1)->G.DC.YUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.YUnit.Z);
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.ZUnit.X,
              (EleArr + EleCntr - 1)->G.DC.ZUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.ZUnit.Z);
      fprintf(fElem, "#EType\tLambda\n");
      fprintf(fElem, "%d\t%lg\n", (EleArr + EleCntr - 1)->E.Type,
              (EleArr + EleCntr - 1)->E.Lambda);
      fprintf(fElem, "#NbBCs\tCPX\tCPY\tCPZ\tValue\n");
      fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%lg\n",
              (EleArr + EleCntr - 1)->BC.NbOfBCs,
              (EleArr + EleCntr - 1)->BC.CollPt.X,
              (EleArr + EleCntr - 1)->BC.CollPt.Y,
              (EleArr + EleCntr - 1)->BC.CollPt.Z,
              (EleArr + EleCntr - 1)->BC.Value);
    }  // if OptElementFiles
       //
    // mark centroid
    if (OptGnuplot && OptGnuplotElements) {
      fprintf(fgpElem, "%g\t%g\t%g\n", (EleArr + EleCntr - 1)->BC.CollPt.X,
              (EleArr + EleCntr - 1)->BC.CollPt.Y,
              (EleArr + EleCntr - 1)->BC.CollPt.Z);
    }  // if OptElementFiles
       // File operations end
  }    // seg loop for wire elements
  ElementEnd[prim] = EleCntr;
  NbElmntsOnPrim[prim] = ElementEnd[prim] - ElementBgn[prim] + 1;

  if (OptPrimitiveFiles) {
    fprintf(fPrim, "Element begin: %d, Element end: %d\n", ElementBgn[prim],
            ElementEnd[prim]);
    fprintf(fPrim, "Number of elements on primitive: %d\n",
            NbElmntsOnPrim[prim]);
    fclose(fPrim);
  }

  if (OptElementFiles) fclose(fElem);
  if (OptGnuplot && OptGnuplotElements) fclose(fgpElem);

  return (0);
}  // end of DiscretizeWire

// NbSegX is considered to be the number of rectangular + one triangular
// elements into which the base of the triangular primitive is divided.
// NbSegZ is the total number of elements into which the height is divided.
int DiscretizeTriangle(int prim, int nvertex, double xvert[], double yvert[],
                       double zvert[], double xnorm, double ynorm, double znorm,
                       int volref1, int volref2, int inttype, double potential,
                       double charge, double lambda, int NbSegX, int NbSegZ) {
  int SurfParentObj, SurfEType;
  double SurfX, SurfY, SurfZ, SurfLX, SurfLZ;
  double SurfElX, SurfElY, SurfElZ, SurfElLX, SurfElLZ;
  DirnCosn3D PrimDirnCosn;  // direction cosine of the current primitive
  char primstr[10];
  char gpElem[256], gpMesh[256];
  FILE *fPrim, *fElem, *fgpPrim, *fgpElem, *fgpMesh;

  // Check inputs
  if ((NbSegX <= 0) || (NbSegZ <= 0)) {
    printf("segmentation input wrong in DiscretizeTriangle ...\n");
    return -1;
  }

  if (OptPrintVertexAndNormal) {
    printf("nvertex: %d\n", nvertex);
    for (int vert = 0; vert < nvertex; ++vert) {
      printf("vert: %d, x: %lg, y: %lg, z: %lg\n", vert, xvert[vert],
             yvert[vert], zvert[vert]);
    }
    printf("Normal: %lg, %lg, %lg\n", xnorm, ynorm, znorm);
  }  // if OptPrintVertexAndNormal

  // necessary for separating filenames
  sprintf(primstr, "%d", prim);

  // in order to avoid warning messages
  fPrim = NULL;
  fElem = NULL;
  fgpMesh = NULL;
  fgpElem = NULL;

  // Compute all the properties of this surface
  // Boundary types from 1 to 7 have been defined
  SurfParentObj = 1;
  SurfEType = inttype;
  if ((SurfEType <= 0) || (SurfEType >= 8)) {
    printf("Wrong SurfEType for prim %d\n", prim);
    exit(-1);
  }
  // Origin of the local coordinate center is at the right angle corner
  SurfX = xvert[1];
  SurfY = yvert[1];
  SurfZ = zvert[1];

  // Find the proper direction cosines first - little more tricky that in the
  // rectangular case
  int flagDC = 0;  // DC has not been ascertained as yet
  // We begin with trial 1: one of the possible orientations
  // Intially, the lengths are necessary
  // lengths of the sides - note that right angle corner is [1]-th element
  SurfLX = sqrt((xvert[0] - xvert[1]) * (xvert[0] - xvert[1]) +
                (yvert[0] - yvert[1]) * (yvert[0] - yvert[1]) +
                (zvert[0] - zvert[1]) * (zvert[0] - zvert[1]));
  SurfLZ = sqrt((xvert[2] - xvert[1]) * (xvert[2] - xvert[1]) +
                (yvert[2] - yvert[1]) * (yvert[2] - yvert[1]) +
                (zvert[2] - zvert[1]) * (zvert[2] - zvert[1]));
  // Direction cosines - note that right angle corner is [1]-th element
  PrimDirnCosn.XUnit.X = (xvert[0] - xvert[1]) / SurfLX;
  PrimDirnCosn.XUnit.Y = (yvert[0] - yvert[1]) / SurfLX;
  PrimDirnCosn.XUnit.Z = (zvert[0] - zvert[1]) / SurfLX;
  PrimDirnCosn.ZUnit.X = (xvert[2] - xvert[1]) / SurfLZ;
  PrimDirnCosn.ZUnit.Y = (yvert[2] - yvert[1]) / SurfLZ;
  PrimDirnCosn.ZUnit.Z = (zvert[2] - zvert[1]) / SurfLZ;
  PrimDirnCosn.YUnit =
      Vector3DCrossProduct(&PrimDirnCosn.ZUnit, &PrimDirnCosn.XUnit);
  if ((fabs(PrimDirnCosn.YUnit.X - xnorm) <= 1.0e-3) &&
      (fabs(PrimDirnCosn.YUnit.Y - ynorm) <= 1.0e-3) &&
      (fabs(PrimDirnCosn.YUnit.Z - znorm) <= 1.0e-3))
    flagDC = 1;
  if (DebugLevel == 202) {
    printf("First attempt: \n");
    PrintDirnCosn3D(PrimDirnCosn);
    printf("\n");
  }

  if (!flagDC)  // if DC search is unsuccessful, try the other orientation
  {
    SurfLX = sqrt((xvert[2] - xvert[1]) * (xvert[2] - xvert[1]) +
                  (yvert[2] - yvert[1]) * (yvert[2] - yvert[1]) +
                  (zvert[2] - zvert[1]) * (zvert[2] - zvert[1]));
    SurfLZ = sqrt((xvert[0] - xvert[1]) * (xvert[0] - xvert[1]) +
                  (yvert[0] - yvert[1]) * (yvert[0] - yvert[1]) +
                  (zvert[0] - zvert[1]) * (zvert[0] - zvert[1]));
    // Direction cosines - note that right angle corner is [1]-th element
    PrimDirnCosn.XUnit.X = (xvert[2] - xvert[1]) / SurfLX;
    PrimDirnCosn.XUnit.Y = (yvert[2] - yvert[1]) / SurfLX;
    PrimDirnCosn.XUnit.Z = (zvert[2] - zvert[1]) / SurfLX;
    PrimDirnCosn.ZUnit.X = (xvert[0] - xvert[1]) / SurfLZ;
    PrimDirnCosn.ZUnit.Y = (yvert[0] - yvert[1]) / SurfLZ;
    PrimDirnCosn.ZUnit.Z = (zvert[0] - zvert[1]) / SurfLZ;
    PrimDirnCosn.YUnit =
        Vector3DCrossProduct(&PrimDirnCosn.ZUnit, &PrimDirnCosn.XUnit);
    if ((fabs(PrimDirnCosn.YUnit.X - xnorm) <= 1.0e-3) &&
        (fabs(PrimDirnCosn.YUnit.Y - ynorm) <= 1.0e-3) &&
        (fabs(PrimDirnCosn.YUnit.Z - znorm) <= 1.0e-3))
      flagDC = 2;
    if (DebugLevel == 202) {
      printf("Second attempt: \n");
      PrintDirnCosn3D(PrimDirnCosn);
      printf("\n");
    }
  }

  if (!flagDC)  // No other possibility, DC search failed!!!
  {
    printf("Triangle DC problem ... returning ...\n");
    // getchar();
    return -1;
  }

  // primitive direction cosine assignments
  PrimDC[prim].XUnit.X = PrimDirnCosn.XUnit.X;
  PrimDC[prim].XUnit.Y = PrimDirnCosn.XUnit.Y;
  PrimDC[prim].XUnit.Z = PrimDirnCosn.XUnit.Z;
  PrimDC[prim].YUnit.X = PrimDirnCosn.YUnit.X;
  PrimDC[prim].YUnit.Y = PrimDirnCosn.YUnit.Y;
  PrimDC[prim].YUnit.Z = PrimDirnCosn.YUnit.Z;
  PrimDC[prim].ZUnit.X = PrimDirnCosn.ZUnit.X;
  PrimDC[prim].ZUnit.Y = PrimDirnCosn.ZUnit.Y;
  PrimDC[prim].ZUnit.Z = PrimDirnCosn.ZUnit.Z;

  // primitive origin - for a triangle, origin is at the right corner
  PrimOriginX[prim] = SurfX;
  PrimOriginY[prim] = SurfY;
  PrimOriginZ[prim] = SurfZ;
  PrimLX[prim] = SurfLX;
  PrimLZ[prim] = SurfLZ;
  if (flagDC == 1) {
    int tmpVolRef1 = VolRef1[prim];
    VolRef1[prim] = VolRef2[prim];
    VolRef2[prim] = tmpVolRef1;
    int tmpEpsilon1 = Epsilon1[prim];
    Epsilon1[prim] = Epsilon2[prim];
    Epsilon2[prim] = tmpEpsilon1;
  }

  double SurfLambda = lambda;
  double SurfV = potential;

  // file output for a primitive
  if (OptPrimitiveFiles) {
    char OutPrim[256];
    strcpy(OutPrim, ModelOutDir);
    strcat(OutPrim, "/Primitives/Primitive");
    strcat(OutPrim, primstr);
    strcat(OutPrim, ".out");
    fPrim = fopen(OutPrim, "w");
    if (fPrim == NULL) {
      neBEMMessage("DiscretizeTriangle - OutPrim");
      return -1;
    }
    fprintf(fPrim, "#prim: %d, nvertex: %d\n", prim, nvertex);
    fprintf(fPrim, "Node1: %lg\t%lg\t%lg\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fPrim, "Node2: %lg\t%lg\t%lg\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fPrim, "Node3: %lg\t%lg\t%lg\n", xvert[2], yvert[2], zvert[2]);
    fprintf(fPrim, "PrimOrigin: %lg\t%lg\t%lg\n", PrimOriginX[prim],
            PrimOriginY[prim], PrimOriginZ[prim]);
    fprintf(fPrim, "Primitive lengths: %lg\t%lg\n", PrimLX[prim], PrimLZ[prim]);
    fprintf(fPrim, "Norm: %lg\t%lg\t%lg\n", xnorm, ynorm, znorm);
    fprintf(fPrim, "#volref1: %d, volref2: %d\n", volref1, volref2);
    fprintf(fPrim, "#NbSegX: %d, NbSegZ: %d (check note!)\n", NbSegX, NbSegZ);
    fprintf(fPrim, "#ParentObj: %d\tEType: %d\n", SurfParentObj, SurfEType);
    fprintf(fPrim, "#SurfX\tSurfY\tSurfZ\tSurfLX\tSurfLZ (Rt. Corner)\n");
    fprintf(fPrim, "%lg\t%lg\t%lg\t%lg\t%lg\n", SurfX, SurfY, SurfZ, SurfLX,
            SurfLZ);
    fprintf(fPrim, "#DirnCosn: \n");
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.XUnit.X,
            PrimDirnCosn.XUnit.Y, PrimDirnCosn.XUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.YUnit.X,
            PrimDirnCosn.YUnit.Y, PrimDirnCosn.YUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.ZUnit.X,
            PrimDirnCosn.ZUnit.Y, PrimDirnCosn.ZUnit.Z);
    fprintf(fPrim, "#SurfLambda: %lg\tSurfV: %lg\n", SurfLambda, SurfV);
  }

  // necessary for gnuplot
  if (OptGnuplot && OptGnuplotPrimitives) {
    char gpPrim[256];
    strcpy(gpPrim, MeshOutDir);
    strcat(gpPrim, "/GViewDir/gpPrim");
    strcat(gpPrim, primstr);
    strcat(gpPrim, ".out");
    fgpPrim = fopen(gpPrim, "w");
    if (fgpPrim == NULL) {
      neBEMMessage("DiscretizeTriangle - OutgpPrim");
      return -1;
    }
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[2], yvert[2], zvert[2]);
    fprintf(fgpPrim, "%g\t%g\t%g\n", xvert[0], yvert[0], zvert[0]);
    fclose(fgpPrim);

    if (prim == 1)
      fprintf(fgnuPrim, " '%s\' w l", gpPrim);
    else
      fprintf(fgnuPrim, ", \\\n \'%s\' w l", gpPrim);
  }

  // file outputs for elements on primitive
  if (OptElementFiles) {
    char OutElem[256];
    strcpy(OutElem, MeshOutDir);
    strcat(OutElem, "/Elements/ElemOnPrim");
    strcat(OutElem, primstr);
    strcat(OutElem, ".out");
    fElem = fopen(OutElem, "w");
    if (fElem == NULL) {
      neBEMMessage("DiscretizeTriangle - OutElem");
      return -1;
    }
  }
  // gnuplot friendly file outputs for elements on primitive
  if (OptGnuplot && OptGnuplotElements) {
    strcpy(gpElem, MeshOutDir);
    strcat(gpElem, "/GViewDir/gpElemOnPrim");
    strcat(gpElem, primstr);
    strcat(gpElem, ".out");
    fgpElem = fopen(gpElem, "w");
    // assert(fgpElem != NULL);
    if (fgpElem == NULL) {
      neBEMMessage("DiscretizeTriangle - OutgpElem");
      if (fElem) fclose(fElem);
      return -1;
    }
    // gnuplot friendly file outputs for elements on primitive
    strcpy(gpMesh, MeshOutDir);
    strcat(gpMesh, "/GViewDir/gpMeshOnPrim");
    strcat(gpMesh, primstr);
    strcat(gpMesh, ".out");
    fgpMesh = fopen(gpMesh, "w");
    if (fgpMesh == NULL) {
      neBEMMessage("DiscretizeTriangle - OutgpMesh");
      fclose(fgpElem);
      if (fElem) fclose(fElem);
      return -1;
    }
  }

  // Compute element positions (CGs) in primitive local coordinate system (PCS).
  // Then map these CGs to the global coordinate system.
  // (xav, 0, zav) is the CG of an element wrt the primitive coordinate system
  // From this, we find the offset of the element from the primitive CG in the
  // global coordinate system by just rotating the displacement vector
  // (0,0,0) -> (xav,0,zav) using the known DCM of the surface.
  // This rotated vector (now in the global system) is then added to the
  // position vector of the primitive CG (always in the global system) to get
  // the CG of the element in the global system.

  // Consult note for clarifications on the discretization algorithm.
  // SurfElLX is true for the lowest row  of elements, but indicative of the
  // other rows
  // Make sure that the larger number is being used to discretize the longer
  // side - can be OverSmart in certain cases - there may be cases where
  // smaller number of elements suffice for a longer side
  if (NbSegX == NbSegZ) {
    SurfElLX = SurfLX / NbSegX;  // element sizes
    SurfElLZ = SurfLZ / NbSegZ;
  } else if (NbSegX > NbSegZ) {
    if (SurfLX > SurfLZ) {
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    } else  // interchange NbSegX and NbSegZ
    {
      int tmp = NbSegZ;
      NbSegZ = NbSegX;
      NbSegX = tmp;
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    }
  }     // NbSegX > NbSegZ
  else  // NbSegX < NbSegZ
  {
    if (SurfLX < SurfLZ) {
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    } else  // interchange NbSegX and NbSegZ
    {
      int tmp = NbSegZ;
      NbSegZ = NbSegX;
      NbSegX = tmp;
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    }
  }

  // Analysis of element aspect ratio.
  // Note that we can afford only to reduce the total number of elements.
  // Otherwise, we'll have to realloc `EleArr' array.
  // Later, we'll make the corrections such that the total number of elements
  // remain close to the originally intended value.
  double AR = SurfElLX / SurfElLZ;  // indicative element aspect ratio
  if (OptPrimitiveFiles) {
    fprintf(fPrim,
            "Using the input, the aspect ratio of the elements on prim: %d\n",
            prim);
    fprintf(fPrim,
            "NbSegX: %d, SurfElLX: %lg, NbSegZ: %d, SurfElLZ: %lg, AR: %lg\n",
            NbSegX, SurfElLX, NbSegZ, SurfElLZ, AR);
  }
  if (AR > 10.0)  // NbSegZ needs to be reduced
  {
    double tmpElLZ = SurfElLX / 10.0;
    NbSegZ = (int)(SurfLZ / tmpElLZ);
    if (NbSegZ <= 0) NbSegZ = 1;
    SurfElLZ = SurfLZ / NbSegZ;
  }
  if (AR < 0.1)  // NbSegX need to be reduced
  {
    double tmpElLX = SurfElLZ * 0.1;
    NbSegX = (int)(SurfLX / tmpElLX);
    if (NbSegX <= 0) NbSegX = 1;
    SurfElLX = SurfLX / NbSegX;
  }
  AR = SurfElLX / SurfElLZ;  // indicative element aspect ratio
  if (OptPrimitiveFiles) {
    fprintf(fPrim, "After analyzing the likely aspect ratio of the elements\n");
    fprintf(fPrim,
            "NbSegX: %d, SurfElLX: %lg, NbSegZ: %d, SurfElLZ: %lg, AR: %lg\n",
            NbSegX, SurfElLX, NbSegZ, SurfElLZ, AR);
  }

  // The top-most right angle triangle is a lone element on the highest row, its
  // properties being determined irrespective of the others. Despite that, this
  // element is being treated as one among the others, especially to facilitate
  // element indexing, file output etc.
  // Each row, thus, has at least one triangle, and possibly several rectangular
  // elements. All the elements in any given row has the same SurfElLZ as has
  // been determined above.
  // First we create the triangular element and then move on to create the
  // rectangular ones.
  double xv0, yv0, zv0, xv1, yv1, zv1, xv2, yv2, zv2;
  ElementBgn[prim] = EleCntr + 1;
  for (int k = 1; k <= NbSegZ; ++k)  // consider the k-th row
  {
    double grad = (SurfLZ / SurfLX);
    double zlopt = (k - 1) * SurfElLZ;
    double zhipt = (k)*SurfElLZ;
    double xlopt = (SurfLZ - zlopt) / grad;  // used again on 21 Feb 2014
    double xhipt = (SurfLZ - zhipt) / grad;

    // the triangular element on the k-th row can now be specified in PCS
    double xtorigin = xhipt;
    double ytorigin = 0.0;
    double ztorigin = zlopt;
    {  // Separate block for position rotation - local2global
      Point3D localDisp, globalDisp;

      localDisp.X = xtorigin;
      localDisp.Y = ytorigin;
      localDisp.Z = ztorigin;
      globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
      SurfElX =
          SurfX + globalDisp.X;  // these are the coords in GCS of origin of
      SurfElY =
          SurfY + globalDisp.Y;  // the triangluar element under consideration
      SurfElZ = SurfZ + globalDisp.Z;
    }  // vector rotation over

    // Assign element values and write in the file
    ++EleCntr;
    if (EleCntr > NbElements) {
      neBEMMessage("DiscretizeTriangle - EleCntr more than NbElements 1!");
      if (fgpElem) fclose(fgpElem);
      if (fgpMesh) fclose(fgpMesh);
      return -1;
    }

    (EleArr + EleCntr - 1)->DeviceNb =
        1;  // At present, there is only one device
    (EleArr + EleCntr - 1)->ComponentNb = SurfParentObj;
    (EleArr + EleCntr - 1)->PrimitiveNb = prim;
    (EleArr + EleCntr - 1)->Id = EleCntr;
    (EleArr + EleCntr - 1)->G.Type = 3;  // triangular here
    (EleArr + EleCntr - 1)->G.Origin.X = SurfElX;
    (EleArr + EleCntr - 1)->G.Origin.Y = SurfElY;
    (EleArr + EleCntr - 1)->G.Origin.Z = SurfElZ;
    // (EleArr+EleCntr-1)->G.LX = SurfElLX;	// previously written as xlopt -
    // xhipt;
    (EleArr + EleCntr - 1)->G.LX =
        xlopt - xhipt;  // back to old ways on 21 Feb 2014
    (EleArr + EleCntr - 1)->G.LZ = SurfElLZ;
    (EleArr + EleCntr - 1)->G.LZ =
        zhipt - zlopt;  // to be on the safe side, 21/2/14
    (EleArr + EleCntr - 1)->G.dA =
        0.5 * (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.LZ;
    // Safe to use the direction cosines obtained for the triangular primitive
    // since they are bound to remain unchanged for the rectangular sub-elements
    (EleArr + EleCntr - 1)->G.DC.XUnit.X = PrimDirnCosn.XUnit.X;
    (EleArr + EleCntr - 1)->G.DC.XUnit.Y = PrimDirnCosn.XUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.XUnit.Z = PrimDirnCosn.XUnit.Z;
    (EleArr + EleCntr - 1)->G.DC.YUnit.X = PrimDirnCosn.YUnit.X;
    (EleArr + EleCntr - 1)->G.DC.YUnit.Y = PrimDirnCosn.YUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.YUnit.Z = PrimDirnCosn.YUnit.Z;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.X = PrimDirnCosn.ZUnit.X;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.Y = PrimDirnCosn.ZUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.Z = PrimDirnCosn.ZUnit.Z;
    (EleArr + EleCntr - 1)->E.Type = SurfEType;
    (EleArr + EleCntr - 1)->E.Lambda = SurfLambda;
    (EleArr + EleCntr - 1)->Solution = 0.0;
    (EleArr + EleCntr - 1)->Assigned = charge;
    // Boundary condition is applied at the barycenter, not at the origin
    // of the element coordinate system (ECS) which is at the right corner
    (EleArr + EleCntr - 1)->BC.NbOfBCs = 1;  // assume one BC per element

    xv0 = (EleArr + EleCntr - 1)->G.Origin.X;
    yv0 = (EleArr + EleCntr - 1)->G.Origin.Y;
    zv0 = (EleArr + EleCntr - 1)->G.Origin.Z;
    xv1 = (EleArr + EleCntr - 1)->G.Origin.X +
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.DC.XUnit.X;
    yv1 = (EleArr + EleCntr - 1)->G.Origin.Y +
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.DC.XUnit.Y;
    zv1 = (EleArr + EleCntr - 1)->G.Origin.Z +
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.DC.XUnit.Z;
    xv2 = (EleArr + EleCntr - 1)->G.Origin.X +
          (EleArr + EleCntr - 1)->G.LZ * (EleArr + EleCntr - 1)->G.DC.ZUnit.X;
    yv2 = (EleArr + EleCntr - 1)->G.Origin.Y +
          (EleArr + EleCntr - 1)->G.LZ * (EleArr + EleCntr - 1)->G.DC.ZUnit.Y;
    zv2 = (EleArr + EleCntr - 1)->G.Origin.Z +
          (EleArr + EleCntr - 1)->G.LZ * (EleArr + EleCntr - 1)->G.DC.ZUnit.Z;
    // assign vertices of the element
    (EleArr + EleCntr - 1)->G.Vertex[0].X = xv0;
    (EleArr + EleCntr - 1)->G.Vertex[0].Y = yv0;
    (EleArr + EleCntr - 1)->G.Vertex[0].Z = zv0;
    (EleArr + EleCntr - 1)->G.Vertex[1].X = xv1;
    (EleArr + EleCntr - 1)->G.Vertex[1].Y = yv1;
    (EleArr + EleCntr - 1)->G.Vertex[1].Z = zv1;
    (EleArr + EleCntr - 1)->G.Vertex[2].X = xv2;
    (EleArr + EleCntr - 1)->G.Vertex[2].Y = yv2;
    (EleArr + EleCntr - 1)->G.Vertex[2].Z = zv2;
    (EleArr + EleCntr - 1)->G.Vertex[3].X = 0.0;
    (EleArr + EleCntr - 1)->G.Vertex[3].Y = 0.0;
    (EleArr + EleCntr - 1)->G.Vertex[3].Z = 0.0;

    if (DebugLevel == 201) {
      printf("Primitive nb: %d\n", (EleArr + EleCntr - 1)->PrimitiveNb);
      printf("Element id: %d\n", (EleArr + EleCntr - 1)->Id);
      printf("Element X, Y, Z: %lg %lg %lg\n",
             (EleArr + EleCntr - 1)->G.Origin.X,
             (EleArr + EleCntr - 1)->G.Origin.Y,
             (EleArr + EleCntr - 1)->G.Origin.Z);
      printf("Element LX, LZ: %lg %lg\n", (EleArr + EleCntr - 1)->G.LX,
             (EleArr + EleCntr - 1)->G.LZ);
      printf("Element (primitive) X axis dirn cosines: %lg, %lg, %lg\n",
             PrimDirnCosn.XUnit.X, PrimDirnCosn.XUnit.Y, PrimDirnCosn.XUnit.Z);
      printf("Element (primitive) Y axis dirn cosines: %lg, %lg, %lg\n",
             PrimDirnCosn.YUnit.X, PrimDirnCosn.YUnit.Y, PrimDirnCosn.YUnit.Z);
      printf("Element (primitive) Z axis dirn cosines: %lg, %lg, %lg\n",
             PrimDirnCosn.ZUnit.X, PrimDirnCosn.ZUnit.Y, PrimDirnCosn.ZUnit.Z);
    }
    // Following are the location in the ECS
    double dxl = (EleArr + EleCntr - 1)->G.LX / 3.0;
    double dyl = 0.0;
    double dzl = (EleArr + EleCntr - 1)->G.LZ / 3.0;
    {  // Separate block for position rotation - local2global
      Point3D localDisp, globalDisp;

      localDisp.X = dxl;
      localDisp.Y = dyl;
      localDisp.Z = dzl;
      if (DebugLevel == 201) {
        printf("Element dxl, dxy, dxz: %lg %lg %lg\n", localDisp.X, localDisp.Y,
               localDisp.Z);
      }

      globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
      (EleArr + EleCntr - 1)->BC.CollPt.X =
          (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
      (EleArr + EleCntr - 1)->BC.CollPt.Y =
          (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
      (EleArr + EleCntr - 1)->BC.CollPt.Z =
          (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      if (DebugLevel == 201) {
        printf("Element global dxl, dxy, dxz: %lg %lg %lg\n", globalDisp.X,
               globalDisp.Y, globalDisp.Z);
        printf("Element BCX, BCY, BCZ: %lg %lg %lg\n",
               (EleArr + EleCntr - 1)->BC.CollPt.X,
               (EleArr + EleCntr - 1)->BC.CollPt.Y,
               (EleArr + EleCntr - 1)->BC.CollPt.Z);
      }
    }  // vector rotation over
    // (EleArr+EleCntr-1)->BC.Value = SurfV; // assigned in BoundaryConditions

    if (OptElementFiles) {
      fprintf(fElem, "##Element Counter: %d\n", EleCntr);
      fprintf(fElem, "#DevNb\tCompNb\tPrimNb\tId\n");
      fprintf(fElem, "%d\t%d\t%d\t%d\n", (EleArr + EleCntr - 1)->DeviceNb,
              (EleArr + EleCntr - 1)->ComponentNb,
              (EleArr + EleCntr - 1)->PrimitiveNb, (EleArr + EleCntr - 1)->Id);
      fprintf(fElem, "#GType\tX\tY\tZ\tLX\tLZ\tdA\n");
      fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\n",
              (EleArr + EleCntr - 1)->G.Type,
              (EleArr + EleCntr - 1)->G.Origin.X,
              (EleArr + EleCntr - 1)->G.Origin.Y,
              (EleArr + EleCntr - 1)->G.Origin.Z, (EleArr + EleCntr - 1)->G.LX,
              (EleArr + EleCntr - 1)->G.LZ, (EleArr + EleCntr - 1)->G.dA);
      fprintf(fElem, "#DirnCosn: \n");
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.XUnit.X,
              (EleArr + EleCntr - 1)->G.DC.XUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.XUnit.Z);
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.YUnit.X,
              (EleArr + EleCntr - 1)->G.DC.YUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.YUnit.Z);
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.ZUnit.X,
              (EleArr + EleCntr - 1)->G.DC.ZUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.ZUnit.Z);
      fprintf(fElem, "#EType\tLambda\n");
      fprintf(fElem, "%d\t%lg\n", (EleArr + EleCntr - 1)->E.Type,
              (EleArr + EleCntr - 1)->E.Lambda);
      fprintf(fElem, "#NbBCs\tCPX\tCPY\tCPZ\tValue\n");
      fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%lg\n",
              (EleArr + EleCntr - 1)->BC.NbOfBCs,
              (EleArr + EleCntr - 1)->BC.CollPt.X,
              (EleArr + EleCntr - 1)->BC.CollPt.Y,
              (EleArr + EleCntr - 1)->BC.CollPt.Z,
              (EleArr + EleCntr - 1)->BC.Value);
    }  // if OptElementFiles

    // mark bary-center and draw mesh
    if (OptGnuplot && OptGnuplotElements) {
      fprintf(fgpElem, "%g\t%g\t%g\n", (EleArr + EleCntr - 1)->BC.CollPt.X,
              (EleArr + EleCntr - 1)->BC.CollPt.Y,
              (EleArr + EleCntr - 1)->BC.CollPt.Z);

      // draw mesh
      // assign vertices of the element
      xv0 = (EleArr + EleCntr - 1)->G.Vertex[0].X;
      yv0 = (EleArr + EleCntr - 1)->G.Vertex[0].Y;
      zv0 = (EleArr + EleCntr - 1)->G.Vertex[0].Z;
      xv1 = (EleArr + EleCntr - 1)->G.Vertex[1].X;
      yv1 = (EleArr + EleCntr - 1)->G.Vertex[1].Y;
      zv1 = (EleArr + EleCntr - 1)->G.Vertex[1].Z;
      xv2 = (EleArr + EleCntr - 1)->G.Vertex[2].X;
      yv2 = (EleArr + EleCntr - 1)->G.Vertex[2].Y;
      zv2 = (EleArr + EleCntr - 1)->G.Vertex[2].Z;

      fprintf(fgpMesh, "%g\t%g\t%g\n", xv0, yv0, zv0);
      fprintf(fgpMesh, "%g\t%g\t%g\n", xv1, yv1, zv1);
      fprintf(fgpMesh, "%g\t%g\t%g\n", xv2, yv2, zv2);
      fprintf(fgpMesh, "%g\t%g\t%g\n\n", xv0, yv0, zv0);
    }  // if OptGnuplotElements

    if (k == NbSegZ)  // no rectangular element on this row
      continue;

    // determine NbSegXOnThisRow and ElLXOnThisRow for the rectagnular elements
    // and then loop.
    double RowLX = xhipt;  // the triangular portion is left outside the slice
    int NbSegXOnThisRow;
    double ElLXOnThisRow;
    if (RowLX <= SurfElLX) {
      NbSegXOnThisRow = 1;
      ElLXOnThisRow = RowLX;
    } else {
      NbSegXOnThisRow = (int)(RowLX / SurfElLX);
      ElLXOnThisRow = RowLX / NbSegXOnThisRow;
    }
    double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
    for (int i = 1; i <= NbSegXOnThisRow; ++i) {
      double xorigin = (i - 1) * ElLXOnThisRow + 0.5 * ElLXOnThisRow;  // PCS
      double yorigin = 0.0;  // centroid of the rectagnular element
      double zorigin = 0.5 * (zlopt + zhipt);
      // printf("k: %d, i: %d, xo: %lg, yo: %lg, zo: %lg\n", k, i,
      // xorigin, yorigin, zorigin);
      // printf("xlopt: %lg, zlopt: %lg, xhipt: %lg, zhipt: %lg\n",
      // xlopt, zlopt, xhipt, zhipt);

      {  // Separate block for vector rotation
        Point3D localDisp, globalDisp;

        localDisp.X = xorigin;
        localDisp.Y = yorigin;
        localDisp.Z = zorigin;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        SurfElX = SurfX + globalDisp.X;  // GCS
        SurfElY = SurfY + globalDisp.Y;
        SurfElZ = SurfZ + globalDisp.Z;
      }  // vector rotation over
      // printf("SurfX: %lg, SurfY: %lg, SurfZ: %lg\n",
      // SurfX, SurfY, SurfZ);
      // printf("SurfElX: %lg, SurfElY: %lg, SurfElZ: %lg\n",
      // SurfElX, SurfElY, SurfElZ);

      // Assign element values and write in the file
      ++EleCntr;
      if (EleCntr > NbElements) {
        neBEMMessage("DiscretizeTriangle - EleCntr more than NbElements 2!");
        return -1;
      }

      (EleArr + EleCntr - 1)->DeviceNb =
          1;  // At present, there is only one device
      (EleArr + EleCntr - 1)->ComponentNb = SurfParentObj;
      (EleArr + EleCntr - 1)->PrimitiveNb = prim;
      (EleArr + EleCntr - 1)->Id = EleCntr;
      (EleArr + EleCntr - 1)->G.Type = 4;  // rectagnular here
      (EleArr + EleCntr - 1)->G.Origin.X = SurfElX;
      (EleArr + EleCntr - 1)->G.Origin.Y = SurfElY;
      (EleArr + EleCntr - 1)->G.Origin.Z = SurfElZ;
      (EleArr + EleCntr - 1)->G.LX = ElLXOnThisRow;
      (EleArr + EleCntr - 1)->G.LZ = SurfElLZ;
      (EleArr + EleCntr - 1)->G.LZ =
          zhipt - zlopt;  // to be on the safe side! 21/2/14
      (EleArr + EleCntr - 1)->G.dA =
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.LZ;
      // Safe to use the direction cosines obtained for the triangular primitive
      // since they are bound to remain unchanged for the triangular
      // sub-elements
      (EleArr + EleCntr - 1)->G.DC.XUnit.X = PrimDirnCosn.XUnit.X;
      (EleArr + EleCntr - 1)->G.DC.XUnit.Y = PrimDirnCosn.XUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.XUnit.Z = PrimDirnCosn.XUnit.Z;
      (EleArr + EleCntr - 1)->G.DC.YUnit.X = PrimDirnCosn.YUnit.X;
      (EleArr + EleCntr - 1)->G.DC.YUnit.Y = PrimDirnCosn.YUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.YUnit.Z = PrimDirnCosn.YUnit.Z;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.X = PrimDirnCosn.ZUnit.X;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.Y = PrimDirnCosn.ZUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.Z = PrimDirnCosn.ZUnit.Z;
      (EleArr + EleCntr - 1)->E.Type = SurfEType;
      (EleArr + EleCntr - 1)->E.Lambda = SurfLambda;
      (EleArr + EleCntr - 1)->Solution = 0.0;
      (EleArr + EleCntr - 1)->Assigned = charge;
      // Boundary condition is applied at the origin for this rectangular
      // element coordinate system (ECS)
      (EleArr + EleCntr - 1)->BC.NbOfBCs = 1;  // assume one BC per element
      // Following are the location in the ECS
      (EleArr + EleCntr - 1)->BC.CollPt.X = (EleArr + EleCntr - 1)->G.Origin.X;
      (EleArr + EleCntr - 1)->BC.CollPt.Y = (EleArr + EleCntr - 1)->G.Origin.Y;
      (EleArr + EleCntr - 1)->BC.CollPt.Z = (EleArr + EleCntr - 1)->G.Origin.Z;
      // find element vertices
      // 1) displacement vector in the ECS is first identified
      // 2) this vector is transformed to the GCS
      // 3) the global displacement vector, when added to the centroid in GCS,
      // gives the node positions in GCS.
      x0 = -0.5 *
           (EleArr + EleCntr - 1)->G.LX;  // xyz displacement of node wrt ECS
      y0 = 0.0;
      z0 = -0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x0;
        localDisp.Y = y0;
        localDisp.Z = z0;  // displacement in GCS
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x0 = (EleArr + EleCntr - 1)->G.Origin.X +
             globalDisp.X;  // xyz position in GCS
        y0 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z0 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x1 = 0.5 * (EleArr + EleCntr - 1)->G.LX;
      y1 = 0.0;
      z1 = -0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x1;
        localDisp.Y = y1;
        localDisp.Z = z1;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x1 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y1 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z1 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x2 = 0.5 * (EleArr + EleCntr - 1)->G.LX;
      y2 = 0.0;
      z2 = 0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x2;
        localDisp.Y = y2;
        localDisp.Z = z2;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x2 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y2 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z2 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x3 = -0.5 * (EleArr + EleCntr - 1)->G.LX;
      y3 = 0.0;
      z3 = 0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x3;
        localDisp.Y = y3;
        localDisp.Z = z3;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x3 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y3 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z3 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      // assign vertices of the element
      (EleArr + EleCntr - 1)->G.Vertex[0].X = x0;
      (EleArr + EleCntr - 1)->G.Vertex[0].Y = y0;
      (EleArr + EleCntr - 1)->G.Vertex[0].Z = z0;
      (EleArr + EleCntr - 1)->G.Vertex[1].X = x1;
      (EleArr + EleCntr - 1)->G.Vertex[1].Y = y1;
      (EleArr + EleCntr - 1)->G.Vertex[1].Z = z1;
      (EleArr + EleCntr - 1)->G.Vertex[2].X = x2;
      (EleArr + EleCntr - 1)->G.Vertex[2].Y = y2;
      (EleArr + EleCntr - 1)->G.Vertex[2].Z = z2;
      (EleArr + EleCntr - 1)->G.Vertex[3].X = x3;
      (EleArr + EleCntr - 1)->G.Vertex[3].Y = y3;
      (EleArr + EleCntr - 1)->G.Vertex[3].Z = z3;

      if (OptElementFiles) {
        fprintf(fElem, "##Element Counter: %d\n", EleCntr);
        fprintf(fElem, "#DevNb\tCompNb\tPrimNb\tId\n");
        fprintf(fElem, "%d\t%d\t%d\t%d\n", (EleArr + EleCntr - 1)->DeviceNb,
                (EleArr + EleCntr - 1)->ComponentNb,
                (EleArr + EleCntr - 1)->PrimitiveNb,
                (EleArr + EleCntr - 1)->Id);
        fprintf(fElem, "#GType\tX\tY\tZ\tLX\tLZ\tdA\n");
        fprintf(
            fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\n",
            (EleArr + EleCntr - 1)->G.Type, (EleArr + EleCntr - 1)->G.Origin.X,
            (EleArr + EleCntr - 1)->G.Origin.Y,
            (EleArr + EleCntr - 1)->G.Origin.Z, (EleArr + EleCntr - 1)->G.LX,
            (EleArr + EleCntr - 1)->G.LZ, (EleArr + EleCntr - 1)->G.dA);
        fprintf(fElem, "#DirnCosn: \n");
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.XUnit.X,
                (EleArr + EleCntr - 1)->G.DC.XUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.XUnit.Z);
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.YUnit.X,
                (EleArr + EleCntr - 1)->G.DC.YUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.YUnit.Z);
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.ZUnit.X,
                (EleArr + EleCntr - 1)->G.DC.ZUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.ZUnit.Z);
        fprintf(fElem, "#EType\tLambda\n");
        fprintf(fElem, "%d\t%lg\n", (EleArr + EleCntr - 1)->E.Type,
                (EleArr + EleCntr - 1)->E.Lambda);
        fprintf(fElem, "#NbBCs\tCPX\tCPY\tCPZ\tValue\n");
        fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%lg\n",
                (EleArr + EleCntr - 1)->BC.NbOfBCs,
                (EleArr + EleCntr - 1)->BC.CollPt.X,
                (EleArr + EleCntr - 1)->BC.CollPt.Y,
                (EleArr + EleCntr - 1)->BC.CollPt.Z,
                (EleArr + EleCntr - 1)->BC.Value);
      }  // if OptElementFiles

      // draw centroid and mesh
      if (OptGnuplot && OptGnuplotElements) {
        fprintf(fgpElem, "%g\t%g\t%g\n", (EleArr + EleCntr - 1)->BC.CollPt.X,
                (EleArr + EleCntr - 1)->BC.CollPt.Y,
                (EleArr + EleCntr - 1)->BC.CollPt.Z);
      }  // if OptGnuplot && OptGnuplotElements

      if (OptGnuplot && OptGnuplotElements) {
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x0, y0, z0);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x1, y1, z1);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x2, y2, z2);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x3, y3, z3);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x0, y0, z0);
      }  // if OptGnuplot && OptGnuplotElements
    }    // for i
  }      // for k
  ElementEnd[prim] = EleCntr;
  NbElmntsOnPrim[prim] = ElementEnd[prim] - ElementBgn[prim] + 1;

  if (OptPrimitiveFiles) {
    fprintf(fPrim, "Element begin: %d, Element end: %d\n", ElementBgn[prim],
            ElementEnd[prim]);
    fprintf(fPrim, "Number of elements on primitive: %d\n",
            NbElmntsOnPrim[prim]);
    fclose(fPrim);
  }

  if (OptGnuplot && OptGnuplotElements) {
    if (prim == 1)
      fprintf(fgnuElem, " '%s\' w p", gpElem);
    else
      fprintf(fgnuElem, ", \\\n \'%s\' w p", gpElem);
    if (prim == 1) {
      fprintf(fgnuMesh, " '%s\' w l", gpMesh);
      fprintf(fgnuMesh, ", \\\n \'%s\' w p ps 1", gpElem);
    } else {
      fprintf(fgnuMesh, ", \\\n \'%s\' w l", gpMesh);
      fprintf(fgnuMesh, ", \\\n \'%s\' w p ps 1", gpElem);
    }

    fclose(fgpElem);
    fclose(fgpMesh);
  }  // if OptGnuplot && OptGnuplotElements

  if (OptElementFiles) fclose(fElem);

  return (0);
}  // end of DiscretizeTriangles

// It may be noted here that the direction cosines of a given primitive are
// the same as those of the elements residing on it. There is only a chage
// of origin associated here.
int DiscretizeRectangle(int prim, int nvertex, double xvert[], double yvert[],
                        double zvert[], double xnorm, double ynorm,
                        double znorm, int volref1, int volref2, int inttype,
                        double potential, double charge, double lambda,
                        int NbSegX, int NbSegZ) {
  int SurfParentObj, SurfEType;
  double SurfX, SurfY, SurfZ, SurfLX, SurfLZ;
  double SurfElX, SurfElY, SurfElZ, SurfElLX, SurfElLZ;
  DirnCosn3D PrimDirnCosn;  // direction cosine of the current primitive
  char primstr[10];
  char gpElem[256], gpMesh[256];
  FILE *fPrim, *fElem, *fgpPrim, *fgpElem, *fgpMesh;

  // Check inputs
  if ((NbSegX <= 0) || (NbSegZ <= 0)) {
    printf("segmentation input wrong in DiscretizeRectangle ...\n");
    return -1;
  }

  if (OptPrintVertexAndNormal) {
    printf("nvertex: %d\n", nvertex);
    for (int vert = 0; vert < nvertex; ++vert) {
      printf("vert: %d, x: %lg, y: %lg, z: %lg\n", vert, xvert[vert],
             yvert[vert], zvert[vert]);
    }
    printf("Normal: %lg, %lg, %lg\n", xnorm, ynorm, znorm);
  }  // if OptPrintVertexAndNormal

  // necessary for separating filenames
  sprintf(primstr, "%d", prim);

  // in order to avoid warning messages
  fPrim = NULL;
  fElem = NULL;
  fgpMesh = NULL;
  fgpElem = NULL;

  // Get volume information, to begin with
  // int shape, material, boundarytype;
  // double eps, potential, charge;
  // neBEMVolumeDescription(volref1, &shape, &material,
  // &eps, &potential, &charge, &boundarytype);

  // compute all the properties of this surface
  SurfParentObj = 1;
  SurfEType = inttype;
  if (SurfEType == 0) {
    printf("Wrong SurfEType for prim %d\n", prim);
    exit(-1);
  }
  // centroid of the local coordinate system
  SurfX = 0.25 * (xvert[0] + xvert[1] + xvert[2] + xvert[3]);
  SurfY = 0.25 * (yvert[0] + yvert[1] + yvert[2] + yvert[3]);
  SurfZ = 0.25 * (zvert[0] + zvert[1] + zvert[2] + zvert[3]);
  // lengths of the sides
  SurfLX = sqrt((xvert[1] - xvert[0]) * (xvert[1] - xvert[0]) +
                (yvert[1] - yvert[0]) * (yvert[1] - yvert[0]) +
                (zvert[1] - zvert[0]) * (zvert[1] - zvert[0]));
  SurfLZ = sqrt((xvert[2] - xvert[1]) * (xvert[2] - xvert[1]) +
                (yvert[2] - yvert[1]) * (yvert[2] - yvert[1]) +
                (zvert[2] - zvert[1]) * (zvert[2] - zvert[1]));
  // Direction cosines
  PrimDirnCosn.XUnit.X = (xvert[1] - xvert[0]) / SurfLX;
  PrimDirnCosn.XUnit.Y = (yvert[1] - yvert[0]) / SurfLX;
  PrimDirnCosn.XUnit.Z = (zvert[1] - zvert[0]) / SurfLX;
  PrimDirnCosn.YUnit.X = xnorm;
  PrimDirnCosn.YUnit.Y = ynorm;
  PrimDirnCosn.YUnit.Z = znorm;
  PrimDirnCosn.ZUnit =
      Vector3DCrossProduct(&PrimDirnCosn.XUnit, &PrimDirnCosn.YUnit);

  // primitive direction cosine assignments
  PrimDC[prim].XUnit.X = PrimDirnCosn.XUnit.X;
  PrimDC[prim].XUnit.Y = PrimDirnCosn.XUnit.Y;
  PrimDC[prim].XUnit.Z = PrimDirnCosn.XUnit.Z;
  PrimDC[prim].YUnit.X = PrimDirnCosn.YUnit.X;
  PrimDC[prim].YUnit.Y = PrimDirnCosn.YUnit.Y;
  PrimDC[prim].YUnit.Z = PrimDirnCosn.YUnit.Z;
  PrimDC[prim].ZUnit.X = PrimDirnCosn.ZUnit.X;
  PrimDC[prim].ZUnit.Y = PrimDirnCosn.ZUnit.Y;
  PrimDC[prim].ZUnit.Z = PrimDirnCosn.ZUnit.Z;

  // primitive origin: also the barcenter for a rectangular element
  PrimOriginX[prim] = SurfX;
  PrimOriginY[prim] = SurfY;
  PrimOriginZ[prim] = SurfZ;
  PrimLX[prim] = SurfLX;
  PrimLZ[prim] = SurfLZ;

  double SurfLambda = lambda;
  double SurfV = potential;

  // file output for a primitive
  if (OptPrimitiveFiles) {
    char OutPrim[256];
    strcpy(OutPrim, ModelOutDir);
    strcat(OutPrim, "/Primitives/Primitive");
    strcat(OutPrim, primstr);
    strcat(OutPrim, ".out");
    fPrim = fopen(OutPrim, "w");
    if (fPrim == NULL) {
      neBEMMessage("DiscretizeRectangle - OutPrim");
      return -1;
    }
    fprintf(fPrim, "#prim: %d, nvertex: %d\n", prim, nvertex);
    fprintf(fPrim, "Node1: %lg\t%lg\t%lg\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fPrim, "Node2: %lg\t%lg\t%lg\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fPrim, "Node3: %lg\t%lg\t%lg\n", xvert[2], yvert[2], zvert[2]);
    fprintf(fPrim, "Node4: %lg\t%lg\t%lg\n", xvert[3], yvert[3], zvert[3]);
    fprintf(fPrim, "PrimOrigin: %lg\t%lg\t%lg\n", PrimOriginX[prim],
            PrimOriginY[prim], PrimOriginZ[prim]);
    fprintf(fPrim, "Primitive lengths: %lg\t%lg\n", PrimLX[prim], PrimLZ[prim]);
    fprintf(fPrim, "Norm: %lg\t%lg\t%lg\n", xnorm, ynorm, znorm);
    fprintf(fPrim, "#volref1: %d, volref2: %d\n", volref1, volref2);
    fprintf(fPrim, "#NbSegX: %d, NbSegZ: %d\n", NbSegX, NbSegZ);
    fprintf(fPrim, "#ParentObj: %d\tEType: %d\n", SurfParentObj, SurfEType);
    fprintf(fPrim, "#SurfX\tSurfY\tSurfZ\tSurfLZ\tSurfLZ\n");
    fprintf(fPrim, "%lg\t%lg\t%lg\t%lg\t%lg\n", SurfX, SurfY, SurfZ, SurfLX,
            SurfLZ);
    // fprintf(fPrim, "#SurfRX: %lg\tSurfRY: %lg\tSurfRZ: %lg\n",
    // SurfRX, SurfRY, SurfRZ);
    fprintf(fPrim, "#DirnCosn: \n");
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.XUnit.X,
            PrimDirnCosn.XUnit.Y, PrimDirnCosn.XUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.YUnit.X,
            PrimDirnCosn.YUnit.Y, PrimDirnCosn.YUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.ZUnit.X,
            PrimDirnCosn.ZUnit.Y, PrimDirnCosn.ZUnit.Z);
    fprintf(fPrim, "#SurfLambda: %lg\tSurfV: %lg\n", SurfLambda, SurfV);
  }  // if OptPrimitiveFiles

  // necessary for gnuplot
  if (OptGnuplot && OptGnuplotPrimitives) {
    char gpPrim[256];
    strcpy(gpPrim, MeshOutDir);
    strcat(gpPrim, "/GViewDir/gpPrim");
    strcat(gpPrim, primstr);
    strcat(gpPrim, ".out");
    fgpPrim = fopen(gpPrim, "w");
    if (fgpPrim == NULL) {
      neBEMMessage("DiscretizeRectangle - OutgpPrim");
      return -1;
    }
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[2], yvert[2], zvert[2]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[3], yvert[3], zvert[3]);
    fprintf(fgpPrim, "%g\t%g\t%g\n", xvert[0], yvert[0], zvert[0]);
    fclose(fgpPrim);

    if (prim == 1)
      fprintf(fgnuPrim, " '%s\' w l", gpPrim);
    else
      fprintf(fgnuPrim, ", \\\n \'%s\' w l", gpPrim);
  }  // if OptGnuplot && OptGnuplotPrimitives

  // file outputs for elements on primitive
  if (OptElementFiles) {
    char OutElem[256];
    strcpy(OutElem, MeshOutDir);
    strcat(OutElem, "/Elements/ElemOnPrim");
    strcat(OutElem, primstr);
    strcat(OutElem, ".out");
    fElem = fopen(OutElem, "w");
    // assert(fElem != NULL);
    if (fElem == NULL) {
      neBEMMessage("DiscretizeRectangle - OutElem");
      return -1;
    }
  }  // if OptElementFiles

  // gnuplot friendly file outputs for elements on primitive
  if (OptGnuplot && OptGnuplotElements) {
    strcpy(gpElem, MeshOutDir);
    strcat(gpElem, "/GViewDir/gpElemOnPrim");
    strcat(gpElem, primstr);
    strcat(gpElem, ".out");
    fgpElem = fopen(gpElem, "w");
    // assert(fgpElem != NULL);
    if (fgpElem == NULL) {
      neBEMMessage("DiscretizeRectangle - OutgpElem");
      if (fElem) fclose(fElem);
      return -1;
    }
    // gnuplot friendly file outputs for elements on primitive
    strcpy(gpMesh, MeshOutDir);
    strcat(gpMesh, "/GViewDir/gpMeshOnPrim");
    strcat(gpMesh, primstr);
    strcat(gpMesh, ".out");
    fgpMesh = fopen(gpMesh, "w");
    if (fgpMesh == NULL) {
      neBEMMessage("DiscretizeRectangle - OutgpMesh");
      fclose(fgpElem);
      return -1;
    }
  }  // if OptGnuplot && OptElements

  // Compute element positions (CGs) in the primitive local coordinate system.
  // Then map these CGs to the global coordinate system.
  // (xav, 0, zav) is the CG of an element wrt the primitive coordinate system
  // From this, we find the offset of the element from the primitive CG in the
  // global coordinate system by just rotating the displacement vector
  // (0,0,0) -> (xav,0,zav) using the known DCM of the surface.
  // This rotated vector (now in the global system) is then added to the
  // position vector of the primitive CG (always in the global system) to get
  // the CG of the element in the global system. Make sure that the larger
  // number is being used to discretize the longer side - can be OverSmart in
  // certain cases - there may be cases where smaller number of elements suffice
  // for a longer side
  if (NbSegX == NbSegZ) {
    SurfElLX = SurfLX / NbSegX;  // element sizes
    SurfElLZ = SurfLZ / NbSegZ;
  } else if (NbSegX > NbSegZ) {
    if (SurfLX > SurfLZ) {
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    } else  // interchange NbSegX and NbSegZ
    {
      int tmp = NbSegZ;
      NbSegZ = NbSegX;
      NbSegX = tmp;
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    }
  }     // NbSegX > NbSegZ
  else  // NbSegX < NbSegZ
  {
    if (SurfLX < SurfLZ) {
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    } else  // interchange NbSegX and NbSegZ
    {
      int tmp = NbSegZ;
      NbSegZ = NbSegX;
      NbSegX = tmp;
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    }
  }

  // Analysis of element aspect ratio.
  // Note that we can afford only to reduce the total number of elements.
  // Otherwise, we'll have to realloc `EleArr' array.
  // Later, we'll make the corrections such that the total number of elements
  // remain close to the originally intended value.
  double AR = SurfElLX / SurfElLZ;  // indicative element aspect ratio
  if (OptPrimitiveFiles) {
    fprintf(fPrim,
            "Using the input, the aspect ratio of the elements on prim: %d\n",
            prim);
    fprintf(fPrim,
            "NbSegX: %d, SurfElLX: %lg, NbSegZ: %d, SurfElLZ: %lg, AR: %lg\n",
            NbSegX, SurfElLX, NbSegZ, SurfElLZ, AR);
  }
  if (AR > 10.0)  // NbSegZ needs to be reduced
  {
    double tmpElLZ = SurfElLX / 10.0;
    NbSegZ = (int)(SurfLZ / tmpElLZ);
    if (NbSegZ <= 0) NbSegZ = 1;
    SurfElLZ = SurfLZ / NbSegZ;
  }
  if (AR < 0.1)  // NbSegX need to be reduced
  {
    double tmpElLX = SurfElLZ * 0.1;
    NbSegX = (int)(SurfLX / tmpElLX);
    if (NbSegX <= 0) NbSegX = 1;
    SurfElLX = SurfLX / NbSegX;
  }
  AR = SurfElLX / SurfElLZ;  // indicative element aspect ratio
  if (OptPrimitiveFiles) {
    fprintf(fPrim, "After analyzing the aspect ratio of the elements\n");
    fprintf(fPrim,
            "NbSegX: %d, SurfElLX: %lg, NbSegZ: %d, SurfElLZ: %lg, AR: %lg\n",
            NbSegX, SurfElLX, NbSegZ, SurfElLZ, AR);
  }

  double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, xav, zav;
  ElementBgn[prim] = EleCntr + 1;
  for (int i = 1; i <= NbSegX; ++i) {
    x1 = -SurfLX / 2.0 +
         (double)(i - 1) * SurfElLX;  // assuming centroid at 0,0,0
    x2 = -SurfLX / 2.0 + (double)(i)*SurfElLX;
    xav = 0.5 * (x1 + x2);

    for (int k = 1; k <= NbSegZ; ++k) {
      z1 = -SurfLZ / 2.0 + (double)(k - 1) * SurfElLZ;
      z2 = -SurfLZ / 2.0 + (double)(k)*SurfElLZ;
      zav = 0.5 * (z1 + z2);

      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = xav;
        localDisp.Y = 0.0;
        localDisp.Z = zav;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        SurfElX = SurfX + globalDisp.X;
        SurfElY = SurfY + globalDisp.Y;
        SurfElZ = SurfZ + globalDisp.Z;
      }  // vector rotation over

      // Assign element values and write in the file
      ++EleCntr;
      if (EleCntr > NbElements) {
        neBEMMessage("DiscretizeRectangle - EleCntr more than NbElements!");
        if (fgpMesh) fclose(fgpMesh);
        return -1;
      }

      (EleArr + EleCntr - 1)->DeviceNb =
          1;  // At present, there is only one device
      (EleArr + EleCntr - 1)->ComponentNb = SurfParentObj;
      (EleArr + EleCntr - 1)->PrimitiveNb = prim;
      (EleArr + EleCntr - 1)->Id = EleCntr;
      (EleArr + EleCntr - 1)->G.Type = 4;  // rectangular here
      (EleArr + EleCntr - 1)->G.Origin.X = SurfElX;
      (EleArr + EleCntr - 1)->G.Origin.Y = SurfElY;
      (EleArr + EleCntr - 1)->G.Origin.Z = SurfElZ;
      (EleArr + EleCntr - 1)->G.LX = SurfElLX;
      (EleArr + EleCntr - 1)->G.LZ = SurfElLZ;
      (EleArr + EleCntr - 1)->G.dA =
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.LZ;
      (EleArr + EleCntr - 1)->G.DC.XUnit.X = PrimDirnCosn.XUnit.X;
      (EleArr + EleCntr - 1)->G.DC.XUnit.Y = PrimDirnCosn.XUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.XUnit.Z = PrimDirnCosn.XUnit.Z;
      (EleArr + EleCntr - 1)->G.DC.YUnit.X = PrimDirnCosn.YUnit.X;
      (EleArr + EleCntr - 1)->G.DC.YUnit.Y = PrimDirnCosn.YUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.YUnit.Z = PrimDirnCosn.YUnit.Z;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.X = PrimDirnCosn.ZUnit.X;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.Y = PrimDirnCosn.ZUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.Z = PrimDirnCosn.ZUnit.Z;
      (EleArr + EleCntr - 1)->E.Type = SurfEType;
      (EleArr + EleCntr - 1)->E.Lambda = SurfLambda;
      (EleArr + EleCntr - 1)->Solution = 0.0;
      (EleArr + EleCntr - 1)->Assigned = charge;
      (EleArr + EleCntr - 1)->BC.NbOfBCs = 1;  // assume one BC per element
      (EleArr + EleCntr - 1)->BC.CollPt.X = (EleArr + EleCntr - 1)->G.Origin.X;
      (EleArr + EleCntr - 1)->BC.CollPt.Y = (EleArr + EleCntr - 1)->G.Origin.Y;
      (EleArr + EleCntr - 1)->BC.CollPt.Z = (EleArr + EleCntr - 1)->G.Origin.Z;

      // find element vertices
      // 1) displacement vector in the ECS is first identified
      // 2) this vector is transformed to the GCS
      // 3) the global displacement vector, when added to the centroid in GCS,
      // gives the node positions in GCS.
      x0 = -0.5 *
           (EleArr + EleCntr - 1)->G.LX;  // xyz displacement of node wrt ECS
      y0 = 0.0;
      z0 = -0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x0;
        localDisp.Y = y0;
        localDisp.Z = z0;  // displacement in GCS
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x0 = (EleArr + EleCntr - 1)->G.Origin.X +
             globalDisp.X;  // xyz position in GCS
        y0 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z0 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x1 = 0.5 * (EleArr + EleCntr - 1)->G.LX;
      y1 = 0.0;
      z1 = -0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x1;
        localDisp.Y = y1;
        localDisp.Z = z1;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x1 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y1 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z1 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x2 = 0.5 * (EleArr + EleCntr - 1)->G.LX;
      y2 = 0.0;
      z2 = 0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x2;
        localDisp.Y = y2;
        localDisp.Z = z2;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x2 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y2 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z2 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x3 = -0.5 * (EleArr + EleCntr - 1)->G.LX;
      y3 = 0.0;
      z3 = 0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x3;
        localDisp.Y = y3;
        localDisp.Z = z3;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x3 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y3 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z3 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      // assign vertices of the element
      (EleArr + EleCntr - 1)->G.Vertex[0].X = x0;
      (EleArr + EleCntr - 1)->G.Vertex[0].Y = y0;
      (EleArr + EleCntr - 1)->G.Vertex[0].Z = z0;
      (EleArr + EleCntr - 1)->G.Vertex[1].X = x1;
      (EleArr + EleCntr - 1)->G.Vertex[1].Y = y1;
      (EleArr + EleCntr - 1)->G.Vertex[1].Z = z1;
      (EleArr + EleCntr - 1)->G.Vertex[2].X = x2;
      (EleArr + EleCntr - 1)->G.Vertex[2].Y = y2;
      (EleArr + EleCntr - 1)->G.Vertex[2].Z = z2;
      (EleArr + EleCntr - 1)->G.Vertex[3].X = x3;
      (EleArr + EleCntr - 1)->G.Vertex[3].Y = y3;
      (EleArr + EleCntr - 1)->G.Vertex[3].Z = z3;

      if (OptElementFiles) {
        fprintf(fElem, "##Element Counter: %d\n", EleCntr);
        fprintf(fElem, "#DevNb\tCompNb\tPrimNb\tId\n");
        fprintf(fElem, "%d\t%d\t%d\t%d\n", (EleArr + EleCntr - 1)->DeviceNb,
                (EleArr + EleCntr - 1)->ComponentNb,
                (EleArr + EleCntr - 1)->PrimitiveNb,
                (EleArr + EleCntr - 1)->Id);
        fprintf(fElem, "#GType\tX\tY\tZ\tLX\tLZ\tdA\n");
        fprintf(
            fElem, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
            (EleArr + EleCntr - 1)->G.Type, (EleArr + EleCntr - 1)->G.Origin.X,
            (EleArr + EleCntr - 1)->G.Origin.Y,
            (EleArr + EleCntr - 1)->G.Origin.Z, (EleArr + EleCntr - 1)->G.LX,
            (EleArr + EleCntr - 1)->G.LZ, (EleArr + EleCntr - 1)->G.dA);
        fprintf(fElem, "#DirnCosn: \n");
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.XUnit.X,
                (EleArr + EleCntr - 1)->G.DC.XUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.XUnit.Z);
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.YUnit.X,
                (EleArr + EleCntr - 1)->G.DC.YUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.YUnit.Z);
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.ZUnit.X,
                (EleArr + EleCntr - 1)->G.DC.ZUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.ZUnit.Z);
        fprintf(fElem, "#EType\tLambda\n");
        fprintf(fElem, "%d\t%lg\n", (EleArr + EleCntr - 1)->E.Type,
                (EleArr + EleCntr - 1)->E.Lambda);
        fprintf(fElem, "#NbBCs\tCPX\tCPY\tCPZ\tValue\n");
        fprintf(fElem, "%d\t%lg\t%lg\t%lg\t%lg\n",
                (EleArr + EleCntr - 1)->BC.NbOfBCs,
                (EleArr + EleCntr - 1)->BC.CollPt.X,
                (EleArr + EleCntr - 1)->BC.CollPt.Y,
                (EleArr + EleCntr - 1)->BC.CollPt.Z,
                (EleArr + EleCntr - 1)->BC.Value);
      }  // if OptElementFiles

      // mark centroid
      if (OptGnuplot && OptGnuplotElements) {
        fprintf(fgpElem, "%g\t%g\t%g\n", (EleArr + EleCntr - 1)->BC.CollPt.X,
                (EleArr + EleCntr - 1)->BC.CollPt.Y,
                (EleArr + EleCntr - 1)->BC.CollPt.Z);
      }  // if OptGnuplot && OptGnuplotElements

      if (OptGnuplot && OptGnuplotElements) {
        fprintf(fgpMesh, "%g\t%g\t%g\n", x0, y0, z0);
        fprintf(fgpMesh, "%g\t%g\t%g\n", x1, y1, z1);
        fprintf(fgpMesh, "%g\t%g\t%g\n", x2, y2, z2);
        fprintf(fgpMesh, "%g\t%g\t%g\n", x3, y3, z3);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x0, y0, z0);
        /*
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n", x0, y0, z0, 1);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n\n", x1, y1, z1, 2);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n", x2, y2, z2, 1);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n\n\n", x2, y2, z2, 2);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n", x0, y0, z0, 1);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n\n", x3, y3, z3, 2);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n", x2, y2, z2, 1);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n\n\n", x2, y2, z2, 2);
        */
      }  // if(OptGnuplot && OptGnuplotElements)
    }    // for k
  }      // for i
  ElementEnd[prim] = EleCntr;
  NbElmntsOnPrim[prim] = ElementEnd[prim] - ElementBgn[prim] + 1;

  if (OptPrimitiveFiles) {
    fprintf(fPrim, "Element begin: %d, Element end: %d\n", ElementBgn[prim],
            ElementEnd[prim]);
    fprintf(fPrim, "Number of elements on primitive: %d\n",
            NbElmntsOnPrim[prim]);
    fclose(fPrim);
  }

  if (OptElementFiles) fclose(fElem);

  if (OptGnuplot && OptGnuplotElements) {
    if (prim == 1)
      fprintf(fgnuElem, " '%s\' w p", gpElem);
    else
      fprintf(fgnuElem, ", \\\n \'%s\' w p", gpElem);
    if (prim == 1) {
      fprintf(fgnuMesh, " '%s\' w l", gpMesh);
      fprintf(fgnuMesh, ", \\\n \'%s\' w p ps 1", gpElem);
    } else {
      fprintf(fgnuMesh, ", \\\n \'%s\' w l", gpMesh);
      fprintf(fgnuMesh, ", \\\n \'%s\' w p ps 1", gpElem);
    }

    fclose(fgpElem);
    fclose(fgpMesh);
  }

  return (0);
}  // end of DiscretizeRectangle

int DiscretizePolygon(int /*prim*/, int /*nvertex*/, double /*xvert*/[], 
                      double /*yvert*/[], double /*zvert*/[], double /*xnorm*/,
                      double /*ynorm*/, double /*znorm*/, 
                      int /*volref1*/, int /*volref2*/, int /*inttype*/, 
                      double /*potential*/, double /*charge*/, 
                      double /*lambda*/, int NbSegX, int NbSegZ) {
  // Check inputs
  if ((NbSegX <= 0) || (NbSegZ <= 0)) {
    printf("segmentation input wrong in DiscretizePolygon ...\n");
    return -1;
  }

  return (0);
}  // end of DiscretizePolygon

int BoundaryConditions(void) {
  // Loop over all the elements of all the primitives.
  // Some of the primitives do not have a boundary condition to be satisfied
  // We'll omit these primitives but before doing that we need to know to
  // which primitive a particular element belongs to.
  // The RHS will also need to be modified to accommodate the presence of
  // floating conductors and charged (known) substances.
  // Looping on primitives (rather than elements) can save precious time!
  for (int ele = 1; ele <= NbElements; ++ele) {
    int prim = (EleArr + ele - 1)->PrimitiveNb;
    // Note this condition is pure geometry, not electric!
    switch (PrimType[prim]) {
      case 2:
      case 3:  // for floating conductor and for dielectric-dielectric intfc
      case 4:  // the BC is zero (may be modified due to known charges, later)
        (EleArr + ele - 1)->BC.Value = ApplPot[prim];
        break;

      default:
        printf("Primitive out of range in BoundaryConditions ... returning\n");
        return (-1);
    }  // switch on PrimType over
  }    // ele loop over

  return (0);
}  // end of BoundaryConditions

#ifdef __cplusplus
} // namespace
#endif
