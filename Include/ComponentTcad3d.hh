// Interpolation in a three-dimensional field map created by Sentaurus Device

#ifndef G_COMPONENT_TCAD_3D_H
#define G_COMPONENT_TCAD_3D_H

#include <string>
#include <vector>

#include "ComponentBase.hh"

namespace Garfield {

class ComponentTcad3d : public ComponentBase {

  public:
    // Constructor
    ComponentTcad3d();
    // Destructor
    ~ComponentTcad3d() {}
    
    void ElectricField(const double x, const double y, const double z,
                       double& ex, double& ey, double& ez, double& v,
                       Medium*& m, int& status);
    void ElectricField(const double x, const double y, const double z,
                       double& ex, double& ey, double& ez,
                       Medium*& m, int& status);    
    
    bool GetMedium(const double x, const double y, const double z,
                   Medium*& medium);

    bool GetVoltageRange(double& vmin, double& vmax);
    bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
                        double& xmax, double& ymax, double& zmax); 
                        
    // Import mesh and field map from files
    bool Initialise(const std::string gridfilename, 
                    const std::string datafilename);

    int  GetNumberOfRegions() const {return nRegions;}
    void GetRegion(const int i, std::string& name, bool& active);
    void SetDriftRegion(const int i);
    void UnsetDriftRegion(const int i);
    // Set/get the medium for a given region
    void SetMedium(const int i, Medium* m);
    bool GetMedium(const int i, Medium*& m) const;

  private:
  
    // Max. number of vertices per element
    static const int nMaxVertices = 7;
  
    // Regions
    int nRegions;
    struct region {
      // Name of region (from Tcad)
      std::string name;
      // Flag indicating if the region is active (i. e. a drift medium)
      bool drift;
      Medium* medium;
    };
    std::vector<region> regions;

    // Vertices
    int nVertices;
    struct vertex {
      // Coordinates [cm]
      double x, y, z;
      // Potential [V] and electric field [V / cm]
      double p, ex, ey, ez;
      // Flag indicating if vertex belongs to more than one region
      bool   isShared;
    };
    std::vector<vertex> vertices;

    // Elements
    int nElements;
    struct element {
      // Indices of vertices
      int vertex[nMaxVertices];
      // Type of element 
      // 1: Segment (line)
      // 2: Triangle
      // 3: Rectangle
      // 4: Polygon
      // 5: Tetrahedron
      // 6: Pyramid
      // 7: Prism
      // 8: Brick
      // 9: Tetrabrick
      // 10: Polyhedron
      // Only types 2 and 5 are supported by this class.
      int type;
      // Associated region
      int region; 
    };
    std::vector<element> elements;

    // Face
    struct face {
      // Indices of edges
      int edge[4];
      int type;
    };
    
    // Voltage range
    double pMin, pMax;

    // Bounding box
    bool hasBoundingBox;
    double xMinBoundingBox, yMinBoundingBox, zMinBoundingBox;
    double xMaxBoundingBox, yMaxBoundingBox, zMaxBoundingBox;

    // Element from the previous call
    int lastElement;
    // Node point weighting factors for interpolation 
    // (local coordinates)
    double w[nMaxVertices];
    
    // Reset the component
    void Reset();
    // Periodicities
    void UpdatePeriodicity();    

    bool CheckTetrahedron(const double x, const double y, const double z, 
                          const int i);
    bool CheckTriangle(const double x, const double y, const double z, 
                       const int i);

    bool LoadGrid(const std::string gridfilename);
    bool LoadData(const std::string datafilename);
    void Cleanup();

};

}
#endif
