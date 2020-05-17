#ifndef G_COMPONENT_GRID_H
#define G_COMPONENT_GRID_H

#include "ComponentBase.hh"

namespace Garfield {

/// Component for interpolating field maps on a regular mesh.

class ComponentGrid : public ComponentBase {
 public:
  /// Constructor
  ComponentGrid();
  /// Destructor
  ~ComponentGrid() {}

  void Clear() override { Reset(); }

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;
  void DelayedWeightingField(const double x, const double y, const double z,
                             const double t, double& wx, double& wy, double& wz,
                             const std::string& label) override;

  void MagneticField(const double x, const double y, const double z, double& bx,
                     double& by, double& bz, int& status) override;

  Medium* GetMedium(const double x, const double y, const double z) override;

  bool GetVoltageRange(double& vmin, double& vmax) override;
  bool GetElectricFieldRange(double& exmin, double& exmax, double& eymin,
                             double& eymax, double& ezmin, double& ezmax);
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override;

  /** Define the grid.
    * \param nx,ny,nz number of nodes along \f$x, y, z\f$.
    * \param xmin,xmax range along \f$x\f$.
    * \param ymin,ymax range along \f$y\f$.
    * \param zmin,zmax range along \f$z\f$.
    */
  bool SetMesh(const unsigned int nx, const unsigned int ny,
               const unsigned int nz, const double xmin, const double xmax,
               const double ymin, const double ymax, const double zmin,
               const double zmax);
  /// Retrieve the parameters of the grid.
  bool GetMesh(unsigned int& nx, unsigned int& ny, unsigned int& nz,
               double& xmin, double& xmax, double& ymin, double& ymax,
               double& zmin, double& zmax) const;
  /** Import electric field and potential values from a file.
    * The file is supposed to contain one line for each grid point starting with
    *   - either two or three floating point numbers,
    *     specifying the coordinates (in cm) of the grid node or
    *   - two or three integers specifying the index of the node,
    *
    * followed by
    *   - two or three floating point numbers for the electric field (in V/cm),
    * and (depending on the value of withPotential and withFlag),
    *   - a floating point number specifying the potential (in V), and
    *   - an integer flag indicating whether the point is in an active region (1) 
    *     or not (0).
    *
    * Format types are:
    *  - "xy", "xyz": nodes are specified by their coordinates
    *  - "ij", "ijk": nodes are specified by their indices
    */
  bool LoadElectricField(const std::string& filename, const std::string& format,
                         const bool withPotential, const bool withFlag,
                         const double scaleX = 1., 
                         const double scaleE = 1., const double scaleP = 1.);

  /// Import (prompt) weighting field from file.
  bool LoadWeightingField(const std::string& filename, const std::string& format,
                          const bool withPotential, 
                          const double scaleX = 1., const double scaleE = 1.,
                          const double scaleP = 1.);
  /// Import delayed weighting field from file.
  bool LoadWeightingField(const std::string& filename, const std::string& format,
                          const double time, const bool withPotential, 
                          const double scaleX = 1., const double scaleE = 1.,
                          const double scaleP = 1.);
  /// Offset coordinates in the weighting field, such that the
  /// same numerical weighting field map can be used for electrodes at
  /// different positions.
  void SetWeightingFieldOffset(const double x, const double y, const double z);


  /// Import magnetic field values from a file.
  bool LoadMagneticField(const std::string& filename, const std::string& format,
                         const double scaleX = 1., const double scaleB = 1.);

  /** Export the electric field and potential of a component to a text file.
    * \param cmp Component object for which to export the field/potential
    * \param filename name of the text file
    * \param format "xy", "xyz", "ij" or "ijk", see @ref LoadElectricField 
    */
  bool SaveElectricField(ComponentBase* cmp, const std::string& filename, 
                         const std::string& format);
  /** Export the weighting field and potential of a component to a text file.
    * \param cmp Component object for which to export the field/potential
    * \param id identifier of the weighting field
    * \param filename name of the text file
    * \param format "xy", "xyz", "ij" or "ijk", see @ref LoadElectricField 
    */
  bool SaveWeightingField(ComponentBase* cmp, const std::string& id,
                          const std::string& filename, 
                          const std::string& format);

  /// Return the field at a given node.
  bool GetElectricField(const unsigned int i, const unsigned int j,
                        const unsigned int k, double& v, double& ex, double& ey,
                        double& ez) const;

  /// Set the medium.
  void SetMedium(Medium* m);
  /// Get the medium.
  Medium* GetMedium() const { return m_medium; }
  
  
  ///new
  // Trapping

  bool ComponentGrid::LoadAttachment(const std::string& fname,
      const std::string& fmt,
      const double scaleX, int col, char particle);

  bool ComponentGrid::ElectronAttachment(const double x, const double y,
      const double z, double& att)

      bool ComponentGrid::HoleAttachment(const double x, const double y,
          const double z, double& att);
	
	
 private:
  Medium* m_medium = nullptr;
  struct Node {
    double fx, fy, fz;   				//< Field
    double v;           				//< Potential
	
  };
 ///new

  /// Electric field values and potentials.
  std::vector<std::vector<std::vector<Node> > > m_efields;
  /// Magnetic field values.
  std::vector<std::vector<std::vector<Node> > > m_bfields;
  /// Prompt weighting field values and potentials.
  std::vector<std::vector<std::vector<Node> > > m_wfields;
  /// Delayed weighting field values and potentials.
  std::vector<std::vector<std::vector<std::vector<Node> > > > m_wdfields;
  std::vector<double> m_wdtimes;
  /// new Attachment map.
  std::vector<std::vector<std::vector<double> > > m_eattachment;
  std::vector<std::vector<std::vector<double> > > m_hattachment;
  /// Active medium flag.
  std::vector<std::vector<std::vector<bool> > > m_active;

  // Dimensions of the mesh
  unsigned int m_nX = 0, m_nY = 0, m_nZ = 0;
  double m_xMin = 0., m_yMin = 0., m_zMin = 0.;
  double m_xMax = 0., m_yMax = 0., m_zMax = 0.;
  double m_dx = 0., m_dy = 0., m_dz = 0.;

  bool m_hasMesh = false;
  bool m_hasPotential = false;
  bool m_hasEfield = false;
  bool m_hasBfield = false;
  bool m_hasWfield = false;
  

  // Are all the cross-sections and concentrations valid and set.
  bool m_validTraps = false;


  // Offset for weighting field
  double m_wField_xOffset = 0.;
  double m_wField_yOffset = 0.;
  double m_wField_zOffset = 0.;

  // Voltage range
  double m_pMin = 0., m_pMax = 0.;

  /// Read/determine mesh parameters from file.
  bool LoadMesh(const std::string& filename, std::string format,
                const double scaleX);
  ///new
  /// Read data from file.
  bool LoadData(const std::string& filename, std::string format,
                const bool withPotential, const bool withFlag
                const double scaleX, const double scaleF, const double scaleP,
                std::vector<std::vector<std::vector<Node> > >& field);

  void Reset() override;
  void UpdatePeriodicity() override;

  /// Look up/interpolate the field at a given point.
  bool GetField(const double x, const double y, const double z,
                const std::vector<std::vector<std::vector<Node> > >& field,
                double& fx, double& fy, double& fz, double& p, bool& active);
  /// Reduce a coordinate to the basic cell (in case of periodicity).
  double Reduce(const double xin, const double xmin, const double xmax,
                const bool simplePeriodic, const bool mirrorPeriodic,
                bool& isMirrored) const;
  /// Set the dimensions of a table according to the mesh.
  void Initialise(std::vector<std::vector<std::vector<Node> > >& fields);
  ///new 
  bool ComponentGrid::LoadData(const std::string& filename, std::string format,
      const double scaleX,
      std::vector<std::vector<std::vector<double> > >& fields, int col);
  void ComponentGrid::Initialise(
      std::vector<std::vector<std::vector<double> > >& fields);
  bool ComponentGrid::GetAttachment(
      const double xi, const double yi, const double zi,
      const std::vector<std::vector<std::vector<double> > >& field, double& att,
      bool& active);
};
}
#endif
