#pragma once

#include "ComponentFieldMap.hh"

namespace Garfield {

/// Component for importing and interpolating Comsol field maps.

class ComponentComsol : public ComponentFieldMap {
 public:
  /// Default constructor.
  ComponentComsol();
  /// Constructor from file names.
  ComponentComsol(const std::string& mesh, const std::string& mplist,
                  const std::string& field, const std::string& unit = "m");
  /// Destructor.
  ~ComponentComsol() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;

  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;
  double DelayedWeightingPotential(const double x, const double y,
                                   const double z, const double t,
                                   const std::string& label) override;
    
    void SetRange(const double xmin,const double xmax,const double ymin,const double ymax,const double zmin,const double zmax){
        m_range.set = true;
        m_range.xmin = xmin;
        m_range.ymin = ymin;
        m_range.zmin = zmin;
        
        m_range.xmax = xmax;
        m_range.ymax = ymax;
        m_range.zmax = zmax;
    }

  Medium* GetMedium(const double x, const double y, const double z) override;

  /** Import a field map.
   * \param header name of the file containing the list of nodes
   * \param mplist name of the file containing the material properties
   * \param field name of the file containing the potentials at the nodes
   * \param unit length unit
   */
  bool Initialise(const std::string& header = "mesh.mphtxt",
                  const std::string& mplist = "dielectrics.dat",
                  const std::string& field = "field.txt",
                  const std::string& unit = "m");
  /// Import the time-independent weighting field maps.
  bool SetWeightingField(const std::string& file, const std::string& label);
  /// Import the time-dependent weighting field maps.
  bool SetDelayedWeightingPotential(const std::string& file,
                                    const std::string& label);
  /// Set the time interval of the time-dependent weighting field.
  void SetTimeInterval(const double mint, const double maxt,
                       const double stept);

 protected:
  void UpdatePeriodicity() override { UpdatePeriodicityCommon(); }

  double GetElementVolume(const unsigned int i) override;
  void GetAspectRatio(const unsigned int i, double& dmin,
                      double& dmax) override;

 private:
  double m_unit = 100.;
  bool m_timeset = false;

  bool GetTimeInterval(const std::string& file);
    
    struct Range {
        bool set  = false;
        
        double xmin = 0;
        double xmax = 0;
        
        double ymin = 0;
        double ymax = 0;
        
        double zmin = 0;
        double zmax = 0;
    };
    
    Range m_range;
    
    bool CheckInRange(const double x, const double y, const double z){
        if(!m_range.set) return true;
        
        if(x<m_range.xmin||x>m_range.xmax||y<m_range.ymin||y>m_range.ymax||z<m_range.zmin||z>m_range.zmax) return false;
            
        return true;
    }
};
}  // namespace Garfield
