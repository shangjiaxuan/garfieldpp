#ifndef G_COMPONENT_PP_H
#define G_COMPONENT_PP_H

#include <string>

#include "Component.hh"
#include "Medium.hh"

#include <TF1.h>
#include <TF2.h>

#define LOG(x) std::cout<<x<<std::endl

namespace Garfield {

/// Component for parallel-plate geometries.

class ComponentParallelPlate : public Component {
 public:
  /// Constructor
  ComponentParallelPlate();
  /// Destructor
  ~ComponentParallelPlate() {}

  /** Define the geometry.
   * \param N amount of layers in the geometry, this includes the gas gaps \f$y\f$.
   * \param d thickness of the layers starting from the bottom to the top lauer along \f$y\f$.
   * \param eps relative permittivities of the layers starting from the bottom to the top lauer along \f$y\f$ . Here, the  gas gaps having a value of 1.
   * \param sigmaIndex Indices of the resistive layers.
   * \param V applied potential difference between the parallel plates.
   */
  void Setup(const int N, std::vector<double> eps,std::vector<double> d,const double V, std::vector<int> sigmaIndex={});

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,double& ey, double& ez, double& v, Medium*& m, int& status) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;

  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  bool GetVoltageRange(double& vmin, double& vmax) override;

  /** Add a pixel electrode.
   * \param x,z position of the center of the electrode in the xz-plane.
   * \param lx width in the along \f$x\f$.
   * \param lz width in the along \f$z\f$.
   * \param label give name using a string.
   */
  void AddPixel(double x, double z, double lx, double lz,
                const std::string& label);
  /// Add strip electrode.
  void AddStrip(double z, double lz, const std::string& label);

  /// Add plane electrode, if you want to read the signal from the cathode set
  /// the second argument to false.
  void AddPlane(const std::string& label, bool anode = true);

  // Setting the medium
  void SetMedium(Medium* medium) { m_medium = medium; }

  // This will calculate the electrode's time-dependent weighting potential on
  // the specified grid.
  void SetWeightingPotentialGrid(const std::string& label, const double xmin,
                                 const double xmax, const double xsteps,
                                 const double ymin, const double ymax,
                                 const double ysteps, const double zmin,
                                 const double zmax, const double zsteps,
                                 const double tmin, const double tmax,
                                 const double tsteps);

  // This will calculate all electrodes time-dependent weighting potential on
  // the specified grid.
  void SetWeightingPotentialGrids(const double xmin, const double xmax,
                                  const double xsteps, const double ymin,
                                  const double ymax, const double ysteps,
                                  const double zmin, const double zmax,
                                  const double zsteps, const double tmin,
                                  const double tmax, const double tsteps);

  Medium* GetMedium(const double x, const double y, const double z) override;
    
    void  GetHFunction(const double k, const double z){
        LOG("m_hIntegrant = "<<m_hIntegrant.Eval(k,z));
    }

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
                      double& xmax, double& ymax, double& zmax) override;
 private:
  static constexpr double m_precision = 1.e-30;
  static constexpr double m_Vw = 1.;
  double m_eps0 = 8.85418782e-3;
  // Voltage difference between the parallel plates.
    double m_V = 0.;
    
    bool m_debuggig = true; ///<debugging switch
    
    int m_N = 0; ///<amount of layers

    double m_upperBoundIntigration = 0;

    std::vector<double> m_eps; ///< list of irelative permitivities or layers
    std::vector<double> m_epsHolder;
    std::vector<double> m_d; ///< list of thickness of layers
    std::vector<double> m_dHolder;
    std::vector<double> m_z; ///< list of indices of conducting layers

    std::vector<int> m_sigmaIndex; ///< list of indices of conducting layers

    TF2 m_hIntegrant;

    TF1 m_wpStripIntegral; ///<Weighting potential integrant for strips
    TF2 m_wpPixelIntegral; ///<Weighting potential integrant for pixels
    
    std::vector<std::vector<std::vector<int>>> m_sigmaMatrix; // sigma_{i,j}^n, where n goes from 1 to N;
    std::vector<std::vector<std::vector<int>>> m_thetaMatrix; // theta_{i,j}^n, where n goes from 1 to N;

    std::vector<std::vector<double>> m_cMatrix; ///<c-matrixl.
    std::vector<std::vector<double>> m_vMatrix; ///<v-matrixl.
    std::vector<std::vector<double>> m_gMatrix; ///<g-matrixl.
    std::vector<std::vector<double>> m_wMatrix; ///<w-matrixl.
    
    int m_currentLayer = 0; ///<Index of the current layer.
    
  Medium* m_medium = nullptr;

  /// Structure that captures the information of the electrodes under study
  struct Electrode {
    std::string label;                     ///< Label.
    int ind = structureelectrode::NotSet;  ///< Readout group.
    double xpos, ypos;                     ///< Coordinates in x/y.
    double lx, ly;                         ///< Dimensions in the x-y plane.
    double flip = 1;                       ///< Dimensions in the x-y plane.
      
     bool m_usegrid =false;                ///< Enabeling grid based calculations.
  };

  enum fieldcomponent { xcomp = 0, ycomp, zcomp };

  /// Possible readout groups
  enum structureelectrode { NotSet = -1, Plane, Strip, Pixel };

  // Vectors storing the readout electrodes.
  std::vector<std::string> m_readout;
  std::vector<Electrode> m_readout_p;

  // Functions that calculate the electric field and potential
  double IntegratePromptField(const Electrode& el, int comp, const double x,
                        const double y, const double z);

  double IntegratePromptPotential(const Electrode& el, const double x,
                                  const double y, const double z);

  void CalculateDynamicalWeightingPotential(const Electrode& el);

  double FindWeightingPotentialInGrid(const Electrode& el, const double x,
                                      const double y, const double z);
    
    // function returning 0 if layer with specific index is conductive.
    double kroneckerDelta(const int index){
        if ( std::find(m_sigmaIndex.begin(), m_sigmaIndex.end(), index) != m_sigmaIndex.end() )
            return 0;
        else
            return 1;
    }
    
    // function construct the sigma matrix needed to calculate the w, v, c and g matrices
    bool Nsigma(int N, std::vector<std::vector<int>>& sigmaMatrix);
    
    // function construct the theta matrix needed to calculate the w, v, c and g matrices
    bool Ntheta(int N,std::vector<std::vector<int>>& thetaMatrix, std::vector<std::vector<int>>& sigmaMatrix);
    
    // function constructing the sigma an theta matrices.
    void constructGeometryMatrices(const int N);
    
    // function connstructing the w, v, c and g matrices needed for constructing the weighting potentials equations.
    void constructGeometryFunction(const int N);
    
    // obtain the index and permitivity of of the layer at hight z.
    bool getLayer(const double z,int& m,double& epsM){
        auto const it = std::upper_bound(m_z.begin(), m_z.end(), z);
        if (it == m_z.begin()||it == m_z.end()) return false;
        double mholer = *(it - 1);
        while(mholer>0&&mholer<1) mholer *=10;
        m = mholer;
        epsM= m_epsHolder[m];
        m++;
        return true;
    }
    
    // build function h needed for the integrant of the weighting potential of a stip and pixel
    void setHIntegrant();
    
    // build integrant of weighting potential of a strip
    void setwpPixelIntegrant();
    
    // build integrant of weighting potential of a pixel
    void setwpStripIntegrant();
    
    // weighting field of a plane in layer with index "indexLayer"
    double constWEFieldLayer(const int indexLayer){
        double invEz = 0;
        for(int i=1; i<=m_N-1;i++){
            invEz+=(m_z[i]-m_z[i-1])/m_epsHolder[i-1];
        }
        return 1/(m_epsHolder[indexLayer-1]*invEz);
    }
    
    // weighting potential of a plane
    double wpPlane(const double z){
        int im =-1; double epsM = -1;
        if(!getLayer(z,im,epsM)) return 0.;
        double v = 1 - (z-m_z[im-1]) * constWEFieldLayer(im);
        for(int i=1; i<=im-1;i++){
            v-=(m_z[i]-m_z[i-1])* constWEFieldLayer(i);
        }
        
        return v;
    }

    // electric field in layer with index "indexLayer"
    double constEFieldLayer(const int indexLayer){
        if(kroneckerDelta(indexLayer)==0) return 0.;
        double invEz = 0;
        for(int i=1; i<=m_N-1;i++){
            invEz+=-(m_z[i]-m_z[i-1])*kroneckerDelta(i)/m_epsHolder[i-1];
        }
        return m_V /(m_epsHolder[indexLayer-1]*invEz);
    }
    
    // function to convert decimal to binary expressed in n digits.
    bool decToBinary(int n,std::vector<int>& binaryNum);
    
    // Rebuilds c, v, g and w matrix.
    void LayerUpdate(const double z, const int im, const double epsM){
        
        if(im != m_currentLayer){
            m_currentLayer = im;
            for(int i = 0; i<im-1; i++) m_eps[i] = m_epsHolder[i];
            m_eps[im-1] =epsM; m_eps[im] =epsM;
            for(int i = im+1; i<m_N; i++) m_eps[i] = m_epsHolder[i-1];
        }
        
        double diff1 = m_z[im]-z;
        double diff2 = z-m_z[im-1];
        
        for(int i = 0; i<im-1; i++) m_d[i] = m_dHolder[i];
        m_d[im-1] =diff2; m_d[im] =diff1;
        for(int i = im+1; i<m_N; i++) m_d[i] = m_dHolder[i-1];
        
        constructGeometryFunction(m_N);
        
    };
    
  void UpdatePeriodicity() override;
  void Reset() override;
};
}  // namespace Garfield
#endif
