#include "Garfield/ComponentParallelPlate.hh"

#include <TF1.h>
#include <TF2.h>

#include <cmath>
#include <limits>
#include <iostream>

#include "Garfield/GarfieldConstants.hh"

#define LOG(x) std::cout<<x<<std::endl

namespace Garfield {

ComponentParallelPlate::ComponentParallelPlate() : Component("ParallelPlate") {}

void ComponentParallelPlate::Setup(const int N, std::vector<double> eps,std::vector<double> d,const double V, std::vector<int> sigmaIndex) {
    
  // Here I switch conventions with the z-axis the direction of drift.
    
    const int Nholder1 = eps.size();
    const int Nholder2 = d.size();
    if(N!=Nholder1||N!=Nholder2){
        LOG("ComponentParallelPlate::Setup:: Inconsistency between the number of layers, permittivities and thicknesses given.");
        return;
    }else if(N<2){
        LOG("ComponentParallelPlate::Setup:: Number of layers must be larger then 1.");
        return;
    }
    
    if(m_debuggig) LOG("ComponentParallelPlate::Setup:: Loading parameters.");
    m_eps = eps;
    //std::for_each(m_eps.begin(), m_eps.end(), [](double &n){ n*=m_eps0; });
    
    m_d = d;
    m_N = N;
    m_V = V;
    
    m_sigmaIndex = sigmaIndex;
    
    m_upperBoundIntigration = *max_element(std::begin(m_d), std::end(m_d));
    
    std::vector<double>  m_zHolder(N+1);
    m_zHolder[0] = 0;
    for(int i = 1; i<=N; i++){
        m_zHolder[i]=m_zHolder[i-1]+m_d[i-1];
        
        if(m_debuggig) LOG("ComponentParallelPlate::Setup:: layer "<<i<<":: z = "<<m_zHolder[i]);
    }
    m_z = m_zHolder;
    
    if(m_debuggig) LOG("ComponentParallelPlate::Setup:: Constructing matrices");
    constructGeometryFunction(m_N);
    
    if(m_debuggig) LOG("ComponentParallelPlate::Setup:: Computing weighting potential functions.");
    setHIntegrant();
    setwpStripIntegrant();
    setwpPixelIntegrant();
    
    LOG("ComponentParallelPlate::Setup:: Geometry with N = "<< N <<" layers set.");
}

bool ComponentParallelPlate::GetBoundingBox(double& x0, double& y0, double& z0,
                                            double& x1, double& y1,
                                            double& z1) {
  // If a geometry is present, try to get the bounding box from there.
    
    // Here I switch conventions back with the y-axis the direction of drift.
  if (m_geometry) {
    if (m_geometry->GetBoundingBox(x0, y0, z0, x1, y1, z1)) return true;
  }
  z0 = -std::numeric_limits<double>::infinity();
  x0 = -std::numeric_limits<double>::infinity();
  z1 = +std::numeric_limits<double>::infinity();
  x1 = +std::numeric_limits<double>::infinity();
  // TODO: check!
  y0 = 0.;
  y1 = m_d[m_N-1];
  return true;
}

double ComponentParallelPlate::IntegratePromptField(const Electrode& el, int comp,
                                              const double x, const double y,
                                              const double z) {
    return 0;
}

double ComponentParallelPlate::IntegratePromptPotential(const Electrode& el,
                                                        const double x,
                                                        const double y,
                                                        const double z) {
  switch (el.ind) {
    case structureelectrode::Plane: {
        return wpPlane(z);
      break;
    }
    case structureelectrode::Pixel: {
        m_wpPixelIntegral.SetParameters(x,y,el.xpos,el.ypos,el.lx,el.ly,z); //(x,y,x0,y0,lx,ly,z)
        int im; double epsm;
        getLayer(z,im,epsm);
        double upLim = 10*std::abs(z-m_z[im]);
        return m_wpPixelIntegral.Integral(0,upLim,0,upLim,1.e-12);
      break;
    }
    case structureelectrode::Strip: {
        m_wpStripIntegral.SetParameters(x,el.xpos,el.lx,z); //(x,x0,lx,z)
        int im; double epsm;
        getLayer(z,im,epsm);
        double upLim = 10*std::abs(z-m_z[im]);
        return m_wpStripIntegral.Integral(0,upLim,1.e-12);
      break;
    }
    default: {
      std::cerr << m_className << "::IntegratePromptPotential:\n"
                << "    Unknown electrode type.\n";
      return 0.;
    }
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double& ex,
                                           double& ey, double& ez, Medium*& m,
                                           int& status) {
    // Here I switch conventions back with the y-axis the direction of drift.
    
    ex = ey = ez = 0;
    
    int im =-1; double epsM = -1;
    if(!getLayer(y,im,epsM)){
        status = -6;
        return;
    }
    
    ey = constEFieldLayer(im);

  m = m_geometry ? m_geometry->GetMedium(x, y, z) : m_medium;

  if (!m) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField: No medium at (" << x << ", "
                << y << ", " << z << ").\n";
    }
    status = -6;
    return;
  }

  if (epsM==1) {
    status = 0;
  } else {
    status = -5;
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double& ex,
                                           double& ey, double& ez, double& v,
                                           Medium*& m, int& status) {
    // Here I switch conventions back with the y-axis the direction of drift.
    
    ex = ey = ez = v = 0;
    
    int im =-1; double epsM = -1;
    if(!getLayer(y,im,epsM)){
        status = -6;
        return;
    }
    ey = constEFieldLayer(im);
    
    // TODO: check sign
    v = -m_V - (z-m_z[im-1]) * constEFieldLayer(im);
    for(int i=1; i<=im-1;i++){
        v-=(m_z[i]-m_z[i-1])* constEFieldLayer(i);
    }

  m = m_geometry ? m_geometry->GetMedium(x, y, z) : m_medium;

  if (!m) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField: No medium at (" << x << ", "
                << y << ", " << z << ").\n";
    }
    status = -6;
    return;
  }

  if (epsM==1) {
    status = 0;
  } else {
    status = -5;
  }
}

bool ComponentParallelPlate::GetVoltageRange(double& vmin, double& vmax) {
  if (m_V == 0) return false;

  if (m_V < 0) {
    vmin = m_V;
    vmax = 0;
  } else {
    vmin = 0;
    vmax = m_V;
  }
  return true;
}

void ComponentParallelPlate::WeightingField(const double x, const double y,
                                            const double z, double& wx,
                                            double& wy, double& wz,
                                            const std::string& label) {
    
    // Here I switch conventions back with the y-axis the direction of drift.
    
  wx = 0;
  wy = 0;
  wz = 0;
    // TODO: implement weighting field
  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      wx = electrode.flip *
        IntegratePromptField(electrode, fieldcomponent::ycomp, z, x, y);
      wy = electrode.flip *
        IntegratePromptField(electrode, fieldcomponent::zcomp, z, x, y);
      wz = electrode.flip *
        IntegratePromptField(electrode, fieldcomponent::xcomp, z, x, y);
    }
  }
}

double ComponentParallelPlate::WeightingPotential(const double x,
                                                  const double y,
                                                  const double z,
                                                  const std::string& label) {
    
    // Here I switch conventions back with the y-axis the direction of drift.
    
  double ret = 0.;

  for (const auto& electrode : m_readout_p) {
    if (electrode.label == label) {
      if (!electrode.m_usegrid) {
        ret += electrode.flip * IntegratePromptPotential(electrode, z, x, y);
      } else {
        ret += FindWeightingPotentialInGrid(electrode, z, x, y);
      }
    }
  }
  return ret;
}

void ComponentParallelPlate::Reset() {
  m_readout.clear();
  m_readout_p.clear();
    
    m_cMatrix.clear();
    m_vMatrix.clear();
    m_gMatrix.clear();
    m_wMatrix.clear();
    
    m_sigmaIndex.clear();
    m_eps.clear();
    m_d.clear();
    m_z.clear();
    
    m_N = 0;
    m_V = 0;
    
    m_medium = nullptr;
}

void ComponentParallelPlate::UpdatePeriodicity() {
  if (m_debug) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Periodicities are not supported.\n";
  }
}

void ComponentParallelPlate::AddPixel(double x, double z, double lx_input,
                                      double lz_input,
                                      const std::string& label) {
    
    // Here I switch conventions back with the y-axis the direction of drift.
    
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it != m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddPixel:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode pixel;
  pixel.label = label;
  pixel.ind = structureelectrode::Pixel;
  pixel.xpos = z;
  pixel.ypos = x;
  pixel.lx = lz_input;
  pixel.ly = lx_input;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(pixel));
  std::cout << m_className << "::AddPixel: Added pixel electrode.\n";
}

void ComponentParallelPlate::AddStrip(double z, double lz_input,
                                      const std::string& label) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it != m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddStrip:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode strip;
  strip.label = label;
  strip.ind = structureelectrode::Strip;
  strip.xpos = z;
  strip.lx = lz_input;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(strip));

  std::cout << m_className << "::AddStrip: Added strip electrode.\n";
}

void ComponentParallelPlate::AddPlane(const std::string& label, bool anode) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it != m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddPlane:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode plate;
  plate.label = label;
  plate.ind = structureelectrode::Plane;

  if (!anode) plate.flip = -1;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(plate));

  std::cout << m_className << "::AddPlane: Added plane electrode.\n";
}

Medium* ComponentParallelPlate::GetMedium(const double x, const double y,
                                          const double z) {
  if (m_geometry) {
    return m_geometry->GetMedium(x, y, z);
  } else if (m_medium) {
    return m_medium;
  }
  return nullptr;
}

bool ComponentParallelPlate::Nsigma(int N, std::vector<std::vector<int>>& sigmaMatrix){
    int nColomb = N-1;
    int nRow = pow(2,N-1);
    // array to store binary number
    std::vector<int> binaryNum(nColomb,0);
    
    for(int i=0; i<nRow;i++){
        if(decToBinary(i,binaryNum)){
            sigmaMatrix.push_back(binaryNum);
            std::reverse(sigmaMatrix[i].begin(),sigmaMatrix[i].end());
            std::for_each(sigmaMatrix[i].begin(), sigmaMatrix[i].end(), [](int &n){ n=1-2*n; });
            
        }
        
    }
    return true;
}

bool ComponentParallelPlate::Ntheta(int N,std::vector<std::vector<int>>& thetaMatrix, std::vector<std::vector<int>>& sigmaMatrix){
    int nColomb = N-1;
    int nRow = pow(2,N-1);
    
    std::vector<int> thetaRow(nColomb,1);
    std::vector<int> thetaRowReset(nColomb,1);
    
    for(int i=0; i<nRow;i++){
        for(int j =0; j<nColomb; j++){
            for(int l =j; l<nColomb;l++) thetaRow[j]*=sigmaMatrix[i][l];
        }
        thetaMatrix.push_back(thetaRow);
        thetaRow=thetaRowReset;
    }
    return true;
}

void ComponentParallelPlate::constructGeometryFunction(const int N){
    
    int nRow = N;int nColomb = pow(2,N-1);
    
    std::vector<std::vector<int>> sigmaMatrix; // sigma_{i,j}^n, where n is defined bellow
    std::vector<std::vector<int>> thetaMatrix; // theta_{i,j}^n
    
    std::vector<std::vector<int>> sigmaMatrix2; // sigma_{i,j}^(N+1-n)
    std::vector<std::vector<int>> thetaMatrix2; // theta_{i,j}^(N+1-n)
    
    std::vector<double> cHold(nColomb,1);
    std::vector<double> vHold(nColomb,0);
    std::vector<double> gHold(nColomb,1);
    std::vector<double> wHold(nColomb,0);
    
    for(int n=1; n<=nRow;n++){
        // building sigma and theta matrices
        Nsigma(n,sigmaMatrix);
        Ntheta(n,thetaMatrix,sigmaMatrix);
        
        Nsigma(N-n+1,sigmaMatrix2);
        Ntheta(N-n+1,thetaMatrix2,sigmaMatrix2);
        
        int ix1 =0; int ix2 =0;
        
        for(int i=0; i<nColomb;i++){
            // cyclic permutation over the rows of sigma
            if(ix1==pow(2,n-1)) ix1 = 0;
            if(ix2==pow(2,N-n)) ix2 = 0;
            // normalization
            cHold[i]*=1/pow(2,n-1);
            gHold[i]*=1/pow(2,N-n);
            // summation for c and v
            if(n>1){
                for(int j = 0;j<n-1;j++){
                    cHold[i]*=(m_eps[j]+sigmaMatrix[ix1][j]*m_eps[j+1])/m_eps[j+1];
                    vHold[i]+=(thetaMatrix[ix1][j]-1)*m_d[j];
                }
            }
            // summation for g and w
            for(int j = 0;j<N-n;j++){
                gHold[i]*=(m_eps[j]+sigmaMatrix2[ix2][j]*m_eps[j+1]/m_eps[j+1]);
                wHold[i]+=(thetaMatrix2[ix2][j]-1)*m_d[j];
            }
            ix1++;
            ix2++;
        }
        
        // store solution for row n
        m_cMatrix.push_back(cHold);
        m_vMatrix.push_back(vHold);
        m_gMatrix.push_back(gHold);
        m_wMatrix.push_back(wHold);
        
        // reset
        std::for_each(cHold.begin(), cHold.end(), [](double &n){ n=1; });
        std::for_each(vHold.begin(), vHold.end(), [](double &n){ n=0; });
        std::for_each(gHold.begin(), gHold.end(), [](double &n){ n=1; });
        std::for_each(wHold.begin(), wHold.end(), [](double &n){ n=0; });
        
        sigmaMatrix.clear();
        sigmaMatrix2.clear();
        thetaMatrix.clear();
        thetaMatrix2.clear();
        
    }
    
}

void ComponentParallelPlate::setHIntegrant(){
    auto hFunction = [=](double* k, double* /*p*/) {
        double kk = k[0];
        double z = k[1];
        
        double hNorm=0;
        double h =0;
        
        int im =-1; double epsM = -1;
        if(!getLayer(z,im,epsM)) return 0.;
      //  im-=2;
        if(pow(2,m_N-im-1)<1) return 0.;
        for(int i=0; i<pow(2,m_N-im-1);i++){
            h+=m_gMatrix[im][i]*sinh(kk*(m_wMatrix[im][i]+m_z[m_N]-z));
        }
        for(int i=0; i<pow(2,m_N-1);i++){
            hNorm+=m_cMatrix[m_N-1][i]*sinh(kk*(m_vMatrix[m_N-1][i]+m_z[m_N]));
        }
        return h*m_eps[0]/(m_eps[m_N-1]*hNorm);
    };
    TF2* hF = new TF2("hFunction", hFunction, 0, 10*m_upperBoundIntigration,0, m_z[m_N], 0);
    
    hF->Copy(m_hIntegrant);
    
    delete hF;
    
}

void ComponentParallelPlate::setwpPixelIntegrant(){
    auto intFunction = [=](double* k, double* p) {
        double kx = k[0];
        double ky = k[1];
        
        double K = sqrt(kx*kx+ky*ky);
        
        double x = p[0];
        double y = p[1];
        double x0 = p[2];
        double y0 = p[3];
        double wx = p[4];
        double wy = p[5];
        double z = p[6];
        
        double sol = cos(kx*(x-x0))*sin(kx*wx/2)*cos(ky*(y-y0))*sin(ky*wy/2)*m_hIntegrant.Eval(K,z)/(kx*ky);
        
        return 4*sol/(Pi*Pi);
    };
    
    TF2* wpPixelIntegrant = new TF2("wpPixelIntegrant", intFunction, 0, 20*m_upperBoundIntigration,0, 20*m_upperBoundIntigration, 7);
    wpPixelIntegrant->SetNpx(10000); // increasing number of points the function is evaluated on
    wpPixelIntegrant->SetNpy(10000);
    wpPixelIntegrant->Copy(m_wpPixelIntegral);
    
    delete wpPixelIntegrant;
}

void ComponentParallelPlate::setwpStripIntegrant(){
    auto intFunction = [=](double* k, double* p) {
        double kk = k[0];
        double x = p[0];
        double x0 = p[1];
        double wx = p[2];
        double z = p[3];
        double sol = cos(kk*(x-x0))*sin(kk*wx/2)*m_hIntegrant.Eval(kk,z)/kk;
        return 2*sol/Pi;
    };
    TF1* wpStripIntegrant = new TF1("wpStripIntegrant", intFunction, 0, 20*m_upperBoundIntigration, 4);
    wpStripIntegrant->SetNpx(1000); // increasing number of points the function is evaluated on
    wpStripIntegrant->Copy(m_wpStripIntegral);
    
    delete wpStripIntegrant;
}

bool ComponentParallelPlate::decToBinary(int n,std::vector<int>& binaryNum)
{
    int L = binaryNum.size();
    // counter for binary array
    int i = 0;
    while (n > 0) {
        if(i+1>L){
            LOG("Size of binary exceeds amount of colomb.");
            return false; // Triggered if binary expression is larger then n.
        }
        // storing remainder in binary array
        binaryNum[i] =n % 2;
        n = n / 2;
        i++;
    }
    return true; // Succesfully
}

void ComponentParallelPlate::SetWeightingPotentialGrid(
    const std::string& label, const double xmin, const double xmax,
    const double xsteps, const double ymin, const double ymax,
    const double ysteps, const double zmin, const double zmax,
    const double zsteps, const double tmin, const double tmax,
    const double tsteps) {
    
    // TODO: Use existing classes for a grid based field map!
 
}

void ComponentParallelPlate::SetWeightingPotentialGrids(const double xmin, const double xmax,
const double xsteps, const double ymin,
const double ymax, const double ysteps,
const double zmin, const double zmax,
const double zsteps, const double tmin,
const double tmax, const double tsteps){
    
    // TODO: Use existing classes for a grid based field map!
   
}

double ComponentParallelPlate::FindWeightingPotentialInGrid(const Electrode& el,
                                                            const double x,
                                                            const double y,
                                                            const double z) {
    
    // TODO: Use existing classes for a grid based field map!
    
    return 0.;
}

}  // namespace Garfield
