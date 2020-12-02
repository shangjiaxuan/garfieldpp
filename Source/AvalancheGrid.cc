#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

#include <TF1.h>

#include "Garfield/AvalancheGrid.hh"

namespace Garfield {
void AvalancheGrid::SetZGrid(Grid& av, const double ztop, const double zbottom, const int zsteps){
    
    av.zsteps = zsteps;
    av.zStepSize = (ztop-zbottom)/zsteps;
    
    for(int i = 0; i < zsteps; i++){
        
        
        av.zgrid.push_back(zbottom+i*av.zStepSize);
    }
    
    if (!m_diffusion) {
        
        av.xgrid.push_back(0.);
        std::vector<std::vector<int>> nh(zsteps,{0});
        av.n = nh;
        av.xsteps=1;
        av.transverseDiffusion={1.};
        
    }
    
    av.gridset = true;
}

void AvalancheGrid::SetXGrid(Grid& av, const double xtop, const double xbottom, const int xsteps){
    av.xsteps = xsteps;
    av.xStepSize = (xtop-xbottom)/xsteps;
    
    for(int i = 0; i < xsteps; i++){
        
        av.xgrid.push_back(xbottom+i*av.xStepSize);
        
    }
    
    std::vector<int> nh0(xsteps,0);
    std::vector<std::vector<int>> nh(av.zsteps,nh0);
    av.n = nh;
    
    DiffusionFactors(av);
}

int AvalancheGrid::GetAvalancheSize(double dx, const int nsize, const double alpha, const double eta){
    
    dx = 10*dx;
    
    int newnsize = 0;
    
    const double k = eta/alpha;
    const double ndx = exp((alpha-eta)*dx);
    
    if(nsize<1e3){
        
        for(int i = 0; i < nsize; i++){
            
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0., 1.);
            double s =dis(gen);
            
            double condition = k*(ndx-1)/(ndx-k);
            
            int Nholder;
            
            if(s<condition){
                Nholder = 0;
            } else{
                Nholder = (int) (1 + log((ndx - k)*(1 - s)/(ndx*(1 - k)))/log(1 - (1 - k)/(ndx - k)));
            }
            newnsize += Nholder;
        }
        
    }else{
        
        const double sigma =sqrt((1 + k)* ndx * (ndx - 1)/(1 - k));
        
        std::default_random_engine generator;
        std::normal_distribution<double> normalDistr(nsize*ndx,sqrt(nsize)*sigma);
        newnsize = normalDistr(generator);
    }
    
    return newnsize;
    
}

void AvalancheGrid::SnapToGrid(Grid& av, const double x, const double z, const double v){
    if(av.gridset==false) {
        
        std::cerr << m_className << "::SnapToGrid:Error: grid is not defined.\n";
        return;
        
    }
    
    av.velocity = (av.velocity*av.N +v)/(av.N +1);
    
    int index =av.zgrid.size()-1;
    
    for(int i = 0; i < av.zsteps; i++){
        
        if (z>=av.zgrid[index]){
            
            break;
            
        }else{
            
            index-=1;
            
        }
        
        if(index<0){
            
            index =0;
            break;
            
        }
    }

    av.gridPosition.push_back(index);
    
    const double nholder =GetAvalancheSize(z-av.zgrid[index],1,13,3.5);
    av.N += nholder;
    
    if(!m_diffusion){

        av.n[index][0]+=nholder;

    }else{
        
        int index2 =0;
        
        for(int i = 0; i < av.xsteps; i++){
            if (x<=av.xgrid[i]){
                
                index2=i;
                
                if(i!=0 && (x-av.xgrid[i]<-x+av.xgrid[i-1])) index2-=1;
                
                break;
            }
            
            if(i==av.xsteps-1){
                index2 =i;
            }
        }

        av.n[index][index2]+=nholder;
    }
}

void  AvalancheGrid::NextAvalancheGridPoint(Grid& av){

    int Nholder =0;
    
    av.run = false;
    
    for (int& pos : av.gridPosition){

        if(pos==0) continue;
        
        av.run=true;
        
        for (int i = 0; i < av.xsteps; i++){
            
            Nholder =av.n[pos][i];
            
            if(Nholder==0) continue;
            
            
            if(m_diffusion){
                
                double holdnsize =0.;
                        
                    if(av.N<m_MaxSize){
                        
                        holdnsize=GetAvalancheSize(av.zStepSize,av.n[pos][i],m_Townsend,m_Attachment);
                        
                    } else{
                        
                        holdnsize = av.n[pos][i];
                        m_Saturated =true;
                        
                        if(m_SaturationTime==-1) m_SaturationTime = av.time +abs(av.zStepSize/av.velocity);
                        
                    }
                
                int chargeRemaining =holdnsize;
                
                for (int j = av.transverseDiffusion.size()-1; j >= 0; j--){
                    
                    if (i-j<0||i+j>av.xsteps) continue;
                    
                    if (j>0){
                        
                        int nxd = (int) (av.transverseDiffusion[j]*holdnsize);
                        
                        av.n[pos-1][i-j]+= nxd;
                        
                        av.n[pos-1][i+j]+= nxd;
                        
                        m_sensor->AddSignal(-(nxd+Nholder)/2, av.time, av.time +av.zStepSize/av.velocity, av.xgrid[i], 0, av.zgrid[pos], av.xgrid[i-j], 0, av.zgrid[pos-1] , false,true);
                        
                        m_sensor->AddSignal(-(nxd+Nholder)/2, av.time, av.time +av.zStepSize/av.velocity, av.xgrid[i], 0, av.zgrid[pos], av.xgrid[i+j], 0, av.zgrid[pos-1] , false,true);
                        
                        chargeRemaining-=2*nxd;
                        
                    }else{
                        
                        av.n[pos-1][i]+= chargeRemaining;
                        
                        m_sensor->AddSignal(-(chargeRemaining+Nholder)/2, av.time, av.time +av.zStepSize/av.velocity, av.xgrid[i], 0, av.zgrid[pos], av.xgrid[i], 0, av.zgrid[pos-1] , false,true);
                    }
                }
                
                av.N+=holdnsize-Nholder;
                
            }else{
                if(av.N<m_MaxSize){
                    av.n[pos-1][i]=GetAvalancheSize(av.zStepSize,av.n[pos][i],m_Townsend,m_Attachment);
                } else{
                    av.n[pos-1][i] = av.n[pos][i];
                    m_Saturated =true;
                    if(m_SaturationTime==-1) m_SaturationTime = av.time +abs(av.zStepSize/av.velocity);
                }
                
                m_sensor->AddSignal(-(av.n[pos][i]+Nholder)/2, av.time, av.time +av.zStepSize/av.velocity, av.xgrid[i], 0, av.zgrid[pos], av.xgrid[i], 0, av.zgrid[pos-1] , false,true);
                
                av.N+=av.n[pos-1][i]-Nholder;
            }
            
            av.n[pos][i]=0;
            
        }
        pos-=1;
    }
    
    av.time +=abs(av.zStepSize/av.velocity);
    
}

void  AvalancheGrid::StartGridAvalanche(const double zmin,const double zmax, const int zsteps,const double xmin,const double xmax, const int xsteps) {
    
    if(!m_avmc || !m_sensor) return;
    
    if(zmin>=zmax||zsteps<=0) return;
    
    int np = m_avmc->GetNumberOfElectronEndpoints();
    
    std::cerr << m_className << "::StartGridAvalanche::Number of initial electrons = "<<np<<".\n";
    
    if (np==0) return;
    
    if (!m_diffusion){
        
        SetZGrid(m_avgrid,zmax,zmin,zsteps);
        
    }else{
        
        SetZGrid(m_avgrid,zmax,zmin,zsteps);
        SetXGrid(m_avgrid,xmax,xmin,xsteps);
        
    }
    
    double x1, y1, z1, t1;
    double x2, y2, z2, t2;
    double e1, e2;
    int status;
    
    double vel =0.;
    for (int i = 0; i < np; ++i) {
        
        m_avmc->GetElectronEndpoint(i, x1, y1, z1, t1, e1,x2, y2, z2, t2, e2, status);
        
        vel = (z2-z1)/(t2-t1);
        
        m_avgrid.time = t2;
        SnapToGrid(m_avgrid,x2,z2,vel);
        
    }
    
    sort( m_avgrid.gridPosition.begin(), m_avgrid.gridPosition.end() );
    m_avgrid.gridPosition.erase( unique( m_avgrid.gridPosition.begin(), m_avgrid.gridPosition.end() ), m_avgrid.gridPosition.end() );
    
    if (m_Velocity!=0) m_avgrid.velocity=m_Velocity;
    
    while(m_avgrid.run==true){
        
        NextAvalancheGridPoint(m_avgrid);
        
    }
    
    if (m_Saturated) std::cerr << m_className << "::StartGridAvalanche::Avalanche maximum size of "<<m_MaxSize<<" electrons reached at "<<m_SaturationTime<<" ns.\n";
    
    std::cerr << m_className << "::StartGridAvalanche::Final avalanche size = "<<m_avgrid.N<<" at t = "<<m_avgrid.time<<" ns.\n";
    
    return;
}

void AvalancheGrid::DiffusionFactors(Grid& av){
    
    if (!av.gridset || av.xStepSize<=0) return;
    
    auto  cdfunctop  = TF1("cdftop","ROOT::Math::normal_cdf(x, [0],[1])", -5, 5);
    
    cdfunctop.SetParameters(m_DiffSigma, 0.0);
    
    double factor =1;
    int index =0;
    
    while(factor>1e-3){
        
        factor = cdfunctop.Eval(0 + av.xStepSize/2 + index*av.xStepSize) - cdfunctop.Eval(0 - av.xStepSize/2 + index*av.xStepSize);
        std::cerr << m_className << "::DiffusionFactors::Transvers diffusion factor: "<<factor<<", top: "<<cdfunctop.Eval(0 - av.xStepSize/2 + index*av.xStepSize)<<", bottom: "<<cdfunctop.Eval(0 + av.xStepSize/2 + index*av.xStepSize)<<".\n";
        av.transverseDiffusion.push_back(factor);
        
        index++;
    }
    
    std::cerr << m_className << "::DiffusionFactors::Transvers diffusion spreads to "<<av.transverseDiffusion.size()<<" points.\n";
    
    return;
}

}
