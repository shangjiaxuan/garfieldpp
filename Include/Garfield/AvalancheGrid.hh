#ifndef G_AVALANCHE_GRID_H
#define G_AVALANCHE_GRID_H

#include <string>
#include <vector>
#include <iostream>

#include "GarfieldConstants.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"

namespace Garfield {

/// Calculate avalanches in an uniform electric field using avalanche statistics.

class AvalancheGrid {
public:
    /// Constructor
    AvalancheGrid() {};
    /// Destructor
    ~AvalancheGrid() {}
    /// Set the sensor.
    void SetSensor(Sensor* sensor){m_sensor=sensor;};
    /// Set the AvalancheMicroscopic.
    void SetAvalancheMicroscopic(AvalancheMicroscopic* avmc){m_avmc=avmc;};
    
    /** Start grid based avalanche simulation.
      *
      * \param zmin,zmax z coordinate range of grid [cm].
      * \param zsteps amount of z-coordinate points in grid.
      * \param xmin,xmax x coordinate range of grid [cm].
      * \param xsteps amount of x-coordinate points in grid.
      */
    void StartGridAvalanche(const double zmin,const double zmax, const int zsteps,const double xmin=0.,const double xmax=0., const int xsteps=0);
    /// Set the electron drift velocity (in cm / ns).
    void SetElectronVelocity(const double vel){m_Velocity=vel;};
    /// Set the electron Townsend coefficient (in 1 / cm).
    void SetElectronTownsend(const double town){m_Townsend=10*town;};
    /// Set the electron attachment coefficient (in 1 / cm).
    void SetElectronAttachment(const double att){m_Attachment=10*att;};
    /// Set the maximum avalanche size (1e7 by default).
    void SetMaxAvalancheSize(const double size){m_MaxSize=size;};
    /// Enable transverse diffusion of electrons with transverse diffusion coefficients (in √cm).
    void EnableDiffusion(const double diffSigma){m_diffusion=true;m_DiffSigma = diffSigma;}
    
private:
    
    double m_Townsend =13.; // [1/mm];
    
    double m_Attachment = 3.5; // [1/mm];
    
    double m_Velocity=0.; // [cm/ns]
    
    double m_MaxSize = 1e7; // Saturations size
    
    bool m_Saturated = false; // Check if avalanche has reached maximum size
    
    double m_SaturationTime = -1.; // Time when the avalanche has reached maximum size
    
    bool m_diffusion = false; // Check if transverse diffusion is enabled.
    
    double m_DiffSigma =0.; // Transverse diffusion coefficients (in √cm).
    
    std::string m_className = "AvalancheGrid";
    
    Sensor* m_sensor = nullptr;
    
    AvalancheMicroscopic* m_avmc= nullptr;
    
    struct Grid{
        
        std::vector<double> zgrid; ///<Grind points of z-coordinate.
        int zsteps=0.; ///<Amount of grid points.
        double zStepSize=0.; ///<Distance between the grind points of z-coordinate.
        
        std::vector<double> xgrid; ///<Grind points of x-coordinate.
        double xStepSize=0.; ///<Amount of grid points.
        int xsteps=0.; ///<Distance between the grind points of x-coordinate.
        
        std::vector<int> gridPosition; ///<Tracking of active z-coordinate grid points.
        
        bool gridset = false; ///<Keeps track if the grid has been defined.
        
        std::vector<std::vector<int>> n; ///<Grid based representation of space-charge.
        int N = 0; ///<Total amount of charge.
        
        std::vector<double> transverseDiffusion; ///< Factors of the charge that go to horizontally neighboring grid points.
        
        double velocity = 0; ///<Velocity of electrons.
        double time; ///<Clock.

        bool run = true; ///<Tracking if the charges are still in the drift gap.
        
    };
    
    Grid m_avgrid;
    // Setting z-coordinate grid.
    void SetZGrid(Grid& av, const double top, const double bottom, const int steps);
    // Setting x-coordinate grid.
    void SetXGrid(Grid& av, const double top, const double bottom, const int steps);
    // Get size of avalanche when going from z to z-dz.
    int GetAvalancheSize(double dz, const int nsize, const double alpha, const double eta);
    // Assign electron to the closest grid point.
    void SnapToGrid(Grid& av, const double x, const double z, const double v);
    // Go to next time step.
    void NextAvalancheGridPoint(Grid& av);
    // Compute the factor of the charge that go to neighboring points through transverse diffusion.
    void DiffusionFactors(Grid& av);
};
}

#endif
