/**
 * parallel_plate.cc
 * General program flow based on example code from the Garfield++ website.
 *
 * Demonstrates the importing of a parallel-plate capacitor field map
 * from Gmsh/Elmer into Garfield++.
 *
*/
#include <iostream>
#include <cmath> 
#include <cstring>
#include <fstream>
#include <TCanvas.h>
#include <TApplication.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentElmer.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewFEMesh.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

    TApplication app("app", &argc, argv);

    // Set relevant geometric parameters.
    const double ext_x = 1.0;        // external box x-width
    const double ext_y = 1.0;        // external box y-width
    const double ext_z = 1.0;        // external box z-width

    // Create a main canvas.
    TCanvas * c1 = new TCanvas();

    // Define the medium.
    MediumMagboltz* gas = new MediumMagboltz();
    gas->SetTemperature(293.15);                  // Set the temperature [K]
    gas->SetPressure(740.);                       // Set the pressure [Torr]
    gas->EnableDrift();                           // Allow for drifting in this medium
    gas->SetComposition("ar", 70., "co2", 30.);   // Specify the gas mixture (Ar/CO2 70:30)

    // Import an Elmer-created parallel plate field map.
    ComponentElmer * elm = new ComponentElmer("parallel_plate/mesh.header","parallel_plate/mesh.elements","parallel_plate/mesh.nodes",
    "parallel_plate/dielectrics.dat","parallel_plate/parallel_plate.result","cm");
    elm->SetMedium(0,gas);

    // Set up a sensor object.
    Sensor* sensor = new Sensor();
    sensor->AddComponent(elm);
    sensor->SetArea(-1*ext_x,-1*ext_y,-1*ext_z,ext_x,ext_y,ext_z);

    // Set up the object for field visualization.
    ViewField * vf = new ViewField();
    vf->SetSensor(sensor);
    vf->SetCanvas(c1);
    vf->SetArea(-1*ext_x,-1*ext_y,ext_x,ext_y);
    vf->SetNumberOfContours(20);
    vf->SetNumberOfSamples2d(30,30);
    vf->SetPlane(0,-1,0,0,0,0);

    // Set up the object for FE mesh visualization.
    ViewFEMesh * vFE = new ViewFEMesh();
    vFE->SetCanvas(c1);
    vFE->SetComponent(elm);
    vFE->SetPlane(0,0,-1,0,0,0);
    vFE->SetFillMesh(true);
    vFE->SetColor(1,kBlue);

    // Create plots.
    vFE->SetArea(-1*ext_x,-1*ext_z,-1*ext_z,ext_x,ext_z,ext_z);
    vf->PlotContour("v");
    //vFE->Plot();

    app.Run(kTRUE);

    return 0;
}
