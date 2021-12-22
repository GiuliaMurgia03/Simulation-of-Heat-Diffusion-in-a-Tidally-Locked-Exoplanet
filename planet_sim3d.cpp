
// g++ -o planet_sim3d.exe planet_sim3d.cpp spherical_grid.cpp spherical_node.cpp

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

#include"spherical_node.h"
#include"spherical_grid.h"


using namespace std;

using namespace planet_code;


int main (){

  //Crust Nature Paper
  //const double density=2700;  //kg/m^3
  //const double K=2.5; // Thermal Conductivity W/m/K

  const double density=5558;  //kg/m^3
  const double K=5; // Thermal Conductivity W/m/K

  const double cp=200; // Heat Capacity J/kg/K

  const double planet_emissivity=0.9;
  const double atmosphere_emissivity=0.9;
  
  const double Lstar=0.00154;   //solar luminosity 
  const double Dist=0.0485;       //planet distance Astronomical Units
  const double albedo=0.3;  //planet albedo

  const double Qdecay=0.01;  //micro-W/m^3 
  const double tdecay=1;     //Gyr

  
  //SPACE GRID
  spherical_grid grid;
  
  //Avoid north and south pole and the center
  grid.set_limits(0, 6880 , 0, 360, 0, 180);
  grid.set_physical_properties(Lstar, Dist, density, planet_emissivity, atmosphere_emissivity, K, cp, albedo, Qdecay, tdecay);

  //Number of Nodes
  int nr=16;
  int nphi=64;
  int ntheta=32;
  grid.set_nodes_number(nr, nphi, ntheta);

  //Initialize grid with the uniform temperature 270 Kelvin
  grid.initialize(1000);
  
  
  //TIME GRID
  
  double dt=0.1; //time step [Myr]
  double tsim; //simulation time [Gyr] 10^9 years
  
  cout<<"Enter simulation time [Gyr]: ";
  cin>>tsim;
  cout<<endl;

  //number of time steps
  int nt=tsim*1000/dt;
  cout<<"number of time steps: "<<nt<<endl;
  double t=0;

  int tsave=100;
  double tt=tsave;
  
  //time loop
  
  while(1) {

    cout<<"Computing Time "<<t<<" Myr of "<<tsim*1000<<" dt="<<dt<<"   \r";

    //Euler
    grid.set_time(t);
    grid.update_temperature_Euler_adaptive(dt);

    //midpoint
    //grid.update_temperature_midpoint_adaptive(dt);
    
    if(t==0 || t>=tt) {
      grid.save("grid_t"+to_string(int(t))+".data");
      tt=t+tsave;
    }
    
    if(t>=tsim*1000) break;
    
    t=t+dt;
  }
  
  cerr<<endl;
  
  grid.save("grid_t"+to_string(int(t))+".data");

  
  return 0;

}



