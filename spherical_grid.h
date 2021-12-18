//Grid Class

#ifndef _SPHERICAL_GRID_
#define _SPHERICAL_GRID_
#include"spherical_node.h"
#include<vector>
#include<cmath>
#include<fstream>

using namespace std;

namespace planet_code {

class spherical_grid {

 private:

  int nr; // number of radial nodes 
  int nphi; //number of longitude nodes
  int ntheta; //number of latitude nodes

  double dr; //radial resolution
  double dphi; //longitude resolution
  double dtheta; //latitude resolution

  double radius_min; //km
  double radius_max; //km

  double theta_min;
  double theta_max;

  double phi_min;
  double phi_max;

  double radconst;

  double radconst2;

  double albedo, density, emissivity, Ktherm, cp, alpha;

  double Lstar, Dist;

  double Qdecay0, Qdecay, tdecay;
  
  double Teq;

  double age;
  
 public:
  
  spherical_grid();
  ~spherical_grid();

   spherical_grid(const spherical_grid& g); //copy constructor

  void set_limits(double r_min, double r_max, double p_min, double p_max, double t_min, double t_max);
  void set_physical_properties(double Lstar, double Dist, double density, double emissivity, double Ktherm, double cp, double albedo, double Qdecay, double tdecay);

  vector<spherical_node> vnodes;

  void initialize(double T0); //Initialize nodes coordinates and set a uniform temoerature T0

  void set_nodes_number (int n_r, int n_phi, int n_theta);

  void calculate_nodes_derivatives();

  double calculate_time_derivative(int k, int i, int j, int index);
  void calculate_radial_derivatives(int k, int i, int j);
  void calculate_longitudinal_derivatives(int k, int i, int j);
  void calculate_latitudinal_derivatives(int k, int i, int j);

  int node_index (int k, int i, int j);

  void update_temperature_Euler_adaptive(double& dt);
  void update_temperature_midpoint_adaptive(double& dt);
  void update_temperature_Euler(double dt);
  void update_temperature_midpoint(double dt);

  void update_radiative_cooling(double dt);

  void set_time(double);
  double get_time();
  
  double get_Tmax();
  double get_Tmin();
  
  double get_dr();
  double get_dphi();
  double get_dtheta();

  double get_radius();
  
  void save(string outfile);
  void load(string infile);
};




}

#endif



 