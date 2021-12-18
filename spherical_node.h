//Node grid class


#ifndef _SPHERICAL_NODE_
#define _SPHERICAL_NODE_


namespace planet_code {
  
class spherical_node {
  
public:

  spherical_node(); //Constructor
  ~spherical_node();//Destructor
  
  spherical_node(const spherical_node& n); //copy constructor
  
  double T; //Temperature
  double Tnew; //New Temperature value
  double Told; //This help with midpoint
  double r; //radius or distance from center 
  double phi; //longitude
  double theta; //latitude

  double dTdt;  //derivative w.r.t. time
  
  double d2Tdr2; //second derivative w.r.t. radius
  double d2Tdphi2; //second derivative w.r.t. longitude
  double d2Tdtheta2; //second derivative w.r.t latitude
  
};

}

#endif
