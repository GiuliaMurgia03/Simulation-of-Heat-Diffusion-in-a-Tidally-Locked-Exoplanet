//Node grid class

#include"spherical_node.h"

using namespace planet_code;

spherical_node::spherical_node(){}

spherical_node::~spherical_node(){}

spherical_node::spherical_node(const spherical_node& n) {
  
    T=n.T;
    Tnew=n.Tnew;

    r=n.r;
    phi=n.phi;
    theta=n.theta;

    dTdt=n.dTdt;
    d2Tdr2=n.d2Tdr2;
    d2Tdphi2=n.d2Tdphi2;
    d2Tdtheta2=n.d2Tdtheta2;
      
  }

