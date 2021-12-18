//Grid Class

#include"spherical_grid.h"
#include<iostream>

using namespace std;

using namespace planet_code;


spherical_grid::spherical_grid() {}

spherical_grid::~spherical_grid() {}



spherical_grid::spherical_grid(const spherical_grid& g) {

  nr=g.nr;
  nphi=g.nphi;
  ntheta=g.ntheta;

  dr=g.dr;
  dphi=g.dphi;
  dtheta=g.dtheta;

  radius_min=g.radius_min;
  radius_max=g.radius_max;

  theta_min=g.theta_min;
  theta_max=g.theta_max;

  phi_min=g.phi_min;
  phi_max=g.phi_max;

  vnodes=g.vnodes;  
}

void spherical_grid::set_limits(double r_min, double r_max,
				double p_min, double p_max,
				double t_min, double t_max) {
  radius_min=r_min;
  radius_max=r_max;
  phi_min=p_min*M_PI/180.0;
  phi_max=p_max*M_PI/180.0;
  theta_min=t_min*M_PI/180.0;
  theta_max=t_max*M_PI/180.0;
  radconst=0;
}

void spherical_grid::set_time(double t) {

  age=t;
  
  Qdecay=Qdecay0*exp(-age/tdecay);


  //cerr<<t<<" "<<tdecay<<" "<<Qdecay<<" "<<Qdecay0<<endl;
}

double spherical_grid::get_time() {
  return age;
}


void spherical_grid::set_physical_properties(double Ls, double D, double d,
					     double e, double K, double c,
					     double a, double qd, double qt) {

  Lstar=Ls*3.828e26;   //Watt
  Dist=D*1.495978707e11; //m

  double solar_constant=(1-albedo)*Lstar/(4*M_PI*pow(Dist,2));
  cout<<"solar_constant="<<solar_constant<<" W/m^2"<<endl;
  
  albedo=a;
  density=d; //kg/m^3
  emissivity=e;
  Ktherm=K; // Thermal Conductivity W/m/K
  cp=c; // Heat Capacity J/kg/K
 
  alpha=K/(cp*density);   //Diffusivity coefficient m^2/s
  cout<<"alpha= "<<alpha<<" m^2/s"<<" or "<<alpha*1e6<<" mm^2/s"<<endl;
  
  double sigma=5.670367e-8; // W/m^2/K^4
  radconst=emissivity*sigma/(cp*density); // m/K^3/s

  radconst=radconst*1e-3*3.1536e13; // km/K^3/Myr

  Teq=pow((1-albedo)*Lstar/(emissivity*sigma*4*M_PI*pow(Dist,2)),1.0/4.0);

  cout<<"Teq="<<Teq<<" K"<<endl;

  tdecay=qt*1000;  //Myr;
  Qdecay0=qd*1e-6*3.1536e13/(cp*density);   // K/Myr
}


//Function that calculate the second derivatives of the temperature with the respect to radius, the angle phi and the angle theta

void spherical_grid::calculate_nodes_derivatives() {

  for(int k=0; k<nr; k++) {
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {
  
	calculate_radial_derivatives(k,i,j);
	calculate_longitudinal_derivatives(k,i,j);
	calculate_latitudinal_derivatives(k,i,j);

      }
    }
  }

}

  
void spherical_grid::calculate_radial_derivatives(int k, int i, int j) {

  int index=node_index(k,i,j);
  double Tcenter= vnodes[index].T;
  
  int index_after, index_before;
  double Tafter, Tbefore;
  
  //boundary conditions
  
  if(k==0) {
    index_before=node_index(k+1,i,j);//Internal Boundary condition: no heat flux, i.e. the first derivative is equal to 0
    index_after=node_index(k+1,i,j);
    Tbefore=vnodes[index_before].T;
    Tafter=vnodes[index_after].T;
  }
  
  else if (k==(nr-1)){

    /*
    //illuminated face is constantly at Teq
    if(vnodes[index].phi<0.5*M_PI || vnodes[index].phi>1.5*M_PI ) {
      
      Tcenter=Teq*pow(fabs(sin(vnodes[index].theta)*cos(vnodes[index].phi)),1.0/4.0);
      vnodes[index].d2Tdr2=0;
      return;
    }
    */
    
    index_before=node_index(k-1,i,j);
    
    Tbefore=vnodes[index_before].T;

    Tafter=Tbefore;
  }

  else {
    
    index_before=node_index(k-1,i,j);
    index_after=node_index(k+1,i,j);
    
    Tbefore=vnodes[index_before].T;
    Tafter=vnodes[index_after].T;
  }
  
  
  //calculate the second derivative of the temperature with the respect to radius r (d2Tdr2)
  double r=vnodes[index].r;

  vnodes[index].d2Tdr2=(Tafter-Tbefore)/(r*dr)+(Tbefore+Tafter-2*Tcenter)/(dr*dr);  
}
 



void spherical_grid::calculate_longitudinal_derivatives(int k, int i, int j) {

	int index=node_index(k,i,j);
	int index_before;
	int index_after;

	//periodic boundary conditions
      
	if(i==0) {

	  index_before=node_index(k,nphi-1,j);
	  index_after=node_index(k,i+1,j);
	}

	else if (i==(nphi-1)){

	  index_before=node_index(k,i-1,j);
	  index_after=node_index(k,0,j);
	}

	else {
	  
	  index_before=node_index(k,i-1,j);
	  index_after=node_index(k,i+1,j);
	}
	
	//calculate the second derivative of the temperature with the respect to the angle phi (d2Tdphi2)
	
	double Tbefore= vnodes[index_before].T;
	double Tcenter= vnodes[index].T;
	double Tafter= vnodes[index_after].T;

	double r=vnodes[index].r;
	double theta=vnodes[index].theta;
	
	vnodes[index].d2Tdphi2= (1.0/pow(r*sin(theta),2))*(Tbefore+Tafter-2*Tcenter)/(dphi*dphi);
}


void spherical_grid::calculate_latitudinal_derivatives(int k, int i, int j) {


	int index=node_index(k,i,j);
	int index_before;
	int index_after;

	//boundary conditions
      
	if(j==0) {
	  index_before=node_index(k,i,j+1); //External Boundary condition: no heat flux, i.e. the first derivative is equal to 0
	  index_after=node_index(k,i,j+1);
	}

	else if (j==(ntheta-1)){
	  index_before=node_index(k,i,j-1); //External Boundary condition: no heat flux, i.e. the first derivative is equal to 0
	  index_after=node_index(k,i,j-1);
	}

	else {
	  index_before=node_index(k,i,j-1);
	  index_after=node_index(k,i,j+1);
	}

	
	//calculate the second derivative of the temperature with the respect to theta r (d2Tdtheta2)

	double Tbefore= vnodes[index_before].T;
	double Tcenter= vnodes[index].T;
	double Tafter= vnodes[index_after].T;
	double r=vnodes[index].r;
	double theta=vnodes[index].theta;

	vnodes[index].d2Tdtheta2=(1.0/(r*r*sin(theta)))*(cos(theta)*(Tafter-Tbefore)/(2*dtheta)+sin(theta)*(Tbefore+Tafter-2*Tcenter)/(dtheta*dtheta));
}

double spherical_grid::calculate_time_derivative(int k, int i, int j, int index) {
  
  
  // calculate the second spatial derivatives of the temperature with the respect to radius, the angle phi and the angle theta
  
  calculate_radial_derivatives(k,i,j);
  calculate_longitudinal_derivatives(k,i,j);
  calculate_latitudinal_derivatives(k,i,j);	  
  
  double dTdt=(alpha*3.1536e7)*(vnodes[index].d2Tdr2+
				vnodes[index].d2Tdphi2+
				vnodes[index].d2Tdtheta2)+Qdecay;

  return dTdt;
}



void spherical_grid::initialize(double T0) {
  
  for(int k=0; k<nr; k++) { 
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {

	spherical_node n;

	n.r=radius_min+k*dr+0.5*dr;
	n.phi=phi_min+i*dphi+0.5*dphi;
	n.theta=theta_min+j*dtheta+0.5*dtheta;
	
	//n.T=T0+100*(radius_max-n.r)/(radius_max-radius_min);
		
	//n.T=T0+100*sin(n.theta);

	//n.T=T0+100*sin(n.phi);
	
	//n.T=T0+150*sin(n.theta)*fabs(cos(0.5*n.phi))+50*(radius_max-n.r)/(radius_max-radius_min);

	n.T=T0;
   
	if (k==(nr-1)){
	  
	  if(n.phi<0.5*M_PI || n.phi>1.5*M_PI ) {
	    
	    n.T=Teq*pow(fabs(sin(n.theta)*cos(n.phi)),1.0/4.0);
	    
	  }

	  if(n.T<T0) n.T=T0;
	}
	
	
	n.Told=n.T;
	
	int index= node_index(k,i,j); //Node position in the vector
	
	vnodes[index]=n;
      }
    }
  }
  

}



int spherical_grid::node_index (int k, int i, int j) {
  
  return k+i*nr+j*nr*nphi;
}


void spherical_grid::set_nodes_number (int n_r, int n_phi, int n_theta) {

  nr=n_r;
  nphi=n_phi;
  ntheta=n_theta;

  vnodes.resize(nr*nphi*ntheta);
  
  cout<<"number of grid nodes: "<<vnodes.size()<<endl;

  //calculate radial resolution
  dr=(radius_max-radius_min)/nr;
  
  cout<<"radius min="<<radius_min<<" max="<<radius_max<<" dr: "<<dr<<endl;

  //calculate angular longitude resolution
  dphi=2*M_PI/nphi; //rad
  cout<<"phi min="<<phi_min*180/M_PI<<" max="<<phi_max*180/M_PI<<" dphi: "<<dphi*180/M_PI<<"°"<<endl;
  
  //calculate angular latitude resolution
  dtheta=(theta_max-theta_min)/(ntheta); //rad
  
  
  cout<<"theta min="<<theta_min*180/M_PI<<" max="<<theta_max*180/M_PI<<" dtheta: "<<dtheta*180/M_PI<<"°"<<endl;
  
}

double spherical_grid::get_radius() {
  return radius_max;
}

double spherical_grid::get_Tmax() {

  double Tmax=vnodes[0].T;

  for(int i=0; i<vnodes.size(); i++) {
    if(vnodes[i].T>Tmax) Tmax=vnodes[i].T;
  }

  return Tmax;
}

double spherical_grid::get_Tmin() {

  double Tmin=vnodes[0].T;

  for(int i=0; i<vnodes.size(); i++) {
    if(vnodes[i].T<Tmin) Tmin=vnodes[i].T;
  }

  return Tmin;
}

double spherical_grid::get_dr() {
  return dr;
}

double spherical_grid::get_dphi() {
  return dphi;
}

double spherical_grid::get_dtheta() {
  return dtheta;
}

void spherical_grid::update_radiative_cooling(double dt) {

  int k=nr-1;
  
  for(int i=0; i<nphi; i++) {
    for(int j=0; j<ntheta; j++) {

      int index= node_index(k,i,j);
      
      if(vnodes[index].phi<0.5*M_PI || vnodes[index].phi>1.5*M_PI ) {

	//illuminated face is constantly at Teq	
	vnodes[index].T=Teq*pow(fabs(sin(vnodes[index].theta)*cos(vnodes[index].phi)),1.0/4.0);
    
      } else {
	
	//dark side

	//cerr<<vnodes[index].T<<" t*="<<(7.0/3.0)*dr/(radconst*pow(vnodes[index].T,3))<<" Myr"<<endl;
	vnodes[index].T=pow(pow(vnodes[index].T,-1.0/3.0)+3*radconst*dt/dr, -1.0/3.0);
	//cerr<<vnodes[index].T<<endl;

	//int pippo; cin>>pippo;
	
      }
       
    } 
  }
}

void spherical_grid::update_temperature_Euler_adaptive(double& dt) {

  update_radiative_cooling(dt);

  double dTdt_max=0.0;
  
  for(int k=0; k<nr; k++) {
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {

	int index=node_index(k,i,j);

	//Calculate the time derivative of the temperature 
	vnodes[index].dTdt=calculate_time_derivative(k, i, j, index);

	if(fabs(vnodes[index].dTdt) > dTdt_max) {
	 dTdt_max=fabs(vnodes[index].dTdt);
	 dt=0.1*vnodes[index].T/dTdt_max;
       	}
	
      }
    }
  }
  
  if(dt>0.01) dt=0.01;
  
  //update all temperature values
  for(auto index=0; index<vnodes.size(); index++) {
    
    //Euler forumla to calculate temperature at t+dt
    vnodes[index].Tnew=vnodes[index].T+vnodes[index].dTdt*dt;

    if(vnodes[index].Tnew<0) {
      
      cerr<<"------------------------------------------"<<endl;
      cerr<<"*** PROBLEM T<0 @ NODE INDEX: "<<index<<endl;
      cerr<<"Tnew: "<<vnodes[index].Tnew<<" dTdt="<<vnodes[index].dTdt<<" T:"<<vnodes[index].T<<endl;
      
      cerr<<"index: "<<index<<" r="<<vnodes[index].r<<" phi="<<vnodes[index].phi<<
	" theta="<<vnodes[index].theta<<endl;
      cerr<<"d2Tdr2="<<vnodes[index].d2Tdr2<<" "
	  <<" d2Tdphi2="<<vnodes[index].d2Tdphi2<<" "
	  <<" d2Tdtheta2="<<vnodes[index].d2Tdtheta2<<" "
	  <<endl;
      cerr<<"------------------------------------------"<<endl;
      
      int pippo; cin>>pippo;
      
      vnodes[index].Tnew=0;
    }
    
    vnodes[index].T=vnodes[index].Tnew;
  }
  
}



void spherical_grid::update_temperature_midpoint_adaptive(double& dt) {

  update_radiative_cooling(dt);

  //save starting temperatures
  for(auto l=0; l<vnodes.size(); l++) {
    vnodes[l].Told=vnodes[l].T;
  }
  
  double dTdt_max=0.0;
  
  for(int k=0; k<nr; k++) {
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {
	
	int index=node_index(k,i,j);
	
	//Calculate the time derivative of the temperature 
	vnodes[index].dTdt=calculate_time_derivative(k, i, j, index);
	
	if(fabs(vnodes[index].dTdt) > dTdt_max) {
	  dTdt_max=fabs(vnodes[index].dTdt);
	  dt=0.1*vnodes[index].T/dTdt_max;
       	}
	
	
      }
    }
  }

  //if(dt>0.01 || dt<=0) dt=0.01;

  //update all temperature values using Euler forumla to calculate temperature at t+0.5*dt
  for(auto index=0; index<vnodes.size(); index++) {
    vnodes[index].T=vnodes[index].Told+vnodes[index].dTdt*0.5*dt;
  }

  for(int k=0; k<nr; k++) {
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {
	
	int index=node_index(k,i,j);
	
	//Calculate the time derivative of the temperature at midpoint
	vnodes[index].dTdt=calculate_time_derivative(k, i, j, index);	
	
      }
    }
  }


  //update all temperature values
  for(auto index=0; index<vnodes.size(); index++) {
    vnodes[index].Tnew= vnodes[index].Told+vnodes[index].dTdt*dt;
  }
  
}



void spherical_grid::update_temperature_Euler(double dt) {
  update_radiative_cooling(dt);

  for(int k=0; k<nr; k++) {
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {

	int index=node_index(k,i,j);

	//Calculate the time derivative of the temperature 
	double dTdt=calculate_time_derivative(k, i, j, index);
	
	//Euler forumla to calculate temperature at t+dt
	vnodes[index].Tnew=vnodes[index].T+dTdt*dt;

	if(vnodes[index].Tnew<0) {

	  cerr<<"------------------------------------------"<<endl;
	  cerr<<"*** PROBLEM T<0 @ NODE INDEX: "<<index<<endl;
	  cerr<<"Tnew: "<<vnodes[index].Tnew<<" dTdt="<<dTdt<<" T:"<<vnodes[index].T<<endl;

	  cerr<<"index: "<<index<<" r="<<vnodes[index].r<<" phi="<<vnodes[index].phi<<
	    " theta="<<vnodes[index].theta<<endl;
	  cerr<<"d2Tdr2="<<vnodes[index].d2Tdr2<<" "
	      <<" d2Tdphi2="<<vnodes[index].d2Tdphi2<<" "
	      <<" d2Tdtheta2="<<vnodes[index].d2Tdtheta2<<" "
	      <<endl;
	  cerr<<"------------------------------------------"<<endl;

	  int pippo; cin>>pippo;
	  
	  vnodes[index].Tnew=0;
	}
	
	
      }
    }
  }

  //update all temperature values
  for(auto l=0; l<vnodes.size(); l++) {
    vnodes[l].T=vnodes[l].Tnew;
  }
}


void spherical_grid::update_temperature_midpoint( double dt) {
  update_radiative_cooling(dt);

  //save starting temperatures
  for(auto l=0; l<vnodes.size(); l++) {
    vnodes[l].Told=vnodes[l].T;
  }

  for(int k=0; k<nr; k++) {
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {

	int index=node_index(k,i,j);

	//Calculate the time derivative of the temperature 
	double dTdt=calculate_time_derivative(k, i, j, index);
	
	//Euler forumla to calculate midpoint temperature at half dt
	vnodes[index].Tnew=vnodes[index].Told+dTdt*(dt/2);
      }
    }
  }
  
  //update all temperature with midpoint values
  for(auto l=0; l<vnodes.size(); l++) {
    vnodes[l].T=vnodes[l].Tnew;
  }
  
  for(int k=0; k<nr; k++) {
    for(int i=0; i<nphi; i++) {
      for(int j=0; j<ntheta; j++) {

	int index=node_index(k,i,j);

	//Temperature derivative at mid point
	double dTdt_midpoint=calculate_time_derivative(k, i, j, index);
	
	//Euler formula to predict Temperature at the next time step using midpoint temperature derivative
	vnodes[index].Tnew= vnodes[index].Told+dTdt_midpoint*dt;



	if(vnodes[index].Tnew<0) {

	  
	  cerr<<"------------------------------------------"<<endl;
	  cerr<<"*** PROBLEM T<0 @ NODE INDEX: "<<index<<endl;
	  cerr<<"Tnew: "<<vnodes[index].Tnew<<" dTdt_midpoint="<<dTdt_midpoint<<" T:"<<vnodes[index].T<<endl;
	  
	  cerr<<"index: "<<index<<" r="<<vnodes[index].r<<" phi="<<vnodes[index].phi<<
	    " theta="<<vnodes[index].theta<<endl;
	  cerr<<"d2Tdr2="<<vnodes[index].d2Tdr2<<" "
	      <<" d2Tdphi2="<<vnodes[index].d2Tdphi2<<" "
	      <<" d2Tdtheta2="<<vnodes[index].d2Tdtheta2<<" "
	      <<endl;
	  cerr<<"------------------------------------------"<<endl;

	  int pippo; cin>>pippo;

	}

	
      }
    }
  }

  //update all temperature values
  for(auto l=0; l<vnodes.size(); l++) {
    vnodes[l].T=vnodes[l].Tnew;
  }
  
}

void spherical_grid::save(string outfile) {
  
  ofstream fout(outfile);
 
  fout<<get_dr()<<" "<<get_dphi()<<" "<<get_dtheta()<<endl;

  fout<<radius_min<<" "<<radius_max<<" "<<phi_min<<" "<<phi_max<<" "<<theta_min<<" "<<theta_max<<" "<<radconst<<endl;

  fout<<density<<" "<<emissivity<<" "<<Ktherm<<" "<<cp<<" "<<alpha<<endl;
  
  
  for(int j=0; j<ntheta; j++) {
    for(int i=0; i<nphi; i++) {
      for(int k=0; k<nr; k++) {
    
      	int index=node_index(k,i,j);
	
	fout<<index<<"\t"<<vnodes[index].r<<"\t"
	    <<vnodes[index].phi<<"\t"
	    <<vnodes[index].theta<<"\t"
	    <<vnodes[index].T<<endl;
	}
      }
   }

  fout.close();


  cout<<"Saved "<<vnodes.size()<<" grid nodes"<<endl;
}




void spherical_grid::load(string infile) {

  vnodes.clear();
  
  ifstream fin(infile);


  fin>>dr>>dphi>>dtheta;
  fin>>radius_min>>radius_max>>phi_min>>phi_max>>theta_min>>theta_max>>radconst;
  fin>>density>>emissivity>>Ktherm>>cp>>alpha;

  cerr<<dr<<" "<<dphi<<" "<<dtheta<<endl;
  
  cerr<<radius_min<<" "<<radius_max<<" "<<phi_min<<" "<<phi_max<<" "<<theta_min<<" "<<theta_max<<" "<<radconst<<endl;

  int i;
  double r, phi, theta, T;

  while(fin>>i>>r>>phi>>theta>>T) {
    
    spherical_node n;
    n.r=r;
    n.phi=phi;
    n.theta=theta;
    n.T=T;

    vnodes.push_back(n);
  }

  fin.close();
  
  cout<<"Read correctly "<<vnodes.size()<<" nodes"<<endl;  
}



