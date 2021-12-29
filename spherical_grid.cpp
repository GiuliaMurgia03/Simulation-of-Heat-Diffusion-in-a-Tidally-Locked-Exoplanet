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
  
  albedo=g.albedo;
  density=g.density;
  
  atmosphere_emissivity=g.atmosphere_emissivity;
  planet_emissivity=g.planet_emissivity;

  Ktherm=g.Ktherm;
  cp=g.cp;
  alpha=g.alpha;
  
  Lstar=g.Lstar;
  Dist=g.Dist;
  
  Qdecay0=g.Qdecay;
  Qdecay=g.Qdecay;
  tdecay=g.tdecay;
  
  Teq=g.Teq;
  Tground=g.Tground;
  Tatm=g.Tatm;

  age=g.age;
}

void spherical_grid::set_limits(double r_min, double r_max,
				double p_min, double p_max,
				double t_min, double t_max) {
  radius_min=r_min;
  radius_max=r_max;
  //Deg to Rad
  phi_min=p_min*M_PI/180.0;
  phi_max=p_max*M_PI/180.0;
  theta_min=t_min*M_PI/180.0;
  theta_max=t_max*M_PI/180.0;
}

void spherical_grid::set_time(double t) {

  age=t;
  Qdecay=Qdecay0*exp(-age/tdecay);
}

double spherical_grid::get_time() {
  return age;
}



void spherical_grid::calculate_average_Tground() {

  double sum=0;
  int n=0;
  
  for(int i=0; i<nphi; i++) {
    for(int j=0; j<ntheta; j++) {
      	
      int index=node_index(nr-1,i,j);
      
      sum=sum+vnodes[index].T;
      n=n+1;
    }
  }

  Tground=sum/n;
}


double spherical_grid::get_Tatm() {
  return Tatm;
}


double spherical_grid::get_Tground() {
  return Tground;
}



void spherical_grid::set_physical_properties(double Ls, double D, double d,
					     double planet_e, double atmosphere_e,
					     double K, double c,
					     double a, double qd, double qt) {

  Lstar=Ls*3.828e26;   //Watt
  Dist=D*1.495978707e11; //m

  double solar_constant=(1-albedo)*Lstar/(4*M_PI*pow(Dist,2));
  cout<<"solar_constant="<<solar_constant<<" W/m^2"<<endl;
  
  albedo=a;  
  density=d; //kg/m^3
  planet_emissivity=planet_e; 
  atmosphere_emissivity=atmosphere_e;

  Ktherm=K; // Thermal Conductivity W/m/K
  cp=c; // Heat Capacity J/kg/K
 
  alpha=K/(cp*density);   //Diffusivity coefficient m^2/s
  cout<<"alpha= "<<alpha<<" m^2/s"<<" or "<<alpha*1e6<<" mm^2/s"<<endl;
  
  sigma=5.670367e-8; // W/m^2/K^4

  Teq=pow((1-albedo)*Lstar/(4*sigma*4*M_PI*pow(Dist,2)),1.0/4.0);
  
  cout<<"Teq="<<Teq<<" K"<<endl;


  Tatm=pow(2-atmosphere_emissivity,-1.0/4.0)*Teq;

  tdecay=qt*1000;  //Gyr to Myr
  Qdecay0=qd*1e-6*3.1536e13/(cp*density);   // K/Myr
}



  
void spherical_grid::calculate_radial_derivatives(int k, int i, int j) {

  int index=node_index(k,i,j);
  double Tcenter= vnodes[index].T;
  
  int index_after, index_before;
  double Tafter, Tbefore;
  
  //Boundary conditions
  
  if(k==0) { 
    index_before=node_index(k+1,i,j);//Internal Boundary condition: no heat flux, i.e. the first derivative is equal to 0
    index_after=node_index(k+1,i,j);
    Tbefore=vnodes[index_before].T;
    Tafter=vnodes[index_after].T;
  }
  
  else if (k==(nr-1)){ //Surface Boundary condition
    
    index_before=node_index(k-1,i,j);
    
    Tbefore=vnodes[index_before].T;
    
    if(vnodes[index].phi<0.5*M_PI || vnodes[index].phi>1.5*M_PI ) {  //Illuminated Side

       Tafter=Tbefore+2*dr*1e3*(sigma/Ktherm)*
	 (4*fabs(sin(vnodes[index].theta)*cos(vnodes[index].phi))*pow(Teq,4.0)
	  +atmosphere_emissivity*pow(Tatm,4)-planet_emissivity*pow(Tcenter,4.0)); 
    }

    else { //Dark Side
      Tafter=Tbefore+2*dr*1e3*(sigma/Ktherm)
	*(atmosphere_emissivity*pow(Tatm,4)-planet_emissivity*pow(Tcenter,4.0));
    }
  }

  else {
    index_before=node_index(k-1,i,j);
    index_after=node_index(k+1,i,j);
    
    Tbefore=vnodes[index_before].T;
    Tafter=vnodes[index_after].T;
  }
  
  
  //Calculate the second derivative of the temperature with the respect to radius r (d2Tdr2)
  double r=vnodes[index].r;
  vnodes[index].d2Tdr2=(Tafter-Tbefore)/(r*dr)+(Tbefore+Tafter-2*Tcenter)/(dr*dr);  
}
 


void spherical_grid::calculate_longitudinal_derivatives(int k, int i, int j) {

	int index=node_index(k,i,j);
	int index_before;
	int index_after;

	//Periodic boundary conditions
     
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
	
	//Calculate the second derivative of the temperature with the respect to the angle phi (d2Tdphi2)
	
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

	//Boundary conditions
      
	if(j==0) {
	  index_before=node_index(k,i,j+1); //External Boundary condition (North Pole): no heat flux, i.e. the first derivative is equal to 0
	  index_after=node_index(k,i,j+1);
	}

	else if (j==(ntheta-1)){
	  index_before=node_index(k,i,j-1); //External Boundary condition (South Pole): no heat flux, i.e. the first derivative is equal to 0
	  index_after=node_index(k,i,j-1);
	}

	else {
	  index_before=node_index(k,i,j-1);
	  index_after=node_index(k,i,j+1);
	}

	
	//Calculate the second derivative of the temperature with the respect to theta r (d2Tdtheta2)

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
	n.T=T0;	
	
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

  //Calculate radial resolution
  dr=(radius_max-radius_min)/nr;
  
  cout<<"radius min="<<radius_min<<" max="<<radius_max<<" dr: "<<dr<<endl;

  //Calculate angular longitude resolution
  dphi=2*M_PI/nphi; //rad
  cout<<"phi min="<<phi_min*180/M_PI<<" max="<<phi_max*180/M_PI<<" dphi: "<<dphi*180/M_PI<<"°"<<endl;
  
  //Calculate angular latitude resolution
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

void spherical_grid::update_temperature_Euler_adaptive(double& dt) {
  
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
  
  if(dt>0.01) dt=0.01; //Maximum dt allowed 0.01 Myr
  
  //Update all temperature values
  for(int index=0; index<vnodes.size(); index++) {
    
    //Euler forumla to calculate temperature at t+dt
    vnodes[index].Tnew=vnodes[index].T+vnodes[index].dTdt*dt;

    if(vnodes[index].Tnew<0) {
      
      cerr<<"The simulation fails due to negative temperature."<<endl;
      exit(EXIT_FAILURE);
    }
    
    vnodes[index].T=vnodes[index].Tnew;
  }


  //Update Tatm
  calculate_average_Tground();
  Tatm=Tground*pow(planet_emissivity/2.0,1.0/4.0);
}



void spherical_grid::save(string outfile) {
  
  ofstream fout(outfile);
  
  fout<<age<<endl;
  fout<<Tatm<<endl;
  fout<<Tground<<endl;
  fout<<Teq<<endl;

  fout<<get_dr()<<" "<<get_dphi()<<" "<<get_dtheta()<<endl;

  fout<<radius_min<<" "<<radius_max<<" "<<phi_min<<" "<<phi_max<<" "<<theta_min<<" "<<theta_max<<endl;

  fout<<density<<" "<<planet_emissivity<<" "<<atmosphere_emissivity<<" "<<Ktherm<<" "<<cp<<" "<<alpha<<endl;
  
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

  calculate_average_Tground();

  fin>>age;
  fin>>Tatm;
  fin>>Tground;
  fin>>Teq;
  fin>>dr>>dphi>>dtheta;
  fin>>radius_min>>radius_max>>phi_min>>phi_max>>theta_min>>theta_max;
  fin>>density>>planet_emissivity>>atmosphere_emissivity>>Ktherm>>cp>>alpha;

  cerr<<dr<<" "<<dphi<<" "<<dtheta<<endl;
  
  cerr<<radius_min<<" "<<radius_max<<" "<<phi_min<<" "<<phi_max<<" "<<theta_min<<" "<<theta_max<<endl;

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



