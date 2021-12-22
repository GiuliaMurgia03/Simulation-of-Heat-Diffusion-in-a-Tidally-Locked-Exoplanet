// g++ -o planet_graphics.exe planet_graphics.cpp spherical_grid.cpp spherical_node.cpp -lGL -lGLU -lglut


//  g++ -w -o planet_graphics.exe planet_graphics.cpp spherical_grid.cpp  spherical_node.cpp -L/System/Library/Frameworks -framework GLUT -framework OpenGL


#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#else
#include <GL/glut.h>
#endif


#include <GL/freeglut.h>

#include<cmath>
#include<iostream>
#include <time.h>
#include<fstream>
#include<vector>
#include<string>
#include <iomanip>
#include <sstream>
#include"spherical_grid.h"

using namespace std;
using namespace planet_code;

double user_phi=0;
double user_theta=0;

double user_height=-6000;
double user_xcenter=0;
double user_ycenter=0;

float xRotated = -45.0, yRotated = 0.0, zRotated = 240.0;

bool user_rotate=false;

spherical_grid grid;

double Tmin, Tmax;

void reshape (int w, int h) {
  
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();

    //https://stackoverflow.com/questions/16571981/gluperspective-parameters-what-do-they-mean
    
    gluPerspective (25.0, (GLdouble)w / (GLdouble)h, 0.5, 12000);
    glMatrixMode   (GL_MODELVIEW);
    glViewport     (0, 0, w, h);
}

void process_key(unsigned char key, int x, int y) {

  if (key == '+') user_height  += 50;  
  if (key == '-') user_height  -= 50;

  if (key == 'r') {

    if(user_rotate){
      user_rotate=false;
    } else {
      user_rotate=true;
    }
  }
  
  glutPostRedisplay();
}

void special_key(int k, int x, int y) {

  switch(k) {

  case GLUT_KEY_F1:  xRotated += 1; break;
  case GLUT_KEY_F2:  xRotated -= 1; break;
  case GLUT_KEY_F3:  yRotated -= 1; break;
  case GLUT_KEY_F4:  yRotated += 1; break;
  case GLUT_KEY_F5:  zRotated -= 1; break;
  case GLUT_KEY_F6:  zRotated += 1; break;

  case GLUT_KEY_UP:  user_ycenter  += 10; break;
  case GLUT_KEY_DOWN:  user_ycenter  -= 10; break;
  case GLUT_KEY_LEFT:  user_xcenter  += 10; break;
  case GLUT_KEY_RIGHT:  user_xcenter  -= 10; break;

  }
  
  glutPostRedisplay();
}



void temp_to_color(double temp, GLdouble& r, GLdouble& g, GLdouble& b) {


  double t1=0;
  double t2=150;
  
  if(temp>=t1 && temp<t2) {       // dark gray to white

    r=0.2;
    g=0.2;
    b=0.2;
    
    return;
  }
    

  t1=150;
  t2=273;
  
  if(temp>=t1 && temp<t2) {       // dark gray to white

    r=0.2+0.8*(temp-t1)/(t2-t1);
    g=0.2+0.8*(temp-t1)/(t2-t1);
    b=0.2+0.8*(temp-t1)/(t2-t1);
    
    return;
  }

  t1=273;
  t2=293;
  
  if(temp>=t1 && temp<t2) {     //cyan
    r=1-(temp-t1)/(t2-t1);
    g=1;
    b=1;
    
    return;
  }
  
  t1=293;
  t2=313;
  
  if(temp>=t1 && temp<t2) {     //blue
    r=0;
    g=1-(temp-t1)/(t2-t1);
    b=1;
    
    return;
  }

  t1=313;
  t2=333;
  
  if(temp>=t1 && temp<t2) {     //brown
    r=0.5*(temp-t1)/(t2-t1);
    g=0.35*(temp-t1)/(t2-t1);
    b=1-(temp-t1)/(t2-t1);
    
    return;
  }
  
  t1=333;
  t2=750;
  
  if(temp>=t1 && temp<t2) {     //brown to yellow 
    r=0.5+0.5*(temp-t1)/(t2-t1);
    g=0.35+0.65*(temp-t1)/(t2-t1);
    b=0;
    
    return;
  }


  t1=750;
  t2=1000;
  
  if(temp>=t1 && temp<t2) {     //yellow to orange
    r=1;
    g=1-0.5*(temp-t1)/(t2-t1);
    b=0;
    
    return;
  }

  t1=1000;
  t2=1250;
  
  if(temp>=t1 && temp<t2) {     //orange to red
    r=1;
    g=0.5-0.5*(temp-t1)/(t2-t1);
    b=0;
    
    return;
  }

  t1=1250;
  t2=1500;
  
  if(temp>=t1 && temp<t2) {     //red to bordeaux
    r=1-0.25*(temp-t1)/(t2-t1);
    g=0.0;
    b=0.1*(temp-t1)/(t2-t1);
    
    return;
  }
  
  r=0.75;
  g=0.0;
  b=0.1;
    
  return;  
}

string double_to_string(double d, int precision) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) <<d;
  return ss.str();
}


void draw_planet (void) {

  //colormap
  glMatrixMode(GL_MODELVIEW);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  
  glEnable(GL_DEPTH_TEST);
  
  glTranslatef(user_xcenter, user_ycenter, user_height); 
  glPushMatrix();

  //atmosphere
  int circle_points = 100;
  float angle = 2.0f * 3.1416f / circle_points;

  GLdouble red, green, blue;
  temp_to_color(grid.get_Tatm(), red, green, blue);

  for(int radius=1100; radius>=1000; radius--) {

    double c=0.75*(1100-radius)/100.0;
    
    glColor3f(red*c, green*c, blue*c);
    
    glBegin(GL_POLYGON);
    double angle1=0.0;
    glVertex2d(radius * cos(0.0) , radius * sin(0.0));

    for (int i=0; i<circle_points; i++) {       
      glVertex2d(radius * cos(angle1), radius *sin(angle1));
      angle1 += angle;
    }
    
    glEnd();
  }
  
  glFlush();

  

  for(auto i=0; i<1024; i++) {

    double T=Tmin+i*(Tmax-Tmin)/1023.0;
    
    GLdouble red, green, blue;
    temp_to_color(T, red, green, blue);
    glColor3f(red, green, blue);                    
    glBegin(GL_LINES);
    glVertex3f(-512+i, -1150, 0);
    glVertex3f(-512+i, -1250, 0);
    glEnd();
  }
  
  glColor3f(1, 1, 1);                    

  /*
  glBegin(GL_LINE_LOOP);
  glVertex3f(-514, -1149, 0);
  glVertex3f(-514+1028, -1149, 0);
  glVertex3f(-514+1028, -1251, 0);
  glVertex3f(-514, -1251, 0);
  glEnd();
  */
  
  
  std::string tmax = double_to_string(Tmax,0)+" K";
  std::string tmin = double_to_string(Tmin,0)+" K";
  std::string age = double_to_string(grid.get_time()/1000.0,2)+" Gyr";
  
  
  glRasterPos2f(512+50, -1225);
  const unsigned char* t = reinterpret_cast<const unsigned char *>(tmax.c_str());
  glutBitmapString(GLUT_BITMAP_HELVETICA_18, t);
  
  glRasterPos2f(-512-175, -1225);
  t = reinterpret_cast<const unsigned char *>(tmin.c_str());
  glutBitmapString(GLUT_BITMAP_HELVETICA_18, t);

  glRasterPos2f(-100, 1225);
  t = reinterpret_cast<const unsigned char *>(age.c_str());
  glutBitmapString(GLUT_BITMAP_HELVETICA_18, t);
 

  
  glPopMatrix();
  glFlush();

  
  
  //planet
   
  //glMatrixMode(GL_MODELVIEW);
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  //glEnable(GL_DEPTH_TEST);

  glTranslatef(user_xcenter, user_ycenter, user_height); 
  glPushMatrix();


  glRotatef(xRotated, 1, 0, 0);
  glRotatef(yRotated, 0, 1, 0);
  glRotatef(zRotated, 0, 0, 1);

  glEnable(GL_DEPTH_TEST);

  double hdr=0.5*grid.get_dr();
  double hdphi=0.5*grid.get_dphi();
  double hdtheta=0.5*grid.get_dtheta();

  double rfac=1000.0/grid.get_radius();
  
  for(int k=0; k<grid.vnodes.size() ; k++){
  
    GLdouble r0=(grid.vnodes[k].r-hdr)*rfac;
    GLdouble r1=(grid.vnodes[k].r+hdr)*rfac;

    GLdouble phi0=grid.vnodes[k].phi-hdphi;
    GLdouble phi1=grid.vnodes[k].phi+hdphi;

    GLdouble theta0=grid.vnodes[k].theta-hdtheta;
    GLdouble theta1=grid.vnodes[k].theta+hdtheta;

    GLdouble sinp0=sin(phi0);
    GLdouble cosp0=cos(phi0);
  
    GLdouble sint0=sin(theta0);
    GLdouble cost0=cos(theta0);
    
    GLdouble sinp1=sin(phi1);
    GLdouble cosp1=cos(phi1);
  
    GLdouble sint1=sin(theta1);
    GLdouble cost1=cos(theta1);
   
    GLdouble T=grid.vnodes[k].T;

    bool plot=true;
    if(grid.vnodes[k].phi>0 && grid.vnodes[k].phi<=M_PI/2 && grid.vnodes[k].theta < 0.5*M_PI) plot=false;
    
    
    if(plot) {

      //top face - External
	
      //bottom left corner
      GLdouble tx1=r1*cosp0*sint0;
      GLdouble ty1=r1*sinp0*sint0;
      GLdouble tz1=r1*cost0;
      
      //bottom right corner
      GLdouble tx2=r1*cosp1*sint0;
      GLdouble ty2=r1*sinp1*sint0;
      GLdouble tz2=r1*cost0;
      
      //top right corner
      GLdouble tx3=r1*cosp1*sint1;
      GLdouble ty3=r1*sinp1*sint1;
      GLdouble tz3=r1*cos(theta1);
      
      //top left corner
      GLdouble tx4=r1*cosp0*sint1;
      GLdouble ty4=r1*sinp0*sint1;
      GLdouble tz4=r1*cos(theta1);
      
      //bottom face - Internal
      //blc
      GLdouble bx1=r0*cosp0*sint0;
      GLdouble by1=r0*sinp0*sint0;
      GLdouble bz1=r0*cost0;
      
      //brc
      GLdouble bx2=r0*cosp1*sint0;
      GLdouble by2=r0*sinp1*sint0;
      GLdouble bz2=r0*cost0;
      
      //trc
      GLdouble bx3=r0*cosp1*sint1;
      GLdouble by3=r0*sinp1*sint1;
      GLdouble bz3=r0*cos(theta1);
      
      //tlc
      GLdouble bx4=r0*cosp0*sint1;
      GLdouble by4=r0*sinp0*sint1;
      GLdouble bz4=r0*cos(theta1);
      
      glDepthFunc(GL_LEQUAL);
      
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      
      GLdouble red, green, blue;
      temp_to_color(T, red, green, blue);
      glColor3f(red, green, blue);                    
      
      //top face (counter-clockwise)
      //glColor3f(1, 0, 0);
      glBegin(GL_POLYGON);
      glVertex3f(tx1, ty1, tz1);	
      glVertex3f(tx2, ty2, tz2);
      glVertex3f(tx3, ty3, tz3);
      glVertex3f(tx4, ty4, tz4);
      //glNormal3f((tx4-bx4)/dr, (ty4-by4)/dr, (tz4-bz4)/dr);
      glEnd();

      //bottom face
      //glColor3f(0, 1, 0);                    
      glBegin(GL_POLYGON);
      glVertex3f(bx1, by1, bz1);
      glVertex3f(bx2, by2, bz2);
      glVertex3f(bx3, by3, bz3);
      glVertex3f(bx4, by4, bz4);      
      glEnd();
      
      //north face
      //glColor3f(0, 0, 1);                    
      glBegin(GL_POLYGON);
      glVertex3f(bx1, by1, bz1);
      glVertex3f(bx2, by2, bz2);
      glVertex3f(tx2, ty2, tz2);
      glVertex3f(tx1, ty1, tz1);    		   	
      glEnd();
      
      //south face
      //glColor3f(0, 1, 1);                    
      glBegin(GL_POLYGON);
      glVertex3f(bx3, by3, bz3);
      glVertex3f(bx4, by4, bz4);
      glVertex3f(tx4, ty4, tz4);
      glVertex3f(tx3, ty3, tz3);
      glEnd();
      
      //east face
      //glColor3f(1, 0, 1);                    
      glBegin(GL_POLYGON);
      glVertex3f(tx2, ty2, tz2);
      glVertex3f(bx2, by2, bz2);
      glVertex3f(bx3, by3, bz3);
      glVertex3f(tx3, ty3, tz3);
      glEnd();
      
      //west face
      //glColor3f(1, 1, 0);                    
      glBegin(GL_POLYGON);
      glVertex3f(bx1, by1, bz1);
      glVertex3f(tx1, ty1, tz1);
      glVertex3f(tx4, ty4, tz4);
      glVertex3f(bx4, by4, bz4);
      glEnd();
      
      
      //grid
      
      glColor3f(0.15, 0.15, 0.15);
      glLineWidth(1.25);
      glEnable(GL_LINE_SMOOTH);
      glDepthFunc(GL_LEQUAL);
      
      
      //top 
      glBegin(GL_LINE_LOOP);
      glVertex3f(tx1, ty1, tz1);
      glVertex3f(tx2, ty2, tz2);
      glVertex3f(tx3, ty3, tz3);
      glVertex3f(tx4, ty4, tz4);
      glEnd();	
      
      //bottom 
      glBegin(GL_LINE_LOOP);
      glVertex3f(bx1, by1, bz1);
      glVertex3f(bx2, by2, bz2);
      glVertex3f(bx3, by3, bz3);
      glVertex3f(bx4, by4, bz4);      
      glEnd();
      
      //north
      glBegin(GL_LINE_LOOP);
      glVertex3f(bx1, by1, bz1);
      glVertex3f(bx2, by2, bz2);
      glVertex3f(tx2, ty2, tz2);
      glVertex3f(tx1, ty1, tz1);
      glEnd();
      
      //south 
      glBegin(GL_LINE_LOOP);
      glVertex3f(bx3, by3, bz3);
      glVertex3f(bx4, by4, bz4);
      glVertex3f(tx4, ty4, tz4);
      glVertex3f(tx3, ty3, tz3);
      glEnd();
      
      //east 
      glBegin(GL_LINE_LOOP);
      glVertex3f(tx2, ty2, tz2);
      glVertex3f(bx2, by2, bz2);
      glVertex3f(bx3, by3, bz3);
      glVertex3f(tx3, ty3, tz3);
      glEnd();
      
      //west 
      glBegin(GL_LINE_LOOP);
      glVertex3f(bx1, by1, bz1);
      glVertex3f(tx1, ty1, tz1);
      glVertex3f(tx4, ty4, tz4);
      glVertex3f(bx4, by4, bz4);
      glEnd();
      
    }
  }
  
  glPopMatrix ();
  glutSwapBuffers();
  glFlush();
}

 
 
void idleFunc (void) {

  if(user_rotate) {
    zRotated += 0.5; 
    glutPostRedisplay();
  }

  glFlush ();
}


int main (int argc, char **argv) {

  //load grid
  if(argc>1) {
    grid.load(argv[1]);
    Tmin=grid.get_Tmin();
    Tmax=grid.get_Tmax();
    cout<<"Grid Tmin="<<Tmin<<" Tmax="<<Tmax<<endl;

    Tmin=0;
    Tmax=1500;
    
    cout<<"Using Tmin="<<Tmin<<" Tmax="<<Tmax<<endl;
  }

  
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE);
  
  glutInitWindowSize(720, 720);
  glutInitWindowPosition(10, 10);
  glutCreateWindow("PLANET CODE");
  
  glClearColor(0.0, 0.0, 0.0, 0.0);
    
  glutDisplayFunc(draw_planet);
  glutReshapeFunc(reshape);
  
  glutIdleFunc(idleFunc);
  
  glutSpecialFunc(special_key);
  glutKeyboardFunc(process_key);
  
  glutMainLoop();

  return 0;
}
