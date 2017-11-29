//
//  Bubble.hpp
//  Bubble Jan 2014
//
//  Created by claire Walsh on 07/01/2014.
//  Copyright (c) 2014 University college london. All rights reserved.
//

#ifndef Bubble_Jan_2014_Bubble_hpp
#define Bubble_Jan_2014_Bubble_hpp
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "parameters.h"
#include <string.h>
#include <sstream>
#include "Dive_profile.hpp"

class Bubble{
private:
    int bubble_number;
    double locations[3];
    std::vector<double> radius;
    // std::vector<double> bubpress;
    std::vector<double> bubconcO2;
    std::vector<double> bubconcN2;
    std::vector<double> dc_dr;
    std::vector<double> dc_drN2;
    std::vector<double> dc_drO2;
    double tauvar;
    std::vector<double> bubmassN2;
    std::vector<double> bubmassO2;
    
    
public:
    //constructors
    Bubble(){};
    Bubble(int a){bubble_number=a;};
    Bubble(int a,int lx,int ly,int lz,double rad,int t,double N2mass,double O2mass, double N2conc, double O2conc){bubble_number=a,locations[0]=lx,locations[1]=ly,locations[2]=lz,radius.assign(t-1,0),radius.push_back(rad),/*bubpress.assign(t,0);*/bubmassN2.assign(t-1,0), bubmassN2.push_back(N2mass), bubmassO2.assign(t-1,0), bubmassO2.push_back(O2mass),bubconcN2.assign(t-1,0), bubconcN2.push_back(N2conc),bubconcO2.assign(t-1,0), bubconcO2.push_back(O2conc);}
    
    //getters and setters
    void set_bubble_number(int n){bubble_number=n;}
    void set_radius(double n){radius.push_back(n);}
    //  void set_mpress(double n){bubpress.push_back(n);}
    void set_locations(double nx, double ny, double nz){
        locations[0]=nx;
        locations[1]=ny;
        locations[2]=nz;}
    void set_dcdrN2(double dcdr){dc_drN2.push_back(dcdr);}
    void set_dcdrO2(double dcdr){dc_drO2.push_back(dcdr);}
    void set_tauvar(double n){tauvar=n;}
    int get_bubble_number(){return bubble_number;}
    double get_locations(int coord){return locations[coord];}
    double get_radius(){return radius.back();}
    double get_radius(int i){return radius[i];}
    //  double get_mpress(){return bubpress.back();}
    // double get_mpress(int i){return bubpress[i];}
    double get_tauvar(){return tauvar;}
    double get_mmassN2(){return bubmassN2.back();}
    double get_mmassO2(){return bubmassO2.back();}
    double get_mmassN2(int i){return bubmassN2[i];}
    double get_mmassO2(int i){return bubmassO2[i];}
    double get_concN2(){return bubconcN2.back();}
    double get_concO2(){return bubconcO2.back();}
    double get_concN2(int i){return bubconcN2[i];}
    double get_concO2(int i){return bubconcO2[i];}
    double get_dcdr(){return (dc_drN2.back()+dc_drO2.back());}
    double get_dcdrN2(){return dc_drN2.back();}
    double get_dcdrO2(){return dc_drO2.back();}
    //
    //initialise the random locations of the bubbles
    void randomise_bubble(int node_number[3]);
    //shift a bubble
    void move();
    //caluculate the pressure in a bubble given the radius and external pressure
    void pressure_inside(int tt);
    void bubbleconc(int tt);
    //rounding function
    int round(double x) {return (int)floor(x+0.5);}
    //calculates dc/dr
    void dcdr(double* conc, double space_step);
    //calculates dr/dt
    void drdt(int t, double timept);
    // diff function calculating the increase in mass flux due to an increase in radius
    double diff(double r){return diffusivity*(1-pow((radius[0]/r),2));}
    //for calculating the conc of spherical points
    //transform to spherical polars
    std::vector<std::vector<double> > spherical_tansform_r1(double tissue_step,double angle);
    std::vector<std::vector<double> > spherical_tansform_r2(double deltar,double angle);
    //get the concentrations of these new points by interplotation
    std::vector<std::vector<double> > interpolate(std::vector<std::vector<std::vector<std::vector <double> > > > conc,std::vector<double> pt);
    //use the new points to get dc/dr
    double dcdr_sphere(std::vector<double>newconcvec_r1_N2,std::vector<double>newconcvec_r1_O2,std::vector<double>newconcvec_r2_N2,std::vector<double>newconcvec_r2_O2,double dr,int tt);    // dr/dt using Gernhardt elasticity
    void drdt_gern(int t,double timestep);
    // dr/dt using lapace
    void drdt_laplace(int t,double timestep);
    //Calculate the mass of gas in the bubble
    double mass_in_bubble(int tt,double timestep);
    // Variable permeability
    void var_perm();
    // bubble 1 zero array
    void compute_bubble_label(uint8*, int tt);
    void compute_boundary_label(uint8*,uint8*, int tt);
    // calculate the distance between a given point (xx,yy,zz) and the bubble edge.
    double surf_to_point(int xx, int yy, int zz);
    double surf_to_point(int xx, int yy, int zz,int tt);
    // caluclates the intersections between a bubble and any line described by the points LinePoint0 and line point1
    std::vector<double> LineSphereIntersections( std::vector<double> linePoint0,  std::vector<double> linePoint1);
    double length_line_square(double theta, double phi,int face);
    double xdist(int i, int j, int k);
    double ydist(int i, int j, int k);
    double zdist(int i, int j, int k);
    int which_face(double i, double j, double k);
};


#endif
