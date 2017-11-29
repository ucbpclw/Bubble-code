//
//  parameters.h
//  Bubble Jan 2014
//
//  Created by claire Walsh on 07/01/2014.
//  Copyright (c) 2014 University college london. All rights reserved.
//

#ifndef Bubble_Jan_2014_parameters_h
#define Bubble_Jan_2014_parameters_h
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#define node_num_x 128   //this a power of 2 so 2,4,8,16,32,64,128,256,512,1024,2048
#define node_num_y 128   // this is a power of 2 also
#define node_num_z 32    //this is a power of 2 also
#define bit_shift_i 7 //should be log2(node_num_x)
#define bit_shift_ij 14 // this needs to be log2(node_num_x * node_num_y) or bitshift_i +bitshift_y
#define U(i,j,k) ((i)|(j)<<bit_shift_i|(k)<<bit_shift_ij)
typedef unsigned char uint8;

const double rbar=1e-4; // non-dimensionalising length in (m),this is also the tissue step
const double pbar=1.0; // non- dimensioalising pressure in bar (atmospheric pressure)
const double RT=8.3144598*310;  // universal gas constant 310K (37C) taken from engineering toolbox // http://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html mol/Nm.
const double MrN2=28;
const double MrO2=32;
const double RT_N2=0.9201;
const double RT_O2=0.8054; // specific gas constant for air at 310K (37C) taken from engineering toolbox // http://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html and converted to m3bar/kg
const double L_N2=0.014;//0.0436;; // Otswald number dimesionless k_h=L/RT set to 0.145 from Lango
const double L_O2=0.027;//0.027;  //set to 0.027 from lango
const double k_h_N2= (L_N2/RT_N2); //Henry's const in kgm^-3bar^-1 from Chappell thesis
const double k_h_O2= (L_O2/RT_O2); //Henry's const in kgm^-3bar^-1 from Chappell thesis
const double cbar= (k_h_N2*pbar); // non-dimensionalising conc. needs to be in kg/m3
const double mbar=k_h_N2*pbar*rbar*rbar*rbar;
const double diffusivity=2.5e-9; // diffusion coefficient m^2/s from Srinivasan et. al (2.2e-12). NB changed to 2.2e-10 on the 19/3/14.changed to 1e-10 on 6/9/14 on experimental data basis
const double diffusivity2=diffusivity/30;; // diffusion coefficient m^2/s from Srinivasan et. al (2.2e-12). NB changed to 2.2e-10 on the 19/3/14.changed to 1e-10 on 6/9/14 on experimental data basis
const double RTbar=RT/(pow(rbar,3)*pbar);
const double pp_fraction_N2=0.8;
const double pp_fraction_O2=0.2;
const double tau=(k_h_N2*diffusivity*RT_N2*RT_O2)/((pow(rbar,2)*((pp_fraction_N2*RT_O2)+(pp_fraction_O2*RT_N2)))); //secs non-dimensional parameters =rbar^2/LD
const double DC_DR_SENS= 0.25;//Distance to next node below which next pointis used in calculating dc_dr
const double mu=40e-5/pbar; // non-dimensional small shear modulus in (bar) taken from Darrell Velegol,Frederick Lanni (2001) 54Pa
const double PI=4*atan(1);
const double sigma=0.07e-5/(rbar*pbar); // non-dimensional surface tension taken from Zhou et. al (it is the surface tension of water (0.07) bodily should be lower more like 0.04 converted to m bar)
const double initial_radius=1;  //non dim initial bubble size
const double min_radius=1; //min radius set as 2x rbar;
const double time_of_dive=6740.0; //time of dive in mins. can be in seconds but need to sort out profile.cpp and change the time loops in the main and radius write files.
const double M= /*5.06625e12*/ 0*(pow(rbar,3)/pbar); //bulk modulus over affected volume, (bar/m^3) given in Srinivasan as 5.0e-6atm/um^3.
const double diffusion_dist=1;   //dc distance
const double diff_region_prop=500.0;  //how thick the diffusion region is as a function of radius.

//********************************************************************
//**************                              ************************
//**************    SQUARED                    ************************
//**************                              ************************
//********************************************************************
inline double SQUARED(double a){
    //TODO: I think you can do this bitwise, which the compiler may not realise..though bitwise on double is a hack
    return a*a;
}

#endif
