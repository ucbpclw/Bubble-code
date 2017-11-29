//
//  tissue.hpp
//  Bubble Jan 2014
//
//  Created by claire Walsh on 07/01/2014.
//  Copyright (c) 2014 University college london. All rights reserved.
//

#ifndef Bubble_Jan_2014_tissue_hpp
#define Bubble_Jan_2014_tissue_hpp
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "parameters.h"
#include <string.h>
#include <sstream>
#include "Bubble.hpp"
#include "Dive_profile.hpp"
class Tissue{
private:
    double tissue_size[3];
    int node_number[3];
    double space_step;
    double time_step;
    double*  concN2A;
    double*  concN2B;
    double*  concO2A;
    double*  concO2B;
    double* conc_sN2;
    double* conc_s1N2;
    double* conc_sO2;
    double* conc_s1O2;
    
    uint8* bub_label;
    uint8* temp_bub_label;
    uint8* boundary_label;
    uint8* temp_boundary_label;
    int number_of_bubbles;
    std::vector<Bubble*> theBubbles;
    std::vector<double> total_bubble_mass;
    std::vector<double> total_mass_tissue;
    std::vector<double> total_tissue_volume;
    std::vector<double> tissue_flux;
    
    double total_bubble_mass_pts; //should be the same as total mass bubble (but calculated on the grid.
    
public:
    // constructor
    Tissue(double tissue_size_x,double tissue_size_y,double tissue_size_z,int number_of_bubbles);
    //getters and setters for private members
    double get_tissue_size_x(){return tissue_size[0];}
    double get_tissue_size_y(){return tissue_size[1];}
    double get_tissue_size_z(){return tissue_size[2];}
    void set_tissue_size(double a,double b,double c){tissue_size[0]=a;tissue_size[1]=b;tissue_size[2]=c;}
    int get_node_number_x(){return node_number[0];}
    int get_node_number_y(){return node_number[1];}
    int get_node_number_z(){return node_number[2];}
    void set_node_numbers(int a,int b,int c){node_number[0]=a;node_number[1]=b;node_number[2]=c;}
    double get_space_step(){return space_step;}
    void set_space_step(double a){space_step=a;}
    double get_time_step(){return time_step;}
    void set_time_step(double a){time_step=a;}
    int get_number_of_bubbles(){return number_of_bubbles;}
    void set_number_of_bubbles(int a){number_of_bubbles=a;}
    void set_tissue_mass(double a){total_mass_tissue.push_back(a);}
    void set_tissue_flux(double a){tissue_flux.push_back(a);}
    double get_tissue_flux(){return tissue_flux.back();}
    
    void delete_conc_grid();
    //    double* get_conc(){return conc;}
    //    void set_conc(double*){conc=a;}
    
    // getters and setters for individual bubble private members
    double get_bub_concN2(int i){double bub_concN2 = theBubbles[i]->get_concN2(); return bub_concN2;}   //for the most recent conc.
    double get_bub_concN2(int i,int j){double bub_concN2=theBubbles[i]->get_concN2(j);return bub_concN2;}  // for the conc at a specific time step
    double get_bub_concO2(int i){double bub_concO2 = theBubbles[i]->get_concO2(); return bub_concO2;}   //for the most recent conc.
    double get_bub_concO2(int i,int j){double bub_concO2=theBubbles[i]->get_concO2(j);return bub_concO2;}  // for the conc at a specific
    double get_bub_massN2(int i){double bub_massN2 = theBubbles[i]->get_mmassN2(); return bub_massN2;}   //for the most recent mass.
    double get_bub_massN2(int i,int j){double bub_massN2=theBubbles[i]->get_mmassN2(j);return bub_massN2;}  // for the mass at a specific time step
    double get_bub_massO2(int i){double bub_massO2 = theBubbles[i]->get_mmassO2(); return bub_massO2;}   //for the most recent mass.
    double get_bub_massO2(int i,int j){double bub_massO2=theBubbles[i]->get_mmassO2(j);return bub_massO2;}  // for the mass at a specific time step
    double get_bub_radius(int i){double rad=theBubbles[i]->get_radius(); return rad;}      //for the most recent radius
    double get_bub_radius(int i,int j){double rad=theBubbles[i]->get_radius(j); return rad;} //for a radius at a specific time step
    double get_bub_locx(int i){double xloc=theBubbles[i]->get_locations(0); return xloc;}
    double get_bub_locy(int i){double yloc=theBubbles[i]->get_locations(1); return yloc;}
    double get_bub_locz(int i){double zloc=theBubbles[i]->get_locations(2); return zloc;}
    // void set_bub_press(int i, double a){theBubbles[i]->set_mpress(a);}
    void set_bub_radius(int i, double a){theBubbles[i]->set_radius(a);}
    void set_bub_locs(int i, double a, double b, double c){theBubbles[i]->set_locations(a,b,c);}
    void set_bub_conc(int i,int t){theBubbles[i]->bubbleconc(t);}
    void set_total_bub_mass(double m){total_bubble_mass.push_back(m);}
    
    
    
    // appy bubble functions on each bubble in the tissue block
    void bub_dcdr();          //mass flux
    void bub_drdt(int tt);     //dr/dt
    void bubble_move(int i){theBubbles[i]->move();}   //move a single bubble
    void bub_dc_dr_spherical(int tt);         // applies dc/dr on each bubble.
    void delete_a_Bubble(int b1){theBubbles.erase(theBubbles.begin()+b1);}
    //Tissue functions
    void Check_overlap();                         //overlap between bubbles
    double euclidean_dist(double x1,double y1,double z1,double x2,double y2,double z2);
    void diffusion(int tt);//finite difference tissue
    void check_radius_folder(std::string filepathandname);
    void write_radius_data(std::string,std::string filepathandname);
    void write_conc_data(int tt);
    void initialise_conc_grid();
    void delete_theBubbles();
    void initialise_bubbles(int z);                   // initialise bubbles in grid without overlaps and at min size
    int von_Neumann();                          // stability criterion
    void param_check(int tt);                  // write the parameters to check progress.
    std::vector<double> interpolate(std::vector<double> pts,int bb); //interpolates concentration for any poisition in the grid based on the concentration of the nearest neighbours
    std::vector<std::vector<double > > nearest_neighbours(std::vector<double> pt,double locx,double locy,double locz, double radius); //finds the nearest neighbours to any point in the tissue
    void bubble_check(std::string filepath,int t,int z);  //checks bubble is not too near the edge or overlapping with others needs t and z (which elasticity is used) so that it can pass them to coalescence if two bubbles are found to overlap
    void bub_drdt_alt(int tt); //with Gernhardt elasticity;
    void bub_drdt_gern(int tt); //with Gernhardt elasticity;
    void bub_drdt_lap(int tt); //with no elasticity
    void diffusion_private(const double * c_O2s,double* c_O2s1,const double * c_N2s,double* c_N2s1,double P_amb); //diffusion dan
    std::vector<double> Inverse_dist_interpolate(int i, int j, int k,double* conc_valsN2,double* conc_valsO2); // does inverse distance weighting for boundary points
    void compute_label_bub_total(int tt);
    void compute_label_boundary_total(int tt);
    void move_bubble(int *edge, int b);   //moves bubbles when they are too close to the edge
    void edge_check();                   //checks if a bubble is too near the tissue edge
    void coalescence(std::string filepathname,int b1,int b2,int t,int z); //coalesces the bubbles b1 and b2 creates a new bubble objects and writes the radii of b1 and b2 to file upto timept t, z is needed so the pressure inside the new bubble can be calculated with the correct pressure
    void update_boundary_pts(int tt);
    std::vector<double> mass_conservation(int timept); //Checks conservation of mass
    std::vector<double> mass_bounding_box(int tt,int bb);
    int find_bubble_index(int b); // this function is to find the position of a bubble in the vector of Bubbles given the bubble number. It is needed as if there is coalescence bubble numbers will not necessarily be in order any more. i.e. bubbles may go 1,2,5,6 i.e. at some point 3 and 4 coalesced and became bubble 6.
};

#endif
