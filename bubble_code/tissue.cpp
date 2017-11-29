//
//  tissue.cpp
//  Bubble Jan 2014
//
//  Created by claire Walsh on 07/01/2014.
//  Copyright (c) 2014 University college london. All rights reserved.
//

#include <iostream>
#include <iostream>
#include <string>
#include <cassert>
#include <math.h>
#include <ctime>
#include <numeric>
#include <cstdlib>
#include "parameters.h"
#include <vector>
#include "Dive_profile.hpp"
#include <stdio.h>
#include <stdlib.h>
#include "tissue.hpp"
//#include "time.h"
//********************************************************************
//**************                              ************************
//**************   construct tissue block     ************************
//**************                              ************************
//********************************************************************

Tissue::Tissue(double size_x,double size_y,double size_z,int bub_num){
    tissue_size[0]=(size_x);
    tissue_size[1]=(size_y);
    tissue_size[2]=(size_z);
    node_number[0]=(int)(size_x/rbar);
    node_number[1]=(int)(size_y/rbar);
    node_number[2]=(int)(size_z/rbar);
    number_of_bubbles=bub_num;
    space_step=rbar;
    
    for(int i=1;i<bub_num+1;i++){ //not this means that bubble number and therefore bubble number in the masks is always 1 higher than the vector position in theBubbles, i.e. bubble 1 is accessed via theBubbles[0];f
        theBubbles.push_back(new Bubble(i));
    }
}


//********************************************************************
//**************                              ************************
//**************   initialise conc grid       ************************
//**************                              ************************
//********************************************************************


void Tissue::initialise_conc_grid(){
    // inititalise the conc. vector, allocating space for the grid
    concN2A=new double[node_number[0]*node_number[1]*node_number[2]];
    concN2B=new double[node_number[0]*node_number[1]*node_number[2]];
    concO2A=new double[node_number[0]*node_number[1]*node_number[2]];
    concO2B=new double[node_number[0]*node_number[1]*node_number[2]];
    
    
    bub_label=new uint8[node_number[0]*node_number[1]*node_number[2]];
    boundary_label=new uint8[node_number[0]*node_number[1]*node_number[2]];
    conc_sN2=concN2A;
    conc_s1N2=concN2B;
    conc_sO2=concO2A;
    conc_s1O2=concO2B;
    
    
    for(int i=0; i<node_number[0]; i++){
        for (int j=0; j<node_number[1]; j++){
            for (int k =0; k<node_number[2]; k++){        // initially all points will be set to initial external pressure other than those within the bubbles
                conc_sN2[U(i,j,k)]=pp_fraction_N2*P_amb[0]; //P_amb is non dimesnional
                conc_sO2[U(i,j,k)]=(k_h_O2/k_h_N2)*pp_fraction_O2*P_amb[0]; //use memset?? sort out time steps 0 and 1
                
                
                
                for(int bb=0; bb<number_of_bubbles; bb++){
                    if (euclidean_dist(i,j,k,get_bub_locx(bb),get_bub_locy(bb),get_bub_locz(bb))<=get_bub_radius(bb)){
                        conc_sN2[U(i,j,k)]=get_bub_concN2(bb,0);
                        conc_sO2[U(i,j,k)]=get_bub_concO2(bb,0);
                        //  std::cout<<conc_sN2[U(i,j,k)]<<" N2 bubble conc \n"<<conc_sO2[U(i,j,k)]<<" O2 bubble conc\n";
                    }
                }
            }
        }
    }
    
    std::cout<<conc_sO2[U(1,1,1)]<<" O2 inital conc tissue (1,1,1) \n"<<conc_sN2[U(1,1,1)]<<" N2 inital conc tissue (1,1,1) \n";
    //  write_conc_data(0);
    
}

//********************************************************************
//**************                              ************************
//**************   delete   conc grid       ************************
//**************                              ************************
//********************************************************************
void Tissue::delete_conc_grid(){
    delete concN2A;
    delete concN2B;
    delete concO2A;
    delete concO2B;
    delete bub_label;
    delete boundary_label;
}

//********************************************************************
//**************                              ************************
//**************      delete the bubbles      ************************
//**************                              ************************
//********************************************************************
void Tissue::delete_theBubbles(){
    for (int i=0; i<number_of_bubbles; i++){
        delete theBubbles[i];            //deletes the bubble vector memory.
    }
}
//********************************************************************
//**************                              ************************
//**************      find bubble index      ************************
//**************                              ************************
//********************************************************************
int Tissue::find_bubble_index(int b_given){
    int bubble_index;
    for(int bb=0; bb<number_of_bubbles; bb++){
        if(theBubbles[bb]->get_bubble_number()==b_given){
            bubble_index=bb;
            break;}
        //elseif{ std::cout<<"that is not a bubble you have a problem\n";
        //  assert(10==1);
    }
    return bubble_index;
}

//********************************************************************
//**************                              ************************
//**************      Von Neumann             ************************
//**************                              ************************
//********************************************************************
//stabiity criteria is (time_step<=space_step'^2 * rbar^2/6*D) from
int Tissue:: von_Neumann()
{
    double crit=pow(space_step,2)/(diffusivity*6.0);
    if(time_step <= crit)
    {//std::cout<<time_step<<" is <= "<<pow(space_step,2)/(diffusivity*6)<<"\n";
        return 1;}
    else {
        std::cout<<time_step<<" is not <= "<<pow(space_step,2)/(diffusivity*6.0)<<"\n";
        return 0;}
}

//********************************************************************
//**************                              ************************
//**************      Initialise bubbles      ************************
//**************                              ************************
//********************************************************************
void Tissue::initialise_bubbles(int z){
    // sets the initial radius, the locations and initial mass
    if(z<0|z>3){std::cout<<"You didnt choose and elasticity method\n";}
    std::cout<<"z is"<<z<<"\n";
    assert(z>0|z<4);
    
    for (int i=0; i<number_of_bubbles; i++){
        set_bub_radius(i, initial_radius);
        
        if (z==3){
            theBubbles[i]->mass_in_bubble(0, time_step); //pressure_inside(0);
            std::cout<<"gent used \n";
        }
        
        // set_bub_mass(i,0);
        std::cout<<"bubble mass initially is "<<get_bub_massN2(i)+get_bub_massO2(i)<<"\n";
    }
    
    
    if (number_of_bubbles==1){
        std::cout<<"only one bubble placing in center m\n";
        //double depthm;
        //std::cin>>depthm;
        //std::cout<<tissue_size[0]<<" "<<depthm<<"\n";
        //assert(depthm<=tissue_size[0]); //you have put the bubble outside the tissue
        double tissue_mid_point[3]={static_cast<double>((node_number[0]-1))/2,
            static_cast<double>((node_number[1]-1))/2,
            static_cast<double>((node_number[2]-1))/2};
        
        set_bub_locs(0,tissue_mid_point[0],tissue_mid_point[1],tissue_mid_point[2]);
        std::cout<<"bubble coordinates are "<<get_bub_locx(0)<<","<<get_bub_locy(0)<<","<<get_bub_locz(0)<<"\n";
    }
    
    else{
        std::cout<<"randomising locations \n";
        for(int i =0; i<number_of_bubbles; i++){
            theBubbles[i]->randomise_bubble(node_number);
            std::cout<<"bubble "<<theBubbles[i]->get_bubble_number()<<" coordinates are "<<get_bub_locx(i)<<","<<get_bub_locy(i)<<","<<get_bub_locz(i)<<"\n";
        }
        Check_overlap();
        
    }
    for (int bb=0; bb<number_of_bubbles; bb++){
        std::cout<<get_bub_concN2(bb)<<" N2 bubble conc initial \n"<<get_bub_concO2(bb)<<" O2 bubble initial conc\n";
    }
}
//********************************************************************
//**************                              ************************
//**************      Check overlap           ************************
//**************                              ************************
//********************************************************************
void Tissue::Check_overlap(){
    double distance;
    int i;
    for(i=0;i<number_of_bubbles;i++){
        for (int j=i+1; j<number_of_bubbles; j++){
            distance=euclidean_dist(get_bub_locx(i),get_bub_locy(i),get_bub_locz(i),get_bub_locx(j),get_bub_locy(j),get_bub_locz(j));
            std::cout<<"distance between bubble "<<i<<" and "<<j<<" is "<<distance<<"\n";
            
            if(distance<((get_bub_radius(i))+(get_bub_radius(j)))){
                std::cout<<"they overlap moving and running again\n";
                theBubbles[i]->randomise_bubble(node_number);
                Check_overlap();
            }
        }
    }
    
    
}

//********************************************************************
//**************                              ************************
//**************      euclidean dist          ************************
//**************                              ************************
//********************************************************************
double Tissue::euclidean_dist(double x1,double y1,double z1,double x2,double y2,double z2){
    
    double distance=sqrt(
                         SQUARED(x1-x2)+
                         
                         SQUARED(y1-y2)+
                         
                         SQUARED(z1-z2));
    
    return distance;}

//********************************************************************
//**************                              ************************
//**************      compute_label_bub_total  ************************
//**************                              ************************
//********************************************************************
//computes the label grid (unit8*) for all tissue nodes.Nodes in bubbles have value bubble number nodes in bub_label. Boundary nodes have value 1 in boundary_label.
void Tissue::compute_label_bub_total(int tt){
    memset(bub_label,0,node_number[0]*node_number[1]*node_number[2]); //sets a block of memory the size of the conc grid to zero
    memset(boundary_label,0,node_number[0]*node_number[1]*node_number[2]);
    for(int bb=0; bb<number_of_bubbles; bb++){
        theBubbles[bb]->compute_bubble_label(bub_label,tt);
        theBubbles[bb]->compute_boundary_label(boundary_label,bub_label,tt);
    }
}


//********************************************************************
//**************                              ************************
//**************      diffusion               ************************
//**************                              ************************
//********************************************************************
void Tissue::diffusion(int tt){
    //std::cout<<unsigned(*(boundary_label+U(8,12,8)))<<"bound label 8,8,8 \n";
    // do all the difficult work
    clock_t timer;
    timer = clock();
    diffusion_private(conc_sN2,conc_s1N2,conc_sO2,conc_s1O2,P_amb[tt]);
    timer=clock()-timer;
    //std::cout<<conc_s1N2[U(0,0,0)]<<" conc at 0 0 0 of tissue\n";
    
    //std::cout<<conc_s1N2[U(10,10,10)]<<" conc at 10 10 10 of tissue\n";
    
    // std::cout<<"time for individual diffusion step is "<<float(timer)/CLOCKS_PER_SEC <<"secs \n";
    //writes the conc_grid every even time step to give a complete set (should commented out when running full simulations)
    //    if (tt % 500 ==0){
    //        std::cout<<"writing conc dat \n";
    //       // std::cout<<conc_s1N2[U(32,32,32)]<<" conc at 32 32 32 of tissue\n";
    //
    //         write_conc_data(tt);
    //    }
    // double surface_area= 2(node_number[0]*node_number[1])+
    //(conc_sO2[U(0,0,0)]-conc_s1O2[U(0,0,0)])*time_step];
    // swap the two conc time points
    double* tempconcO2=conc_sO2; //old time point is moved to a temporary
    conc_sO2=conc_s1O2; //old conc time point is replaced by s+1 time point
    conc_s1O2=tempconcO2; //new time point is set to old this will be overwirtten next time through it would be nice to set it to zero to check but as this is done using pointers that would involve looping over the entire grid again which seems a waste of time??
    
    double* tempconcN2=conc_sN2; //old time point is moved to a temporary
    conc_sN2=conc_s1N2; //old conc time point is replaced by s+1 time point
    conc_s1N2=tempconcN2; //new time point is set to old this will be overwirtten next time through through it would be nice to set it to zero to check but as this is done using pointers that would involve looping over the entire grid again which seems a waste of time??
}

void Tissue::diffusion_private(const double* c_sN2,double* c_s1N2,const double* c_sO2,double* c_s1O2,double P_amb){
    //Finite difference FTCS diffusion for the tissue block uses the bubble points function
    //to implement no mass flux at the bubble boundaries if the bubble is less than min size.
    // Apply the boundary conditions to conc_s1
    
    // if(pp_fraction_N2*P_amb != c_sN2[U(0,0,0)]){
    
    
    
    
    // clock_t timer_total=0;
    // clock_t timer_bubpoint;
    
    //FTCS method on interior (ie. non-boundary)//
    double factor2 = (diffusivity * time_step)/(rbar*rbar); // prefactor in diffusion approximation
    
    for (int i=1; i<node_number[0]-1; i++) for (int j=1; j<node_number[1]-1; j++) for (int k=1; k<node_number[2]-1; k++){
        double delta_x1=space_step/rbar, delta_x2=space_step/rbar, delta_y1=space_step/rbar, delta_y2=space_step/rbar, delta_z1=space_step/rbar, delta_z2=space_step/rbar;
        //if you're in a bubble ignore
        if (bub_label[U(i,j,k)]){
            // std::cout<<i<<j<<k<<" is in bubble"<<unsigned(bub_label[U(i,j,k)])<<"\n";
            c_s1N2[U(i,j,k)]=get_bub_concN2(find_bubble_index(bub_label[U(i,j,k)]));
            c_s1O2[U(i,j,k)]=get_bub_concO2(find_bubble_index(bub_label[U(i,j,k)]));
        }
        // std::cout<<c_s1O2[U(i,j,k)]<<" , "<<c_s1N2[U(i,j,k)]<<" line 384\n";}
        //std::cout<<"conc test "<<c_s1N2[U(32,32,32)]<<"\n";} //this will jump to outside the for loop it is inside (in this case that is the whole diffusion loop)
        
        else{
            double conc_valsN2[7];
            double conc_valsO2[7];
            conc_valsN2[0] = c_sN2[U(i,j,k)];
            conc_valsN2[1] = c_sN2[U(i+1,j,k)];
            conc_valsN2[2] = c_sN2[U(i-1,j,k)];
            conc_valsN2[3] = c_sN2[U(i,j+1,k)];
            conc_valsN2[4] = c_sN2[U(i,j-1,k)];
            conc_valsN2[5] = c_sN2[U(i,j,k+1)];
            conc_valsN2[6] = c_sN2[U(i,j,k-1)];
            
            conc_valsO2[0] = c_sO2[U(i,j,k)];
            conc_valsO2[1] = c_sO2[U(i+1,j,k)];
            conc_valsO2[2] = c_sO2[U(i-1,j,k)];
            conc_valsO2[3] = c_sO2[U(i,j+1,k)];
            conc_valsO2[4] = c_sO2[U(i,j-1,k)];
            conc_valsO2[5] = c_sO2[U(i,j,k+1)];
            conc_valsO2[6] = c_sO2[U(i,j,k-1)];
            
            //std::cout<<conc_valsN2[0]<<","<<i<<j<<k<<" \n "<<conc_valsN2[1]<<" , "<<conc_valsN2[2]<<" , "<<conc_valsN2[3]<<" , "<<conc_valsN2[4]<<" , "<<conc_valsN2[5]<<" , "<<conc_valsN2[6]<<"\n ";
            
            //if a nearest neighbour is in the bubble and the bubble is at min size. The conc. of the nearest neighbour is changed to be the same as the boundary point. This will give a no flux boundary condition.
            //   std::cout<<"stopping point "<<get_bub_radius(boundary_label[U(i,j,k)]-1)<<"\n";
            if(boundary_label[U(i,j,k)] && get_bub_radius(find_bubble_index(boundary_label[U(i,j,k)]))==min_radius){
                // no_flux_conditions
                if(bub_label[U(i+1,j,k)])
                    conc_valsN2[1]=conc_valsN2[0]; conc_valsO2[1]=conc_valsO2[0];
                if(bub_label[U(i-1,j,k)])
                    conc_valsN2[2]=conc_valsN2[0]; conc_valsO2[2]=conc_valsO2[0];
                if(bub_label[U(i,j+1,k)])
                    conc_valsN2[3]=conc_valsN2[0]; conc_valsO2[3]=conc_valsO2[0];
                if(bub_label[U(i,j-1,k)])
                    conc_valsN2[4]=conc_valsN2[0]; conc_valsO2[4]=conc_valsO2[0];
                if(bub_label[U(i,j,k+1)])
                    conc_valsN2[5]=conc_valsN2[0]; conc_valsO2[5]=conc_valsO2[0];
                if(bub_label[U(i,j,k-1)])
                    conc_valsN2[6]=conc_valsN2[0];  conc_valsO2[6]=conc_valsO2[0];
            }
            
            //    std::cout<<"distances are "<<delta_x2<<","<<delta_x1<<","<<delta_y2<<","<<delta_x1<<","<<delta_z2<<","<<delta_z1<<" for point "  <<i<<","<<j<<","<<k<<"\n";
            
            //  finite difference accounting for true distance of points near the bubble boundary.
            if(boundary_label[U(i,j,k)] and bub_label[U(i,j,k)]==0){
                std::vector<double> bound_concs;
                bound_concs= Inverse_dist_interpolate(i, j, k, conc_valsN2, conc_valsO2);
                c_s1N2[U(i,j,k)]=bound_concs[0];
                c_s1O2[U(i,j,k)]=bound_concs[1];
                goto endforloop;
                
            }
            //expressions for finite difference approximations with unequal distance. Simpifies to standard centered FD approx. for delta1=delta2 not delta1=xi-xi-1, delta2= xi-xi+1;
            
            
            double d2cN2_dx2=2*(delta_x2*conc_valsN2[2]-(delta_x1+delta_x2)*conc_valsN2[0]+delta_x1*conc_valsN2[1])/(delta_x1*delta_x2*(delta_x1+delta_x2));
            double   d2cN2_dy2=2*(delta_y2*conc_valsN2[4]-(delta_y1+delta_y2)*conc_valsN2[0]+delta_y1*conc_valsN2[3])/(delta_y1*delta_y2*(delta_y1+delta_y2));
            double  d2cN2_dz2=2*(delta_z2*conc_valsN2[6]-(delta_z1+delta_z2)*conc_valsN2[0]+delta_z1*conc_valsN2[5])/(delta_z1*delta_z2*(delta_z1+delta_z2));
            
            double  d2cO2_dx2=2*(delta_x2*conc_valsO2[2]-(delta_x1+delta_x2)*conc_valsO2[0]+delta_x1*conc_valsO2[1])/(delta_x1*delta_x2*(delta_x1+delta_x2));
            double   d2cO2_dy2=2*(delta_y2*conc_valsO2[4]-(delta_y1+delta_y2)*conc_valsO2[0]+delta_y1*conc_valsO2[3])/(delta_y1*delta_y2*(delta_y1+delta_y2));
            double  d2cO2_dz2=2*(delta_z2*conc_valsO2[6]-(delta_z1+delta_z2)*conc_valsO2[0]+delta_z1*conc_valsO2[5])/(delta_z1*delta_z2*(delta_z1+delta_z2));
            
            
            //This line is the actual equation for conc at the new time point
            c_s1N2[U(i,j,k)]=factor2*(d2cN2_dx2 + d2cN2_dy2 + d2cN2_dz2)+conc_valsN2[0];
            c_s1O2[U(i,j,k)]=factor2*(d2cO2_dx2 + d2cO2_dy2 + d2cO2_dz2)+conc_valsO2[0];
        } // end of the else
        // std::cout<<"co2 "<<c_s1O2[U(i,j,k)]<<" cN2 "<<c_s1N2[U(i,j,k)]<<"\n";
    endforloop:;
    }//end for i-j-k
    
    
    
    // update the boundary conditions and caluclate flux throught the tissue outer boundary, flux is caluclated as (outer - inner)/dr so outward flux from tissue to surround is -ve
    double tissueflux=0;
    // FOR NO FLUX AT Left FACE UNCOMMENT
    //    for (int j=0; j<node_number[1]; j++) for (int k =0; k<node_number[2]; k++){      //Left face no flux
    //        c_s1N2[U(0,j,k)] = c_s1N2[U(1,j,k)];
    //        c_s1O2[U(0,j,k)] = c_s1O2[U(1,j,k)];}
    
    //
    for (int j=0; j<node_number[1]; j++) for (int k =0; k<node_number[2]; k++){//left face henry's Law
        tissueflux += factor2*((c_sN2[U(0,j,k)]-c_sN2[U(1,j,k)])+(c_sO2[U(0,j,k)]-c_sO2[U(1,j,k)])); //tissue flux is dc/dr *D*area*timestep (area is 1 dr is also 1)
        c_s1N2[U(0,j,k)] = pp_fraction_N2*P_amb;
        c_s1O2[U(0,j,k)] =(k_h_O2/k_h_N2)*pp_fraction_O2*P_amb;
        tissueflux +=((c_s1N2[U(0,j,k)]-c_s1N2[U(0,j,k)])+(c_s1O2[U(0,j,k)]-c_s1O2[U(0,j,k)]));} //for boundary points
    
    
    // FOR NO FLUX AT right FACE UNCOMMENT
    //    for (int j=0; j<node_number[1]; j++) for (int k =0; k<node_number[2]; k++){      //right face no flux
    //        c_s1N2[U(node_number[0]-1,j,k)] = c_s1N2[U(node_number[0]-2,j,k)];
    //        c_s1O2[U(node_number[0]-1,j,k)] = c_s1O2[U(node_number[0]-2,j,k)];}
    //
    for (int j=0; j<node_number[1]; j++) for (int k =0; k<node_number[2]; k++){
        //right face Henry's law
        tissueflux += factor2*((c_sN2[U(node_number[0]-1,j,k)]-c_sN2[U(node_number[0]-2,j,k)])+(c_sO2[U(node_number[0]-1,j,k)]-c_sO2[U(node_number[0]-2,j,k)])); //tissue flux is dc/dr *D*area*timestep (area is 1 dr is also 1)
        c_s1N2[U(node_number[0]-1,j,k)] = pp_fraction_N2*P_amb;
        c_s1O2[U(node_number[0]-1,j,k)]=(k_h_O2/k_h_N2)*pp_fraction_O2* P_amb;
        tissueflux += (c_s1N2[U(node_number[0]-1,j,k)]-c_sN2[U(node_number[0]-1,j,k)])+(c_s1O2[U(node_number[0]-1,j,k)]-c_sO2[U(node_number[0]-1,j,k)]);}
    
    // FOR NO FLUX AT TOP FACE UNCOMMENT
    //    for (int i=0; i<node_number[0]; i++) for (int k =0; k<node_number[2]; k++){      //top face no flux
    //        c_s1N2[U(i,0,k)] = c_s1N2[U(i,1,k)];
    //        c_s1O2[U(i,0,k)] = c_s1O2[U(i,1,k)];}
    
    //
    for (int i=0; i<node_number[0]; i++) for (int k =0; k<node_number[2]; k++){
        //top face Henry's law
        tissueflux += factor2*((c_sN2[U(i,0,k)]-c_sN2[U(i,1,k)])+(c_sO2[U(i,0,k)]-c_sO2[U(i,1,k)]));//tissue flux is dc/dr *D*area*timestep (area is 1 dr is also 1);
        c_s1N2[U(i,0,k)] = pp_fraction_N2*P_amb;
        c_s1O2[U(i,0,k)]=(k_h_O2/k_h_N2)*pp_fraction_O2* P_amb;
        tissueflux += (c_s1N2[U(i,0,k)]-c_sN2[U(i,0,k)])+(c_s1O2[U(i,0,k)]-c_sO2[U(i,0,k)]);}
    
    
    // FOR NO FLUX AT BOTTOM FACE UNCOMMENT
    for (int i=0; i<node_number[0]; i++) for (int k =0; k<node_number[2]; k++){      //bottom face no flux
        c_s1N2[U(i,node_number[1]-1,k)] = c_s1N2[U(i,node_number[1]-2,k)];
        c_s1O2[U(i,node_number[1]-1,k)] = c_s1O2[U(i,node_number[1]-2,k)];}
    
    //    for (int i=0; i<node_number[0]; i++) for (int k =0; k<node_number[2]; k++){      //bottom face henry's law
    //  tissueflux += factor2*((c_sN2[U(i,node_number[1]-1,k)]-c_sN2[U(i,node_number[1]-2,k)])+(c_sO2[U(i,node_number[1]-1,k)]-c_sO2[U(i,node_number[1]-2,k)])); //tissue flux is dc/dr *D*area*timestep (area is 1 dr is also 1)
    //        c_s1N2[U(i,node_number[1]-1,k)] = pp_fraction_N2*P_amb;
    //        c_s1O2[U(i,node_number[1]-1,k)] = (k_h_O2/k_h_N2)*pp_fraction_O2* P_amb;
    //        tissueflux += (c_s1N2[U(i,node_number[1]-1,k)]-c_sN2[U(i,node_number[1]-1,k)])+(c_s1O2[U(i,node_number[1]-1,k)]-                    c_sO2[U(i,node_number[1]-1,k)])); }
    //
    
    //FOR NO FLUX AT front FACE UNCOMMENT
    //    for (int i=0; i<node_number[0]; i++) for (int j =0; j<node_number[2]; j++){      //front face no flux
    //        c_s1N2[U(i,j,0)] = c_s1N2[U(i,j,1)];
    //        c_s1O2[U(i,j,0)] = c_s1O2[U(i,j,1)];}
    
    for (int i=0; i<node_number[0]; i++) for (int j =0; j<node_number[1]; j++){       //front face henry's law
        tissueflux += factor2*((c_sN2[U(i,j,0)]-c_sN2[U(i,j,1)])+(c_sO2[U(i,j,0)]-c_sO2[U(i,j,1)]));//tissue flux is dc/dr *D*area*timestep (area is 1 dr is also 1);
        c_s1N2[U(i,j,0)] = pp_fraction_N2*P_amb;
        c_s1O2[U(i,j,0)] = (k_h_O2/k_h_N2)*pp_fraction_O2*P_amb;
        tissueflux += (c_s1N2[U(i,j,0)]-c_sN2[U(i,j,0)])+(c_s1O2[U(i,j,0)]-c_sO2[U(i,j,0)]);}
    
    
    // FOR NO FLUX AT BACK FACE UNCOMMENT
    //    for (int i=0; i<node_number[0]; i++) for (int j =0; j<node_number[1]; j++){      //back face no flux
    //        c_s1N2[U(i,j,node_number[2]-1)] = c_s1N2[U(i,j,node_number[2]-2)];
    //        c_s1O2[U(i,j,node_number[2]-1)] = c_s1O2[U(i,j,node_number[2]-2)];}
    
    for (int i=0; i<node_number[0]; i++) for (int j =0; j<node_number[1]; j++){       //back face henry's law
        tissueflux += factor2*((c_sN2[U(i,j,node_number[2]-1)]-c_sN2[U(i,j,node_number[2]-2)])+(c_sO2[U(i,j,node_number[2]-1)]-c_sO2[U(i,j,node_number[2]-2)])); //tissue flux is dc/dr *D*area*timestep (area is 1 dr is also 1)
        c_s1N2[U(i,j,node_number[2]-1)] = pp_fraction_N2*P_amb;
        c_s1O2[U(i,j,node_number[2]-1)] = (k_h_O2/k_h_N2)*pp_fraction_O2*P_amb;
        tissueflux += (c_s1N2[U(i,j,node_number[2]-1)]-c_sN2[U(i,j,node_number[2]-1)])+(c_s1O2[U(i,j,node_number[2]-1)]-c_sO2[U(i,j,node_number[2]-1)]);}
    
    
    
    set_tissue_flux(tissueflux);
    // std::cout<<get_tissue_flux()<<" tissue flux \n";
    // std::cout<<c_s1N2[U(32,32,32)]<<" conc at in bub of tissue\n";
    
}//end function
//********************************************************************
//**************                              ************************
//*************  interpolate boundary points  ************************
//**************                              ************************
//********************************************************************
std::vector<double> Tissue::Inverse_dist_interpolate(int i, int j, int k,double* conc_valsN2,double* conc_valsO2){
    
    double delta_x1=space_step/rbar, delta_x2=space_step/rbar, delta_y1=space_step/rbar, delta_y2=space_step/rbar, delta_z1=space_step/rbar, delta_z2=space_step/rbar;
    int indx=find_bubble_index(boundary_label[U(i,j,k)]);
    double distsx=theBubbles[indx]->xdist(i,j,k);
    double distsy=theBubbles[indx]->ydist(i,j,k);
    double distsz=theBubbles[indx]->zdist(i,j,k);
    
    // std::cout<<distsx<<","<<distsy<<","<<distsz<<" dists \n";
    
    if(distsx<0.0 && distsx>-1.0) delta_x2=fabs(distsx); // it is negative x direction to the bubble so delta2 should change
    if(distsx>0.0 && distsx<1.0) delta_x1=fabs(distsx);
    if(distsy<0.0 && distsy>-1.0) delta_y2=fabs(distsy); // it is negative x direction to the bubble so delta2 should change
    if(distsy>0.0 && distsy<1.0) delta_y1=fabs(distsy);
    if(distsz<0.0 && distsz>-1.0) delta_z2=fabs(distsz); // it is negative x direction to the bubble so delta2 should change
    if(distsz>0.0 && distsz<1.0) delta_z1=fabs(distsz);
    // std::cout<<"distances are "<<delta_x1<<","<<delta_x2<<","<<delta_y1<<","<<delta_y2<<","<<delta_z1<<","<<delta_z2<<" for point "  <<i<<","<<j<<","<<k<<"\n";
    double weightsum=1/pow(delta_x1,2) + 1/pow(delta_x2,2) + 1/pow(delta_y1,2) + 1/pow(delta_y2,2) + 1/pow(delta_z1,2)+ 1/pow(delta_z2,2);
    double weights[6]={1.0/pow(delta_x1,2), 1/pow(delta_x2,2), 1/pow(delta_y1,2), 1/pow(delta_y2,2), 1/pow(delta_z1,2), 1/pow(delta_z2,2)};
    double valueN2=0;
    double valueO2=0;
    for(int p=0;p<6;p++){
        valueN2 += weights[p]*conc_valsN2[p+1];
        valueO2 += weights[p]*conc_valsO2[p+1];
        //std::cout<<weights[p]<<",";
    }
    // std::cout<<"\n"<<valueN2<<" N2 \n"<<valueO2<<" O2 "<<"\n"<<weightsum<<"sumweight \n";
    
    double c_s1N2=(valueN2)/weightsum;
    double c_s1O2=(valueO2)/weightsum;
    
    std::vector<double> concsN2O2;
    concsN2O2.push_back(c_s1N2); concsN2O2.push_back(c_s1O2);
    
    return concsN2O2;
}

//********************************************************************
//**************                              ************************
//**************   write conc data            ************************
//**************                              ************************
//********************************************************************
//writes the full conc_grid mulitplying by cbar to get real values make sure this isnt running if you are doing a full simulation.
void Tissue::write_conc_data(int t){
    
    //
    //    /// **** open a file for writing the data to ******
    //    //   ---------------------------------------------
    std::ofstream outfile;
    outfile.open("/Users/clairewalsh/Dropbox/Modelling/dan_May_2014/dan_May_2014/outputs/mass_conservation/test1conc.txt", std::ios::app);
    if(!outfile.is_open()){std::cout<<"error opening conc_writing_file.txt /n";}
    std::cout << "Writing to the file" << std::endl;
    //
    //    // ****** write conc data into the file. ********
    //
    outfile<<"N2 grid"<<std::endl;
    for(int i=0; i<node_number[0]; i++){
        for(int j=0; j<node_number[1];j++){
            int k=node_number[2]/2;
            outfile << conc_sN2[U(i,j,k)]*cbar<<",";
        }
        outfile<<std::endl;
    }
    // outfile<<"O2 grid"<<std::endl;
    //for(int i=0; i<node_number[0]; i++){
    //  for(int j=0; j<node_number[1];j++){
    //    int k=node_number[2]/2;
    //  outfile << conc_sO2[U(i,j,k)]*cbar<<"  "<<i<<"  "<<j<<"  "<<k<< std::endl;
    
    //}
    //}
    
    
    // ****** close the opened file. **********
    outfile.close();
}




//********************************************************************
//**************                              ************************
//**************     alternative dc/dr        ************************
//**************                              ************************
//********************************************************************
//applied shperical dc/dr method on each bubble in the tissue grid
void Tissue::bub_dc_dr_spherical(int tt){
    
    for(int bb =0; bb<number_of_bubbles; bb++){
        if(P_amb[tt]>=P_amb[tt-1]& get_bub_radius(bb)==min_radius){
            //  std::cout<<"imposing no flux dcdr==0\n";
            theBubbles[bb]->set_dcdrN2(0.0);
            theBubbles[bb]->set_dcdrO2(0.0);
        }
        else{
            //double diffusion_dist=5e-5/rbar; //it should be whatever dr is used in spherical transform//theBubbles[bb]->get_radius()/diff_region_prop; //sets the diffusion region distance.
            double angle = (PI/4.0);
            std::vector<std::vector<double> > pts_r1=theBubbles[bb]->spherical_tansform_r1(diffusion_dist,angle); //performs co-ordinate transformation
            std::vector<std::vector<double> > pts_r2=theBubbles[bb]->spherical_tansform_r1(2*diffusion_dist,angle);
            std::vector<double> newconcvecr1N2;
            std::vector<double> newconcvecr1O2;
            std::vector<double> newconcvecr2N2;
            std::vector<double> newconcvecr2O2;
            assert(pts_r1.size()==pts_r2.size());
            for(int i=0; i<pts_r1.size(); i++){
                std::vector<double> newconcr1=interpolate(pts_r1[i],bb); // sets the concentrations of these new points based on
                std::vector<double> newconcr2=interpolate(pts_r2[i],bb); //interpolation from the nearest neighbours
                
                newconcvecr1N2.push_back(newconcr1[0]);
                newconcvecr1O2.push_back(newconcr1[1]);//new points with their concentrations
                newconcvecr2N2.push_back(newconcr2[0]);
                newconcvecr2O2.push_back(newconcr2[1]);//new points with their concentrations
            }
            
            theBubbles[bb]->dcdr_sphere(newconcvecr1N2,newconcvecr1O2,newconcvecr2N2,newconcvecr2O2,1.0,tt);  //calculates dc/dr based on these new points and the internal bubble concentration
            pts_r1.clear();
            newconcvecr1N2.clear();
            newconcvecr1O2.clear();
            newconcvecr2N2.clear();
            newconcvecr2O2.clear();
        }
        
    }
}

//********************************************************************
//**************                              ************************
//**************     interpolate dc/dr        ************************
//**************                              ************************
//********************************************************************
//this function finds the concentration of a point given in spherical co-ordinates by interpolation of the surrounding points. It takes input from function spherical transform which is in the bubble class.
std::vector<double> Tissue::interpolate(std::vector<double> pts,int bb){
    double newconc_N2;
    double newconc_O2;
    
    if(floor(pts[0])==pts[0]&&floor(pts[1])==pts[1]&&floor(pts[2])==pts[2]){
        newconc_N2=conc_sN2[U(int (pts[0]),int (pts[1]),int (pts[2]))];
        newconc_O2=conc_sO2[U(int (pts[0]),int (pts[1]),int (pts[2]))];}//for the points that already lie on a grid point no need to find their conc by interpolation
    else{
        std::vector<std::vector< double > > surrounding_pts=nearest_neighbours(pts,get_bub_locx(bb),get_bub_locy(bb), get_bub_locz(bb), get_bub_radius(bb)); // looks at as many of the 8 nearest neighbours that arent in the bubble
        std::vector<double> weights;
        std::vector<double> valueN2;
        std::vector<double> valueO2;
        assert(surrounding_pts.size()<=8);
        for (int i=0; i<surrounding_pts.size(); i++){
            //weighted averager using shepards Method with exponent power 2
            weights.push_back(1/pow(euclidean_dist(surrounding_pts[i][0],surrounding_pts[i][1],surrounding_pts[i][2],pts[0],pts[1],pts[2]),2));
            valueN2.push_back(conc_sN2[U(int(surrounding_pts[i][0]),int(surrounding_pts[i][1]),int(surrounding_pts[i][2]))]*weights[i]);
            valueO2.push_back(conc_sO2[U(int(surrounding_pts[i][0]),int(surrounding_pts[i][1]),int(surrounding_pts[i][2]))]*weights[i]);
            //sum of the values where value is conc*weight)
            
        }
        
        while (surrounding_pts.size()<8){ //if this is the point nearest the bubble and it does not have 8 nearest neighbours another point on the bubble edge is used in the interpolation it will have bubble conc. i.e Pb*kh and will be a distance of 1 from the r1 point. (this should only happen on the R1 point)
            weights.push_back(1.0);
            valueN2.push_back(get_bub_concN2(bb));
            valueO2.push_back(get_bub_concO2(bb));
            std::vector<double>bubble_pt={0,0,0};
            surrounding_pts.push_back(bubble_pt);
        }
        
        double  weight_tot= std::accumulate(weights.begin(), weights.end(), 0.0);
        double valueN2_tot=std::accumulate(valueN2.begin(), valueN2.end(), 0.0);
        double valueO2_tot=std::accumulate(valueO2.begin(), valueO2.end(), 0.0);
        
        newconc_N2=valueN2_tot/weight_tot;
        newconc_O2=valueO2_tot/weight_tot; // conc is two sums divided by each other
        //
        // std::cout<<"("<<pts[0]<<","<<pts[1]<<","<<pts[2]<<"),  "<<valueN2_tot<<" "<<weight_tot<<"\n"<<valueO2_tot<<"\n";
    }
    
    
    std::vector<double> newconc_N2_O2;
    newconc_N2_O2.insert( newconc_N2_O2.end(),  {newconc_N2,newconc_O2} );
    return newconc_N2_O2;
}


//********************************************************************
//**************                              ************************
//**************    nearest neighbour         ************************
//**************                              ************************
//********************************************************************

std::vector<std::vector<double > > Tissue::nearest_neighbours(std::vector<double> pt,double locx,double locy,double locz, double radius){
    
    std::vector<std::vector<double > > surrounding_pts;
    std::vector<double> row;
    // std::cout<<pt[0]<<","<<pt[1]<<","<<pt[2]<<" this is the point going in \n";
    
    int lowerx;
    int lowery;
    int lowerz;
    int higherx;
    int highery;
    int higherz;
    
    
    if(floor(pt[0])==pt[0]){lowerx=pt[0]-1;}  //if one of the coords is exactly on the grid point then we arbitrarily take the lower 8 nearest neighbours
    else{lowerx=floor(pt[0]);}
    if(floor(pt[1])==pt[1]){lowery=pt[1]-1;}
    else{lowery=floor(pt[1]);}
    if(floor(pt[2])==pt[2]){lowerz=pt[2]-1;}
    else{lowerz=floor(pt[2]);}
    if(ceil(pt[0])==pt[0]){higherx=pt[0];}
    else{higherx=ceil(pt[0]);}
    if(ceil(pt[1])==pt[1]){highery=pt[1];}
    else{highery=ceil(pt[1]);}
    if(ceil(pt[2])==pt[2]){higherz=pt[2];}
    else{higherz=ceil(pt[2]);}
    
    
    //std::cout<<lowerx<<" "<<higherx<<" lower and higher x \n";
    //std::cout<<lowery<<" "<<highery<<" lower and higher y \n";
    //std::cout<<lowerz<<" "<<higherz<<" lower and higher z \n";
    
    
    
    //all the nearest neighbours have been found above. They are only returned in the surrounding points vector for interpolation if they are not in the bubble
    
    row.push_back(lowerx);
    row.push_back(lowery);
    row.push_back(lowerz);
    // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<"surrounding pts\n";
    if(euclidean_dist(row[0],row[1],row[2],locx,locy,locz)>=radius){
        surrounding_pts.push_back(row);
        // std::cout<<surrounding_pts.size()<<" at line 150 \n";
    }
    row.clear();
    
    row.push_back(higherx);
    row.push_back(highery);
    row.push_back(higherz);
    // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<"\n";
    if(euclidean_dist(row[0],row[1],row[2],locx,locy,locz)>=radius){
        surrounding_pts.push_back(row);
        // std::cout<<surrounding_pts.size()<<" at line 161 \n";
    }
    row.clear();
    
    
    row.push_back(higherx);
    row.push_back(lowery);
    row.push_back(higherz);
    // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<"\n";
    if(euclidean_dist(row[0],row[1],row[2],locx,locy,locz)>=radius){
        surrounding_pts.push_back(row);
        //std::cout<<surrounding_pts.size()<<" at line 171 \n";
    }
    
    row.clear();
    
    row.push_back(higherx);
    row.push_back(lowery);
    row.push_back(lowerz);
    // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<"\n";
    if(euclidean_dist(row[0],row[1],row[2],locx,locy,locz)>=radius){
        surrounding_pts.push_back(row);
        // std::cout<<surrounding_pts.size()<<" at line 183 \n";
    }
    
    row.clear();
    
    row.push_back(higherx);
    row.push_back(highery);
    row.push_back(lowerz);
    // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<"\n";
    if(euclidean_dist(row[0],row[1],row[2],locx,locy,locz)>=radius){
        surrounding_pts.push_back(row);
        // std::cout<<surrounding_pts.size()<<" at line 194 \n";
    }
    
    row.clear();
    
    row.push_back(lowerx);
    row.push_back(highery);
    row.push_back(higherz);
    //  std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<"\n";
    if(euclidean_dist(row[0],row[1],row[2],locx,locy,locz)>=radius){
        surrounding_pts.push_back(row);
        //     std::cout<<surrounding_pts.size()<<" at line 204 \n";
    }
    
    row.clear();
    
    row.push_back(lowerx);
    row.push_back(lowery);
    row.push_back(higherz);
    // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<"\n";
    if(euclidean_dist(row[0],row[1],row[2],locx,locy,locz)>=radius){
        surrounding_pts.push_back(row);
        //  std::cout<<surrounding_pts.size()<<" at line 215 \n";
        
    }
    
    row.clear();
    
    row.push_back(lowerx);
    row.push_back(highery);
    row.push_back(lowerz);
    //std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<"\n";
    if(euclidean_dist(row[0],row[1],row[2],locx,locy,locz)>=radius){
        surrounding_pts.push_back(row);
        //    std::cout<<surrounding_pts.size()<<" at line 225 \n";
    }
    row.clear();
    
    // std::cout<<surrounding_pts.size()<<"\n";
    //
    //        for(int i=0; i<surrounding_pts.size(); i++){
    //            std::cout<<surrounding_pts[i][0]<<" "<<surrounding_pts[i][1]<<" "<<surrounding_pts[i][2]<<" these are the surrounding points not in the bubble  "<<"\n";
    //        }
    
    // if (surrounding_pts.size()< 8){std::cout<<"next to bubble";}
    return surrounding_pts;
    
}

//********************************************************************
//**************                              ************************
//**************          bub_dr_dt           ************************
//**************                              ************************
//********************************************************************
//applies dr/dt on each bubble in the tissue grid gent elasticity
void Tissue:: bub_drdt(int tt){
    for (int i=0; i<number_of_bubbles; i++){
        
        theBubbles[i]->drdt(tt, get_time_step());}
    // std::cout<<theBubbles[0]->get_radius()<<" radius new \n";
    // compute_label_bub_total(tt);
    
    
    
    
}
//********************************************************************
//**************                              ************************
//**************          bub_dr_dt_gern       ************************
//**************                              ************************
//********************************************************************
//applies dr/dt on each bubble in the tissue grid (Gernhard elasticity)
void Tissue:: bub_drdt_gern(int tt){
    
    for (int i=0; i<number_of_bubbles; i++){
        
        theBubbles[i]->drdt_gern(tt, get_time_step());
        
    }
    //compute_label_bub_total(tt);
}

//********************************************************************
//**************                              ************************
//**************          bub_dr_dt_lap       ************************
//**************                              ************************
//********************************************************************
//applies dr/dt on each bubble in the tissue grid (Gernhard elasticity)
void Tissue:: bub_drdt_lap(int tt){
    for (int i=0; i<number_of_bubbles; i++){
        theBubbles[i]->drdt_laplace(tt, get_time_step());
    }
    // compute_label_bub_total(tt);
    
    
    
}
//********************************************************************
//**************                              ************************
//**************       check radius folder    ************************
//**************                              ************************
//********************************************************************
void Tissue:: check_radius_folder(std::string filepathandname){
    std::ofstream outfile;
    outfile.open(filepathandname, std::ios::app);
    if(!outfile.is_open()){std::cout<<"error opening radius output file /n"; assert(1==0);}
    outfile.close();
    
}
//********************************************************************
//**************                              ************************
//**************       write radius data       ************************
//**************                              ************************
//********************************************************************
//The data is written out into a txt file with columns: time step, bubble number, radius, bubble concentration. (as you see it in the code output without the external pressure).
void Tissue::write_radius_data(std::string dive_name, std::string filepathandname){
    std::cout << "Writing radii to the file" << std::endl;
    // **** open a file for writing the data to ******
    //   ---------------------------------------------
    std::ofstream outfile;
    outfile.open(filepathandname, std::ios::app);
    if(!outfile.is_open()){std::cout<<"error opening radius output file /n";}
    std::cout << "Writing radii to the file" << std::endl;
    
    // ****** write radius data into the file. ********
    outfile<<"using ro=r_initial \n Dive profile is "<<dive_name<<"\n"<<"D= "<<diffusivity<<"\n"<<"dt= "<<get_time_step()<<"\n dx= "<<get_space_step()<<"\n ro= "<<initial_radius<<"\n sigma (N/m)= "<<sigma*(100000.0*rbar*pbar)<<"\n k_h_O2= "<<k_h_O2<<"\n k_h_N2= "<<k_h_N2<<"\n RT_N2"<<RT_N2<<"\n RT_O2"<<RT_O2<<"\n tissue size "<<get_tissue_size_x()<<","<<get_tissue_size_y()<<","<<get_tissue_size_z()<<"\n mu (bar) = "<<mu*pbar<<"\n dc"<<diffusion_dist<<std::endl;
    for(int bb=0; bb<number_of_bubbles; bb++){
        outfile<<"\n bubble location"<<get_bub_locx(bb)<<","<<get_bub_locy(bb)<<","<<get_bub_locz(bb)<<std::endl<<std::endl<<
        "time_step"<<"\t"<<"Bubble number"<<"\t"<<"Radius"<<"\t"<<"\t"<<"bubble mass"<<"\t"<<"tissue mass"<<"\t"<<"tissue flux"<<std::endl;
        for(int tt=0; tt<time_of_dive/get_time_step() /*((time_of_dive*60)/get_time_step())*/;tt++){
            outfile <<tt*get_time_step()<<"\t"<<theBubbles[bb]->get_bubble_number()<<"\t"<<theBubbles[bb]->get_radius(tt)*rbar<<"\t"<<total_bubble_mass[tt]<<"\t"<<total_mass_tissue[tt]<<"\t"<<tissue_flux[tt]<<std::endl;
        }
    }
    
    // ****** close the opened file. **********
    outfile.close();
}



////*************************************************************************************************************
////***************************              param check                                    *********************
////*************************************************************************************************************
void Tissue::param_check(int tt){
    //std::cout<<"time step = "<<tt<<"\n";
    for(int bb=0; bb<number_of_bubbles; bb++){
        std::cout<<"time      "<<"bubble "<<"   radius   "<<"  bubbble pressure   "<<"  external pressure       "<<"     bubble_mass            "<<"       tissue mass              "<<"\n"<<tt<<"          "<<theBubbles[bb]->get_bubble_number()<<"         "<<theBubbles[bb]->get_radius()<<"      "<<theBubbles[bb]->get_concN2()+theBubbles[bb]->get_concO2()<< "                  "<<P_amb[tt]<<"                  "<<theBubbles[bb]->get_mmassN2()+theBubbles[bb]->get_mmassO2()<<"              "<<total_mass_tissue[tt]<<"\n";
    }
}




//********************************************************************
//**************                              ************************
//**************    bubble_check              ************************
//**************                              ************************
//********************************************************************
void Tissue::bubble_check(std::string filepath,int t,int z){
    for(int i=0;i<number_of_bubbles;i++){
        for (int j=i+1; j<number_of_bubbles; j++){
            double distance=euclidean_dist(get_bub_locx(i),get_bub_locy(i),get_bub_locz(i),get_bub_locx(j),get_bub_locy(j),get_bub_locz(j));
            
            if(distance<(get_bub_radius(i)+get_bub_radius(j)+1)){
                std::cout<<"bubbles overlap writing current rad datat\n";
                coalescence(filepath, i, j, t, z);
                // exit(EXIT_SUCCESS);
            }
        }
    }
    edge_check();
    
}
//********************************************************************
//**************                              ************************
//**************    edge_check                ************************
//**************                              ************************
//********************************************************************
//checks the bubbles isn't too near an edge and shift the center if it is.
void Tissue::edge_check(){
    int edge[3]={0, 0, 0};
    for (int b=0; b<number_of_bubbles; b++){
        
        if (ceil(get_bub_radius(b))+get_bub_locx(b)+2>=node_number[0]){edge[0]=1;}
        else if (get_bub_locx(b)-ceil(get_bub_radius(b))-2<=0){edge[0]=2;}
        
        if (ceil(get_bub_radius(b))+get_bub_locy(b)+2>=node_number[1]){edge[1]=1;}
        else if (get_bub_locy(b)-ceil(get_bub_radius(b))-2<=0){edge[1]=2;}
        
        if (ceil(get_bub_radius(b))+get_bub_locz(b)+2>=node_number[2]){edge[2]=1;}
        else if (get_bub_locz(b)-ceil(get_bub_radius(b))-2<=0){edge[2]=2;}
        //else edge=0;
        
        // std::cout<<"edge is "<<edge<<" "<<get_bub_locx(b)-1-get_bub_radius(b)<<" "<<get_bub_locx(b)<<" "<<-get_bub_radius(b)<<"\n";
        if(edge[0]!=0 || edge[1]!=0 || edge[2]!=0){
            
            std::cout<<"Bubble "<<theBubbles[b]->get_bubble_number()<<" is too near the edge and needs to be moved edge array is "<<edge[0]<<" "<<edge[1]<<" "<<edge[2]<<"\n";
            move_bubble(edge,b);
            //edge_check();
            edge[0]=0; edge[1]=0; edge[2]=0;}
    }
    
    
}

//********************************************************************
//**************                              ************************
//**************   move bubble                ************************
//**************                              ************************
//********************************************************************
//checks the bubbles isn't too near an edge and shift the center if it is.
void Tissue::move_bubble(int *edge, int b){
    //moving in the x direction
    if(edge[0]==1){set_bub_locs(b, get_bub_locx(b)-1,get_bub_locy(b),get_bub_locz(b));}
    else if(edge[0]==2){set_bub_locs(b, get_bub_locx(b)+1,get_bub_locy(b),get_bub_locz(b));}
    //moving in the y direction
    if(edge[1] ==1){set_bub_locs(b, get_bub_locx(b),get_bub_locy(b)-1,get_bub_locz(b));}
    else if(edge[1] ==2){set_bub_locs(b, get_bub_locx(b),get_bub_locy(b)+1,get_bub_locz(b));}
    //moving in the z direction
    if(edge[2] ==1){set_bub_locs(b, get_bub_locx(b),get_bub_locy(b),get_bub_locz(b)-1);}
    else if(edge[2] ==2){set_bub_locs(b, get_bub_locx(b),get_bub_locy(b),get_bub_locz(b)+1);}
    
    //else{std::cout<<"EUSTON WE HAVE A PROBLEM!!!!!!!!!!!!!!!!!! (With move bubble)\n";}
    
    std::cout<<"new bubble location is "<<get_bub_locx(b)<<" "<<get_bub_locy(b)<<" "<<get_bub_locz(b)<<"\n";
}

//********************************************************************
//**************                              ************************
//**************   coalescence                ************************
//**************                              ************************
//********************************************************************
//new bubble with location of one of the old ones
//give this a new radius
//get rid of one of the two old ones
//re-allocate bubble numbers
//re-allocate conc. points inside and outside
//
void Tissue::coalescence(std::string filepathname,int b1,int b2,int t,int z){
    std::cout<<"\n\nCoalescence of bubbles "<<theBubbles[b1]->get_bubble_number()<<" and "<<theBubbles[b2]->get_bubble_number()<<"\n\n";
    double radius1cubed=pow(get_bub_radius(b1),3);
    double radius2cubed=pow(get_bub_radius(b2),3);
    double N2mass=get_bub_massN2(b1,t)+get_bub_massN2(b2,t);
    double O2mass=get_bub_massO2(b1,t)+get_bub_massO2(b2,t);
    double new_rad=pow(pow(get_bub_radius(b1),3.0)+pow(get_bub_radius(b2),3.0),1/3.0);
    double N2conc=N2mass*L_N2/((4.0*PI/3.0)*pow(new_rad,3));
    double O2conc=O2mass*L_O2/((4.0*PI/3.0)*pow(new_rad,3));
    
    int smaller_bub;
    int larger_bub;
    if(theBubbles[b1]->get_radius()<theBubbles[b2]->get_radius()){smaller_bub=b1,larger_bub=b2;}
    else{smaller_bub=b2,larger_bub=b1;}
    
    // which ever bubble is larger that's the one that becomes the location of the new coalessed bubble.
    
    
    double newLocation[3]={get_bub_locx(larger_bub),get_bub_locy(larger_bub),get_bub_locz(larger_bub)};
    
    
    
    //    write
    //Write out the radii of the two that coalesed up to the current time point.
    std::ofstream outfile;
    outfile.open(filepathname, std::ios::app);
    if(!outfile.is_open()){std::cout<<"error opening radius output file /n";}
    std::cout << "Writing coalesced radii to the file" << std::endl;
    
    
    for(int bb=0; bb<number_of_bubbles; bb++){
        outfile<<"\n bubble location"<<get_bub_locx(bb)<<","<<get_bub_locy(bb)<<","<<get_bub_locz(bb)<<std::endl<<std::endl<<
        "time_step"<<"\t"<<"Bubble number"<<"\t"<<"Radius"<<"\t"<<"bubble pressure"<<"\t"<<"bubble mass"<<"\t"<<"tissue mass"<<"\t"<<"tissue flux"<<std::endl;
        for(int tt=0; tt<time_of_dive/get_time_step()/*((time_of_dive*60)/get_time_step())*/;tt++){
            outfile <<tt*get_time_step()<<"\t"<<theBubbles[bb]->get_bubble_number()<<"\t"<<theBubbles[bb]->get_radius(tt)*rbar<<"\t"<<"\t"<<total_bubble_mass[tt]<<"\t"<<total_mass_tissue[tt]<<"\t"<<tissue_flux[tt]<<std::endl;
        }
    }
    
    
    
    // ****** write radius data into the file. ********
    outfile<<"D= "<<diffusivity<<"\n"<<"dt= "<<get_time_step()<<"\n dx= "<<get_space_step()<<"\n ro= "<<initial_radius<<"\n sigma= "<<sigma<<"\n k_h_O2= "<<k_h_O2<<"\n k_h_N2= "<<k_h_N2<<"\n RT_N2= "<<RT_N2<<"\n RT_O2= "<<RT_O2<<"\n tissue size "<<get_tissue_size_x()<<","<<get_tissue_size_y()<<","<<get_tissue_size_z()<<"\n mu = "<<mu<<"\n dc="<<diffusion_dist<<std::endl;
    {
        outfile<<"\n bubble location for bubble "<<theBubbles[b1]->get_bubble_number()<<" "<<get_bub_locx(b1)<<","<<get_bub_locy(b1)<<","<<get_bub_locz(b1)<<std::endl<<
        "\n bubble location for bubble "<<theBubbles[b2]->get_bubble_number()<<" "<<get_bub_locx(b2)<<","<<get_bub_locy(b2)<<","<<get_bub_locz(b2)<<std::endl;
        for (int tt=0; tt<t;tt++){
            outfile <<tt*get_time_step()<<"\t"<<b1<<"\t"<<theBubbles[b1]->get_radius(tt)*rbar<<"\t"<<total_bubble_mass[tt]<<"\t"<<total_mass_tissue[tt]<<"\t"<<tissue_flux[tt]<<std::endl;
        }
        
        for (int tt=0; tt<t;tt++){
            outfile   <<tt*get_time_step()<<"\t"<<b2<<"\t"<<theBubbles[b2]->get_radius(tt)*rbar<<"\t"<<total_bubble_mass[tt]<<"\t"<<total_mass_tissue[tt]<<"\t"<<tissue_flux[tt]<<std::endl;
        }
    }
    
    
    // ****** close the opened file. **********
    outfile.close();
    
    theBubbles.push_back(new Bubble (number_of_bubbles+1 , newLocation[0] , newLocation[1] ,newLocation[2],new_rad,t,N2mass,O2mass,N2conc,O2conc));
    //create a new bubble all the vectors in this new bubble are made same length as current ones with zeros for previous time points.  The new bubble has bubble number one greater than the original number of bubbles i.e. there were originally 4 bubbles; 2 and 4 coalesed bubble 5 is created.
    std::cout<<theBubbles[number_of_bubbles]->get_radius()<<" radius back\n";
    std::cout<<theBubbles[number_of_bubbles]->get_radius(t)<<" t radius \n";
    std::cout<<"z is "<<z<<"\n";
    // Finding list of point that were in the two old bubbles and are now not or were not in and now are
    //memset(temp_bub_label,0,node_number[0]*node_number[1]*node_number[2]);
    //    theBubbles[b1]->compute_bubble_label(temp_bub_label);
    //    theBubbles[b2]->compute_bubble_label(temp_bub_label);
    theBubbles[number_of_bubbles]->compute_bubble_label(bub_label,t);
    
    // getting a concentration value for those points which lie outside the new bubble but were inside the old smaller bubble. They are set to the average of the points which surrounded the bubble and are not in the new bubble now on the previous time point.
    std::vector<std::vector<double> > pts=theBubbles[smaller_bub]->spherical_tansform_r1(diffusion_dist,PI/4); //performs co-ordinate transformation
    std::vector<std::vector<double> >surroundingconcvec;
    double tot_conc_N2;
    double tot_conc_O2;
    for(int i=0; i<pts.size(); i++){
        if(euclidean_dist(pts[i][0], pts[i][1], pts[i][2], get_bub_locx(number_of_bubbles), get_bub_locy(number_of_bubbles),get_bub_locz(number_of_bubbles))<=get_bub_radius(number_of_bubbles))
            continue;
        std::vector<double> surroundingconc=interpolate(pts[i],smaller_bub); // sets the concentrations of these new points based on interpolation from the nearest neighbours
        surroundingconcvec.push_back(surroundingconc);   //new points with their concentrations
        tot_conc_N2 += surroundingconc[0];
        tot_conc_O2 += surroundingconc[1];
    }
    
    double avrg_conc_N2 = tot_conc_N2/surroundingconcvec[0].size();
    double avrg_conc_O2 = tot_conc_O2/surroundingconcvec[1].size();
    
    // If point are in the new bubble and weren't previously they are set to the inner bubble conc. those that were in a bubble but are not now are set to the average conc of the surrounding points.
    for (int i=1; i<node_number[0]-1; i++) for (int j=1; j<node_number[1]-1; j++) for (int k=1; k<node_number[2]-1; k++){
        
        if(bub_label[U(i,j,k)]==number_of_bubbles+1){
            conc_sN2[U(i,j,k)]=theBubbles[number_of_bubbles]->get_concN2();
            conc_sO2[U(i,j,k)]=theBubbles[number_of_bubbles]->get_concO2();
        }
        else if(bub_label[U(i,j,k)]==smaller_bub){ //we only need to consider the smaller bubble as the new coalesced bubble is always in the location of the previously larger bubble and will always be slightly bigger so larger bubble doesn't have points that are no longer inside.
            conc_sN2[U(i,j,k)]=avrg_conc_N2;
            conc_sO2[U(i,j,k)]=avrg_conc_O2;
        }
    }
    
    
    
    delete_a_Bubble(b2); //delete the two bubbles that have coalesed
    delete_a_Bubble(b1);
    number_of_bubbles=theBubbles.size(); //reduce the total number of bubbles by 1 (two were removed and one created)
    std::cout<<"after all coalescence the new number of bubbles is now "<<number_of_bubbles<<"\n";
    
    
}

//********************************************************************
//**************                              ************************
//**************   update boundary_pts mass   ************************
//**************                              ************************
//********************************************************************
void Tissue::update_boundary_pts(int tt){
    
    std::vector<double> mass_defecit;
    
    for(int bb=0;bb<number_of_bubbles;bb++){
        if(P_amb[tt]>=P_amb[tt-1]& get_bub_radius(bb)==min_radius){
            //  std::vector<double> flux_solution=mass_bounding_box(tt, bb);
            
            std::vector<double> full_solution=mass_conservation(tt);
        }
        
        else{
            mass_defecit=mass_bounding_box(tt, bb); //call mass of bounding box to calculate the mass difference between half way solution and final solution.
            double mass_added_per_pt=mass_defecit[0]/mass_defecit[2];
            // std::cout<<mass_defecit[2]<<" prop \n";
            //  std::cout<<mass_defecit[0]<<" mass difference between boundaing boxes\n";
            
            //  add the deficit to the boundary points
            for(int i=floor(theBubbles[bb]->get_locations(0))-theBubbles[bb]->get_radius();i<=ceil(theBubbles[bb]->get_locations(0))+theBubbles[bb]->get_radius();i++){
                for(int j=floor(theBubbles[bb]->get_locations(1))-theBubbles[bb]->get_radius();j<=ceil(theBubbles[bb]->get_locations(1))+theBubbles[bb  ]->   get_radius();j++){
                    for(int k=floor(theBubbles[bb]->get_locations(2))-theBubbles[bb]->get_radius();k<=ceil(theBubbles[bb]->get_locations(2))+theBubbles[bb]->get_radius();k++){
                        if(boundary_label[U(i,j,k)]!=0){
                            double prop=(theBubbles[bb]->surf_to_point(i,j,k));
                            //  std::cout<<conc_sN2[U(i,j,k)]<<" conc before,  ";
                            conc_sN2[U(i,j,k)]+=(mass_added_per_pt*(1-prop))*pp_fraction_N2;
                            conc_sO2[U(i,j,k)]+=(mass_added_per_pt*(1-prop))*pp_fraction_O2;
                        }
                    }
                }
            }
        }
        //
        //recalculate the mass conservation and push it back into the appropriate member variable of tissue
        std::vector<double> full_solution=mass_conservation(tt);
        // std::cout<<tissue_mass_current<<" halfway "<<full_solution[0]<<" full solution\n";
        //     total_bubble_mass.push_back(full_solution[1]);
        //      total_mass_tissue.push_back(full_solution[0]);
        //  std::cout<<total_mass_tissue[tt]<<" mass post correction \n";
    }
}
//


//********************************************************************
//**************                              ************************
//**************   conservation mass  2        ************************
//**************                              ************************
//********************************************************************
std::vector<double> Tissue::mass_bounding_box(int tt,int bb){
    double mass_tissue=0;
    double mass_bubbles=0;
    int num_pts_tissue=0;
    int num_pts_bub=0;
    double tissue_volume=0;
    int bound_pt_no_tiss=0;
    int bound_pt_no_bub=0;
    double pts=0;
    double mass_tissue2=0;
    double mass_bubbles2=0;
    int num_pts_tissue2=0;
    int num_pts_bub2=0;
    double tissue_volume2=0;
    int bound_pt_no_tiss2=0;
    int bound_pt_no_bub2=0;
    double pts2=0;
    
    // double prop=1.0;
    
    int rad=ceil(get_bub_radius(bb,tt-1));
    //std::cout<<bubble_number<<" bubble number is \n";
    for(int i=floor(get_bub_locx(bb))-rad;i<=ceil(get_bub_locx(bb))+rad;i++) for (int j=floor(get_bub_locy(bb))-rad;j<=ceil(get_bub_locy(bb))+rad;j++) for (int k=floor(get_bub_locz(bb))-rad;k<=ceil(get_bub_locz(bb))+rad;k++){
        pts+=1;
        if(bub_label[U(i,j,k)]==0){
            num_pts_tissue+=1;
            if(boundary_label[U(i,j,k)]!=0){
                // std::cout<<unsigned(boundary_label[U(i, j, k)])<<" boundary label of point "<<i<<","<<j<<","<<k<<"\n";
                double prop=(theBubbles[bb]->surf_to_point(i,j,k,tt-1)); // if it is in the tissue it will always be less than half the grid box (otherwise the center will be in the bubble)
                if(prop<1.0){
                    bound_pt_no_tiss+=1;
                    // std::cout<<prop<<" prop\n";
                    mass_tissue += ((1-prop)*1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                    tissue_volume+=1-prop;
                    //  std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";
                }
                else{mass_tissue += (1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                    //  tissue_volume+=1.0;
                }
                // std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";}
            }
            
            else{mass_tissue += (1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                //  tissue_volume+=1.0;
                //  std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";
            }
            
        }
        
        else if(bub_label[U(i,j,k)]){
            
            if(boundary_label[U(i,j,k)]!=0){
                double prop=(theBubbles[bb]->surf_to_point(i,j,k,tt-1));
                if(prop<1.0){
                    bound_pt_no_bub+=1;
                    
                    std::vector<double> dists{(theBubbles[bb]->xdist(i,j,k)),theBubbles[bb]->ydist(i,j,k),theBubbles[bb]->zdist(i,j,k)};
                    //these are temporary concentrations that will be asigned to the boundary bubble points, they take the concentration of which ever node they are nearest to.
                    double massconcN2;
                    double massconcO2;
                    int index = 0;
                    
                    for(int i = 1; i < 3; i++)
                    {
                        if(fabs(dists[i]) < fabs(dists[index]))index = i;
                    }
                    
                    
                    if(index==0){
                        if(dists[0]<0.0){
                            massconcN2=conc_sN2[U(i+1,j,k)];
                            massconcO2=conc_sO2[U(i+1,j,k)];
                        }
                        
                        else{
                            massconcN2=conc_sN2[U(i-1,j,k)];
                            massconcO2=conc_sO2[U(i-1,j,k)];
                        }
                    }
                    
                    else if(index==1){
                        if(dists[1]<0.0){
                            massconcN2=conc_sN2[U(i,j+1,k)];
                            massconcO2=conc_sO2[U(i,j+1,k)];
                        }
                        else{
                            massconcN2=conc_sN2[U(i,j-1,k)];
                            massconcO2=conc_sO2[U(i,j-1,k)];
                        }
                    }
                    else if(index==2){
                        if(dists[2]<0.0){
                            massconcN2=conc_sN2[U(i,j,k+1)];
                            massconcO2=conc_sO2[U(i,j,k+1)];
                        }
                        else{massconcN2=conc_sN2[U(i,j,k-1)];
                            massconcO2=conc_sO2[U(i,j,k-1)];
                        }
                    }
                    
                    
                    mass_tissue += ((prop)*1)*(massconcN2 + massconcO2);
                    
                    
                    tissue_volume+=prop;
                    //std::cout<<prop<<" prop\n";
                }
            }
            
            num_pts_bub+=1;
        }
    }
    
    
    compute_label_bub_total(tt);
    
    for(int i=floor(get_bub_locx(bb))-rad;i<=ceil(get_bub_locx(bb))+rad;i++) for (int j=floor(get_bub_locy(bb))-rad;j<=ceil(get_bub_locy(bb))+rad;j++) for (int k=floor(get_bub_locz(bb))-rad;k<=ceil(get_bub_locz(bb))+rad;k++){
        pts2+=1;
        if(bub_label[U(i,j,k)]==0){
            num_pts_tissue2+=1;
            if(boundary_label[U(i,j,k)]!=0){
                // std::cout<<unsigned(boundary_label[U(i, j, k)])<<" boundary label of point "<<i<<","<<j<<","<<k<<"\n";
                double prop=(theBubbles[bb]->surf_to_point(i,j,k)); // if it is in the tissue it will always be less than half the grid box (otherwise the center will be in the bubble)
                if(prop<1.0){
                    bound_pt_no_tiss2+=1;
                    //std::cout<<prop<<" prop\n";
                    mass_tissue2 += ((1-prop)*1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                    tissue_volume2+=1-prop;
                    //  std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";
                }
                else{mass_tissue2 += (1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                    //  tissue_volume2+=1.0;}
                }
                // std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";}
            }
            
            else{mass_tissue2 += (1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                // tissue_volume2+=1.0;
                //  std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";
            }
            
        }
        
        else if(bub_label[U(i,j,k)]){
            
            if(boundary_label[U(i,j,k)]!=0){
                double prop=(theBubbles[bb]->surf_to_point(i,j,k));
                if(prop<1.0){
                    bound_pt_no_bub2+=1;
                    std::vector<double> dists{(theBubbles[bb]->xdist(i,j,k)),theBubbles[bb]->ydist(i,j,k),theBubbles[bb]->zdist(i,j,k)};
                    //these are temporary concentrations that will be asigned to the boundary bubble points, they take the concentration of which ever node they are nearest to.
                    double massconcN2;
                    double massconcO2;
                    int index = 0;
                    
                    for(int i = 1; i < 3; i++)
                    {
                        if(fabs(dists[i]) < fabs(dists[bb]))index = i;
                    }
                    
                    
                    if(index==0){
                        if(dists[0]<0.0){
                            massconcN2=conc_sN2[U(i+1,j,k)];
                            massconcO2=conc_sO2[U(i+1,j,k)];
                        }
                        
                        else{
                            massconcN2=conc_sN2[U(i-1,j,k)];
                            massconcO2=conc_sO2[U(i-1,j,k)];
                        }
                    }
                    
                    else if(index==1){
                        if(dists[1]<0.0){
                            massconcN2=conc_sN2[U(i,j+1,k)];
                            massconcO2=conc_sO2[U(i,j+1,k)];
                        }
                        else{
                            massconcN2=conc_sN2[U(i,j-1,k)];
                            massconcO2=conc_sO2[U(i,j-1,k)];
                        }
                    }
                    else if(index==2){
                        if(dists[2]<0.0){
                            massconcN2=conc_sN2[U(i,j,k+1)];
                            massconcO2=conc_sO2[U(i,j,k+1)];
                        }
                        else{massconcN2=conc_sN2[U(i,j,k-1)];
                            massconcO2=conc_sO2[U(i,j,k-1)];
                        }
                    }
                    
                    
                    mass_tissue2 += ((prop)*1)*(massconcN2 + massconcO2);
                    
                    
                    tissue_volume2+=prop;
                    //std::cout<<prop<<" prop\n";
                }
            }
            
            num_pts_bub2+=1;
        }
    }
    
    //  std::cout<<num_pts_bub2+num_pts_tissue2<<" "<<num_pts_bub+num_pts_tissue<<" num pts\n";
    //  total_mass_tissue.push_back(mass_tissue);
    //  total_tissue_volume.push_back(tissue_volume2);
    //   std::cout<<tissue_volume2 - tissue_volume<<" vols 2 and 1 "<<(4.0/3.0*PI)*pow(get_bub_radius(bb,tt),3)-((4.0/3.0*PI)*pow(get_bub_radius(bb,tt-1),3))<<" bubble volume diff\n";
    //    std::cout<<(get_bub_massN2(bb, tt-1)+get_bub_massO2(bb,tt-1))- (get_bub_massN2(bb, tt)+get_bub_massO2(bb,tt))<<" bub mass diff \nflux "<<theBubbles[bb]->get_dcdr()*4*PI*pow(theBubbles[bb]->get_radius(tt-1),2)*(diffusivity2/(rbar*rbar))*time_step<<"\n"<<mass_tissue-mass_tissue2<<" tissue diff \n";
    double mass_difference= (mass_tissue2 + get_bub_massN2(bb, tt)+get_bub_massO2(bb,tt))-(mass_tissue+get_bub_massN2(bb, tt-1)+get_bub_massO2(bb,tt-1));
    //   std::cout<<mass_difference<<" difference\n";
    double prop_bound_pts=static_cast<double>(bound_pt_no_tiss)/(static_cast<double>(bound_pt_no_tiss)+static_cast<double>(bound_pt_no_bub));
    
    
    
    std::vector<double> masses={-mass_difference,static_cast<double>(bound_pt_no_tiss2+bound_pt_no_bub2),double(tissue_volume)};
    return masses;
}


//********************************************************************
//**************                              ************************
//**************   conservation mass  2        ************************
//**************                              ************************
//********************************************************************

std::vector<double> Tissue::mass_conservation(int tt){
    double mass_tissue=0;
    double mass_bubbles=0;
    int num_pts_tissue=0;
    int num_pts_bub=0;
    double tissue_volume=0;
    int bound_pt_no_tiss=0;
    int bound_pt_no_bub=0;
    double pts=0;
    // double prop=1.0;
    
    for(int i=0; i<node_number[0]; i++) for (int j=0; j<node_number[1]; j++) for (int k=0; k<node_number[2]; k++){
        pts+=1;
        if(bub_label[U(i,j,k)]==0){
            num_pts_tissue+=1;
            if(boundary_label[U(i,j,k)]!=0){
                // std::cout<<unsigned(boundary_label[U(i, j, k)])<<" boundary label of point "<<i<<","<<j<<","<<k<<"\n";
                double prop=(theBubbles[find_bubble_index(boundary_label[U(i,j,k)])]->surf_to_point(i,j,k)); // if it is in the tissue it will always be less than half the grid box (otherwise the center will be in the bubble)
                if(prop<1.0){
                    bound_pt_no_tiss+=1;
                    //std::cout<<prop<<" prop\n";
                    mass_tissue += ((1-prop)*1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                    tissue_volume+=1-prop;
                    //  std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";
                }
                else{mass_tissue += (1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                    tissue_volume+=1.0;}
                // std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";}
            }
            
            else{mass_tissue += (1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
                tissue_volume+=1.0;
                //  std::cout<<conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]<<" "<<i<<" "<<j<<" "<<k<<"\n";
            }
            
        }
        
        else if(bub_label[U(i,j,k)]){
            
            if(boundary_label[U(i,j,k)]!=0){
                double prop=(theBubbles[find_bubble_index(boundary_label[U(i,j,k)])]->surf_to_point(i,j,k));
                if(prop<1.0){
                    bound_pt_no_bub+=1;
                    int indx=find_bubble_index(boundary_label[U(i,j,k)]);
                    std::vector<double> dists{(theBubbles[indx]->xdist(i,j,k)),theBubbles[indx]->ydist(i,j,k),theBubbles[indx]->zdist(i,j,k)};
                    //these are temporary concentrations that will be asigned to the boundary bubble points, they take the concentration of which ever node they are nearest to.
                    double massconcN2;
                    double massconcO2;
                    int index = 0;
                    
                    for(int i = 1; i < 3; i++)
                    {
                        if(fabs(dists[i]) < fabs(dists[index]))index = i;
                    }
                    
                    
                    if(index==0){
                        if(dists[0]<0.0){
                            massconcN2=conc_sN2[U(i+1,j,k)];
                            massconcO2=conc_sO2[U(i+1,j,k)];
                        }
                        
                        else{
                            massconcN2=conc_sN2[U(i-1,j,k)];
                            massconcO2=conc_sO2[U(i-1,j,k)];
                        }
                    }
                    
                    else if(index==1){
                        if(dists[1]<0.0){
                            massconcN2=conc_sN2[U(i,j+1,k)];
                            massconcO2=conc_sO2[U(i,j+1,k)];
                        }
                        else{
                            massconcN2=conc_sN2[U(i,j-1,k)];
                            massconcO2=conc_sO2[U(i,j-1,k)];
                        }
                    }
                    else if(index==2){
                        if(dists[2]<0.0){
                            massconcN2=conc_sN2[U(i,j,k+1)];
                            massconcO2=conc_sO2[U(i,j,k+1)];
                        }
                        else{massconcN2=conc_sN2[U(i,j,k-1)];
                            massconcO2=conc_sO2[U(i,j,k-1)];
                        }
                    }
                    
                    
                    mass_tissue += ((prop)*1)*(massconcN2 + massconcO2);
                    
                    
                    tissue_volume+=prop;
                    //std::cout<<prop<<" prop\n";
                }
            }
            
            num_pts_bub+=1;
        }
    }
    
    
    
    
    
    total_mass_tissue.push_back(mass_tissue);
    total_tissue_volume.push_back(tissue_volume);
    
    double bubble_mass_ind[number_of_bubbles];
    for (int bb=0;bb<number_of_bubbles; bb++){
        bubble_mass_ind[bb]=theBubbles[bb]->get_mmassN2(tt)+theBubbles[bb]->get_mmassO2(tt);
        mass_bubbles += bubble_mass_ind[bb];
    }
    
    double prop_bound_pts=static_cast<double>(bound_pt_no_tiss)/(static_cast<double>(bound_pt_no_tiss)+static_cast<double>(bound_pt_no_bub));
    
    total_bubble_mass.push_back(mass_bubbles);
    
    //    if(tt==!0){
    //       std::cout<<(mass_tissue/mbar)-theBubbles[0]->get_dcdr()<<" tissue mass plus dc/dr\n";
    //    }
    //  std::cout<<tt<<","<<num_pts_tissue<<","<<num_pts_bub<<","<<pts<<","<<tissue_volume<<" mass con\n";
    
    
    std::vector<double> masses={mass_tissue,mass_bubbles,static_cast<double>(bound_pt_no_tiss+bound_pt_no_bub),double(prop_bound_pts)};
    return masses;
}




////********************************************************************
////**************                              ************************
////**************   update boundary_pts mass   ************************
////**************                              ************************
////********************************************************************
//void Tissue::update_boundary_pts(int tt){
//
//    std::vector<double> halfway_solution=mass_conservation(tt); //call mass conservation initially to find out the mass difference.
//
//    double tissue_mass_current=halfway_solution[0];
//    double no_bound_pts=halfway_solution[2];
//    double prop_bound_pts=halfway_solution[3];
//
//    for(int bb=0;bb<number_of_bubbles; bb++){
//        //   //work out the mass deficit if any
//        // std::cout<<no_bound_pts<<" boundary pts " <<prop_bound_pts<<" prop in tissue\n";
//        double mass_dc_dr=theBubbles[bb]->get_dcdr()*4*PI*pow(theBubbles[bb]->get_radius(tt-1),2)*(diffusivity/(rbar*rbar))*time_step;
//        double bub_change=(total_bubble_mass[tt-1] + mass_dc_dr)-(halfway_solution[1]+mass_dc_dr);
//        double gas_difference=(((-tissue_mass_current+total_mass_tissue[tt-1])-bub_change))/no_bound_pts;
//      //  std::cout<<mass_dc_dr<<" dc_dr mass,   "<<total_bubble_mass[tt-1]<<" bubble mass,   "<<bub_change<<" change in bubble mass,   "<<total_mass_tissue[tt-1]<<" tissue mass,   "<<gas_difference<<" change in tissue \n";
//
//       // add the deficit to the boundary points
//        for(int i=floor(theBubbles[bb]->get_locations(0))-theBubbles[bb]->get_radius();i<=ceil(theBubbles[bb]->get_locations(0))+theBubbles[bb]->get_radius();i++){
//            for(int j=floor(theBubbles[bb]->get_locations(1))-theBubbles[bb]->get_radius();j<=ceil(theBubbles[bb]->get_locations(1))+theBubbles[bb]->get_radius();j++){
//                for(int k=floor(theBubbles[bb]->get_locations(2))-theBubbles[bb]->get_radius();k<=ceil(theBubbles[bb]->get_locations(2))+theBubbles[bb]->get_radius();k++){
//                    if(boundary_label[U(i,j,k)]!=0){
//
//                        conc_sN2[U(i,j,k)]+=(gas_difference*pp_fraction_N2);
//                        conc_sO2[U(i,j,k)]+=gas_difference*pp_fraction_O2;
//
//                        if(bub_label!=0){
//                            std::vector<double> dists{(theBubbles[bb]->xdist(i,j,k)),theBubbles[bb]->ydist(i,j,k),theBubbles[bb]->zdist(i,j,k)};
//                            //these are temporary concentrations that will be asigned to the boundary bubble points, they take the concentration of which ever node they are nearest to.
//                            int index = 0;
//
//                            for(int i = 1; i < 3; i++)
//                            {
//                                if(fabs(dists[i]) < fabs(dists[index]))index = i;
//                            }
//
//
//                            if(index==0){
//                                if(dists[0]<0.0){
//                                    conc_sN2[U(i+1,j,k)]+=gas_difference*pp_fraction_N2;
//                                    conc_sO2[U(i+1,j,k)]+=gas_difference*pp_fraction_O2;
//                                }
//
//                                else{
//                                    conc_sN2[U(i-1,j,k)]+=gas_difference*pp_fraction_N2;
//                                    conc_sO2[U(i-1,j,k)]+=gas_difference*pp_fraction_O2;
//                                }
//                            }
//
//                            else if(index==1){
//                                if(dists[1]<0.0){
//                                    conc_sN2[U(i,j+1,k)]+=gas_difference*pp_fraction_N2;
//                                    conc_sO2[U(i,j+1,k)]+=gas_difference*pp_fraction_O2;
//                                }
//                                else{
//                                    conc_sN2[U(i,j-1,k)]+=gas_difference*pp_fraction_N2;
//                                    conc_sO2[U(i,j-1,k)]+=gas_difference*pp_fraction_O2;
//                                }
//                            }
//                            else if(index==2){
//                                if(dists[2]<0.0){
//                                    conc_sN2[U(i,j,k+1)]+=gas_difference*pp_fraction_N2;
//                                    conc_sO2[U(i,j,k+1)]+=gas_difference*pp_fraction_O2;
//                                }
//                                else{conc_sN2[U(i,j,k-1)]+=gas_difference*pp_fraction_N2;
//                                    conc_sO2[U(i,j,k-1)]+=gas_difference*pp_fraction_O2;
//                                }
//                            }
//
//
//                        }
//
//                    }
//                }
//            }
//        }
//    }
//
//    //recalculate the mass conservation and push it back into the appropriate member variable of tissue
//
//    std::vector<double> full_solution=mass_conservation(tt);
//    // std::cout<<tissue_mass_current<<" halfway "<<full_solution[0]<<" full solution\n";
//    total_bubble_mass.push_back(full_solution[1]);
//    total_mass_tissue.push_back(full_solution[0]);
//}
//
//
//
////********************************************************************
////**************                              ************************
////**************   conservation mass          ************************
////**************                              ************************
////********************************************************************
//std::vector<double> Tissue::mass_conservation(int tt){
//    double mass_tissue=0;
//    double mass_bubbles=0;
//    int num_pts_tissue=0;
//    int num_pts_bub=0;
//    double tissue_volume=0;
//    int bound_pt_no_tiss=0;
//    int bound_pt_no_bub=0;
//    double pts=0;
//    // double prop=1.0;
//    for(int i=0; i<node_number[0]; i++) for (int j=0; j<node_number[1]; j++) for (int k=0; k<node_number[2]; k++){
//        pts+=1;
//        if(bub_label[U(i,j,k)]==0){
//            num_pts_tissue+=1;
//            if(boundary_label[U(i,j,k)]!=0){
//                // std::cout<<unsigned(boundary_label[U(i, j, k)])<<" boundary label of point "<<i<<","<<j<<","<<k<<"\n";
//                double prop=(theBubbles[find_bubble_index(boundary_label[U(i,j,k)])]->surf_to_point(i,j,k)); // if it is in the tissue it will always be less than half the grid box (otherwise the center will be in the bubble)
//                if(prop<1.0){
//                    bound_pt_no_tiss+=1;
//                    //std::cout<<prop<<" prop\n";
//                    mass_tissue += ((1-prop)*1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
//                    tissue_volume+=1-prop;
//                }
//                else{mass_tissue += (1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
//                    tissue_volume+=1.0;}
//            }
//
//            else{mass_tissue += (1)*(conc_sN2[U(i, j, k)] + conc_sO2[U(i, j, k)]);
//                tissue_volume+=1.0;
//            }
//        }
//
//        else if(bub_label[U(i,j,k)]){
//
//            if(boundary_label[U(i,j,k)]!=0){
//                double prop=(theBubbles[find_bubble_index(boundary_label[U(i,j,k)])]->surf_to_point(i,j,k));
//                if(prop<1.0){
//                    bound_pt_no_bub+=1;
//                    int indx=find_bubble_index(boundary_label[U(i,j,k)]);
//                    std::vector<double> dists{(theBubbles[indx]->xdist(i,j,k)),theBubbles[indx]->ydist(i,j,k),theBubbles[indx]->zdist(i,j,k)};
//                    //these are temporary concentrations that will be asigned to the boundary bubble points, they take the concentration of which ever node they are nearest to.
//                    double massconcN2;
//                    double massconcO2;
//                    int index = 0;
//
//                    for(int i = 1; i < 3; i++)
//                    {
//                        if(fabs(dists[i]) < fabs(dists[index]))index = i;
//                    }
//
//
//                    if(index==0){
//                        if(dists[0]<0.0){
//                            massconcN2=conc_sN2[U(i+1,j,k)];
//                            massconcO2=conc_sO2[U(i+1,j,k)];
//                        }
//
//                        else{
//                            massconcN2=conc_sN2[U(i-1,j,k)];
//                            massconcO2=conc_sO2[U(i-1,j,k)];
//                        }
//                    }
//
//                    else if(index==1){
//                        if(dists[1]<0.0){
//                            massconcN2=conc_sN2[U(i,j+1,k)];
//                            massconcO2=conc_sO2[U(i,j+1,k)];
//                        }
//                        else{
//                            massconcN2=conc_sN2[U(i,j-1,k)];
//                            massconcO2=conc_sO2[U(i,j-1,k)];
//                        }
//                    }
//                    else if(index==2){
//                        if(dists[2]<0.0){
//                            massconcN2=conc_sN2[U(i,j,k+1)];
//                            massconcO2=conc_sO2[U(i,j,k+1)];
//                        }
//                        else{massconcN2=conc_sN2[U(i,j,k-1)];
//                            massconcO2=conc_sO2[U(i,j,k-1)];
//                        }
//                    }
//
//
//                    mass_tissue += ((prop)*1)*(massconcN2 + massconcO2);
//
//
//                    tissue_volume+=prop;
//                    //std::cout<<prop<<" prop\n";
//                }
//            }
//
//            num_pts_bub+=1;
//        }
//    }
//
//
//
//
//
//    // total_mass_tissue.push_back(mass_tissue);
//    total_tissue_volume.push_back(tissue_volume);
//
//    double bubble_mass_ind[number_of_bubbles];
//    for (int bb=0;bb<number_of_bubbles; bb++){
//        bubble_mass_ind[bb]=theBubbles[bb]->mass_in_bubble(tt);
//        mass_bubbles += bubble_mass_ind[bb];
//    }
//
//    double prop_bound_pts=static_cast<double>(bound_pt_no_tiss)/(static_cast<double>(bound_pt_no_tiss)+static_cast<double>(bound_pt_no_bub));
//    // total_bubble_mass.push_back(mass_bubbles);
//    //    if(tt==!0){
//    //        std::cout<<(mass_tissue/mbar)-theBubbles[0]->get_dcdr()<<" tissue mass plus dc/dr\n";
//    //    }
//    // std::cout<<tt<<","<<num_pts_tissue<<","<<num_pts_bub<<","<<pts<<","<<tissue_volume<<" mass con\n";
//
//    std::vector<double> masses={mass_tissue,mass_bubbles,static_cast<double>(bound_pt_no_tiss+bound_pt_no_bub),double(prop_bound_pts)};
//    return masses;
//}


























