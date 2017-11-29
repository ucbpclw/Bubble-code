//
//  main.cpp
//  Bubble Jan 2014
//
//  Created by claire Walsh on 07/01/2014.
//  Copyright (c) 2014 University college london. All rights reserved.
//

#include <iostream>
#include <istream>
#include <math.h>
#include <vector>
#include <fstream>
#include "parameters.h"
#include "Dive_profile.hpp"
#include "tissue.hpp"
#include "Bubble.hpp"
#include <string>
#include <cassert>
#include <time.h>
std::vector<double> P_amb;
int main(int argc, const char * argv[]){
    unsigned int time_r=(unsigned int)time(NULL);    // seed the random number generator 5 used instead of time_r for repeatable simulations
    srand(51/*time_r*/);
    
    //Decide on a tissue time and space step
    
    Tissue block(node_num_x*rbar,node_num_y*rbar,node_num_z*rbar,3); //constructor of a tissue block with first term size x, second size y, 3rd size z, 4th number of bubbles.
    // Check that the bit shifts are correct for accessing the conc. grid
    assert(log2(block.get_node_number_x())==bit_shift_i);
    assert(log2(block.get_node_number_x()*block.get_node_number_y())==bit_shift_ij);
    std::cout<<"cbar is "<<cbar<<"\n";
    std::cout<<"grid size is squares is "<<block.get_node_number_x()<<","<<block.get_node_number_y()<<","<<block.get_node_number_z()<<"\n";
    std::cout<<block.get_space_step()<<" is the space step \n";
    std::cout<<node_num_x*rbar<<","<<node_num_y*rbar<<","<<node_num_z*rbar<<" is the tissue size in m \n";
    block.set_time_step(0.025);
    std::cout<<block.get_time_step()<<" time step (s) \n";
    std::cout<<diffusivity<< "diffusivity\n";
    clock_t timer;
    //check Stability criterion are met
    assert(block.von_Neumann());
    
    
    //create and read in the dive profile
    std::cout<<"what units is the pressure profile in? for meters depth type 1 \n for psi type 2 \n for bar type 3\n";
    int units;
    std::cin>>units;
    Dive_profile diveprofile("sens_20min");
    
    std::string const filepath("../outputs/FILENAME.txt");
    block.check_radius_folder(filepath);
    diveprofile.read_in_profile(block.get_time_step(),units);
    int z=3;
    //std::cout<<"select which method you are using 1 for lapace, 2 for Gernhadt, 3 for Gent \n";
    //std::cin>>z;
    block.initialise_bubbles(z);
    block.initialise_conc_grid();
    block.compute_label_bub_total(0);
    std::vector<double> initial_masses=block.mass_conservation(0);
    block.set_tissue_mass(initial_masses[0]);
    block.set_total_bub_mass(initial_masses[1]);
    //  std::cout<<"time step \t \t \t radius\n";
    for (int t=1; t<time_of_dive/block.get_time_step() /*floor(time_of_dive*60/(block.get_time_step()))*/; t++){
        // std::cout<<"time step "<<t<<"\t";
        timer = clock();
        block.diffusion(t);
        //        timer = clock() - timer;
        //        std::cout<<"time for diffusion was "<< float(timer)/CLOCKS_PER_SEC <<"secs \n";
        // block.bub_dcdr();   //JP's original dc/dr function
        //timer =clock();
        block.bub_dc_dr_spherical(t); //dc/dr across bubble surface in spherical co-ordinates
        //timer=clock()-timer;
        //std::cout<<"time for dc/dr was "<< float(timer)/CLOCKS_PER_SEC <<"secs \n";
        
        if (z==1){
            block.bub_drdt_lap(t);        //dr/dt with no elasticity
        }
        else if (z==2){
            block.bub_drdt_gern(t);       //dr/dt with Gern elasticity
        }
        else if( z==3){
            // timer=clock();
            block.bub_drdt(t);            //dr/dt with gent elasticity
        }
        
        else {std::cout<<"you have not picked an elasticity method try again \n";
            assert(z==3);};
        
        
        //timer=clock()-timer;
        //std::cout<<"time for dr/dt was "<< float(timer)/CLOCKS_PER_SEC <<"secs \n";
        //output the total mass in bubble and tisue
        // block.mass_conservation(t);
        block.update_boundary_pts(t);//check the bubble conc and radius (now called in mass conservation.
        //tissue diffusion
        if(t%200==0){
            //std::cout<<t<<"\t\t\t"<<block.get_bub_radius(0)<<"\n";
            block.param_check(t);
        }
        
        block.bubble_check(filepath, t,z); //Check bubble is not too near the edge or overlaps with others
    }
    
    block.write_radius_data(diveprofile.get_dive_name(),filepath);
    //delete the bubbles
    block.delete_theBubbles();
    block.delete_conc_grid();
    return 0;
}

