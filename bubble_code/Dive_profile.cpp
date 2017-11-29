//
//  Dive_profile.cpp
//  Bubble Jan 2014
//
//  Created by claire Walsh on 07/01/2014.
//  Copyright (c) 2014 University college london. All rights reserved.
//

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "parameters.h"
#include "Dive_profile.hpp"
#include <cassert>


//********************************************************************
//**************                              ************************
//**************      read_in_ profile        ************************
//**************                              ************************
//********************************************************************

void Dive_profile::read_in_profile(double timept,int units){
    std::ifstream input_file ("../dive_profiles/PROFILE.txt");
    
    if(!input_file.is_open()){std::cout<<"error opening pressure input file \n";}
    std::string lineData;
    
    while(getline(input_file, lineData))
    {
        double d;
        std::vector<double> column;
        std::stringstream lineStream(lineData);
        std::cout<<lineData<<"\n";
        while (lineStream >> d)
            column.push_back(d);
        
        // std::cout<<column[0]<<" "<<column[1] <<"\n";
        
        P_amb_read_in.push_back(column);
    }
    // std::cout<<P_amb_read_in.size()<<"\n";
    
    std::cout<<"dive length is "<<P_amb_read_in[P_amb_read_in.size()-1][0]<<"\n";
    std::cout<<"length of dive profile is "<<P_amb_read_in.size()<<"\n";
    assert(P_amb_read_in[P_amb_read_in.size()-1][0]==time_of_dive);
    
    
    setP_amb(timept,units);
}
//********************************************************************
//**************                              ************************
//**************      setP_amb               ************************
//**************                              ************************
//********************************************************************


void Dive_profile::setP_amb(double timestep,int units){
    
    
    //int maxtime=floor(time_of_dive*60)/(timestep);
    //uncomment for seconds
    int maxtime=time_of_dive/timestep;
    
    std::cout<<"time step "<<timestep<<" maxtime is "<<maxtime<<"\n";
    double pressure;
    for (int i=0; i<=maxtime-1; i++){
        
        // pressure=linear_interpolate(i*(timestep/60));
        // uncomment below for seconds
        pressure=linear_interpolate(i*(timestep));
        if(units==1)
        { P_amb.push_back(((pressure*0.10693064)+1.01232)/pbar);} //convert m depth to pressure in bar
        else if (units==2)
        {P_amb.push_back(((pressure*0.0689475729)+1.01232)/pbar);} //pressure is in psi convert to bar
        else if (units==3)
        {P_amb.push_back(pressure/pbar);} //non-dimensinalise P_amb
        
    }
    
    if(units==1)
    { P_amb[maxtime]=((1.01232)/pbar);} //convert m depth to pressure in bar
    else if (units==2)
    {P_amb[maxtime]=((1.01232)/pbar);} //pressure is in psi convert to bar
    else if (units==3)
    {P_amb[maxtime]=(1/pbar);} //pressure is in bar already this
    std::cout<<"P_amb set\n";
    //    for (int i=0; i<=maxtime;i++){
    //    std::cout<<P_amb[i]<<"\n";
    //    }
    //  writing to file
    //    std::ofstream outfile;
    //    outfile.open("/Users/clairewalsh/Dropbox/Modelling/dan_May_2014/dan_May_2014/outputs/sensitvity_analysis/sensitivity_analysis_june_2016/profile.txt", std::ios::app);
    //    if(!outfile.is_open()){std::cout<<"error opening pressure output file \n";}
    //    std::cout << "Writing pressure profile to the file" << std::endl;
    //
    //    for(int i=0; i<=maxtime; i++){
    //        outfile <<P_amb[i]<<std::endl;
    //    }
    //    // ****** close the opened file. **********
    //    outfile.close();
}


//********************************************************************
//**************                              ************************
//**************      linear interpolate      ************************
//**************                              ************************
//********************************************************************


double Dive_profile::linear_interpolate(double timept){
    int i=0;
    double ya;
    double yb;
    double xa;
    double xb;
    double pressure;
    assert(timept>=0);
    //std::cout<<P_amb_read_in[P_amb_read_in.size()-1][0]/60<<" max \n";
    //std::cout<<timept<<" p_amb final term read in is "<<P_amb_read_in[P_amb_read_in.size()-1][0]<<"\n";
    
    assert (timept<=P_amb_read_in[P_amb_read_in.size()-1][0]);// check that the time point you've provided is within the dive time
    while(i < P_amb_read_in.size()){
        //  std::cout<<"P_amb is "<<P_amb_read_in[i][0]<<"\n";
        if(P_amb_read_in[i][0]==timept)
        {
            pressure=P_amb_read_in[i][1];
            //  std::cout<<"time is exactly equal to known point pressure is "<<pressure<<"\n";
            break;
        }
        
        else if(P_amb_read_in[i][0]<timept && P_amb_read_in[i+1][0]>timept)
        {   ya=P_amb_read_in[i][1];
            yb=P_amb_read_in[i+1][1];
            xa=P_amb_read_in[i][0];
            xb=P_amb_read_in[i+1][0];
            //std::cout<<"times either side of "<<timept<<" are "<<xa<<" , "<<xb<<"\n";
            pressure=(ya + (yb - ya) * (timept - xa) / (xb - xa));
            //std::cout<<"pressure caluclated to be "<<pressure<<"\n";
            break;
            
        }
        i++;
    }
    
    //std::cout<<"pressure is "<<pressure<<" time pt is "<<timept<<"\n";
    return pressure;
}



