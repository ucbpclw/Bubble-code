//
//  Dive_profile.hpp
//  Bubble Jan 2014
//
//  Created by claire Walsh on 07/01/2014.
//  Copyright (c) 2014 University college london. All rights reserved.
//

#ifndef Bubble_Jan_2014_Dive_profile_hpp
#define Bubble_Jan_2014_Dive_profile_hpp
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "parameters.h"
#include <string.h>
#include <sstream>
extern  std::vector<double> P_amb;


class Dive_profile{
private:
    std::vector<std::vector< double> > P_amb_read_in;
    std::string dive_name;
    
public:
    //constructor
    Dive_profile(std::string a){dive_name=a;};
    //////////////////////////////////////////
    void read_in_profile(double,int);
    double linear_interpolate(double);
    void setP_amb(double,int);
    std::string get_dive_name(){return dive_name;};
    
};

#endif
