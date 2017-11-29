//
//  Bubble.cpp
//  Bubble Jan 2014
//
//  Created by claire Walsh on 07/01/2014.
//  Copyright (c) 2014 University college london. All rights reserved.
//

#include <iostream>
#include "Bubble.hpp"
#include <string>
#include <cassert>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include "parameters.h"
#include <vector>
#include "Dive_profile.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <limits>
#include <iomanip>

//********************************************************************
//**************                              ************************
//**************      randomise bubbles       ************************
//**************                              ************************
//********************************************************************
// set the initial bubble sizes and positions the min and max are for positioning the bubble to ensure it is not so close the the edge that dc/dr can't be calculated
void Bubble::randomise_bubble(int* node_number){
    std::cout<<"number of node is "<<*node_number<<","<<*(node_number+1)<<","<<*(node_number+2)<<"\n";
    std::vector<double> temploc;
    int maxint[3]={static_cast<int>(node_number[0] - 3- min_radius),static_cast<int>((node_number[1] -3- min_radius)),static_cast<int>((node_number[2]- 3 -min_radius))};
    int minint[3]={static_cast<int>(3 + min_radius),static_cast<int>(3 + min_radius),static_cast<int>(3 + min_radius)};
    
    double min[3]={static_cast<double>(minint[0]),static_cast<double>(minint[1]),static_cast<double>(minint[2])};// convert to doubles so the bubble locations can be non-grid centered
    double max[3]={static_cast<double>(maxint[0]),static_cast<double>(maxint[1]),static_cast<double>(maxint[2])};
    
    temploc={min[0] + (rand() / ( RAND_MAX / (max[0]-min[0]) ) ),min[1] + (rand() / ( RAND_MAX / (max[1]-min[1]) ) ),min[2] + (rand() / ( RAND_MAX / (max[2]-min[2]) ) )};
    
    //            std::cout<<min + (rand() / div)<<"\n";
    
    set_locations(temploc[0], temploc[1], temploc[2]);
    //     std::cout<<locations[0]<<","<<locations[1]<<","<<locations[2]<<"\n";
    
    
}

//********************************************************************
//**************                              ************************
//**************      compute bubble label    ************************
//**************                              ************************
//********************************************************************
// points in bubble are asigned bubble number in the mask otherwise they are zero
void Bubble::compute_bubble_label(uint8* label,int tt){
    // std::cout<<radius.back()<<" rad at start of compute bubble label\n";
    int rad=ceil(radius.back());
    int pts=0;
    //std::cout<<bubble_number<<" bubble number is \n";
    for(int i=floor(locations[0])-rad;i<=ceil(locations[0])+rad;i++){
        for(int j=floor(locations[1])-rad;j<=ceil(locations[1])+rad;j++){
            for(int k=floor(locations[2])-rad;k<=ceil(locations[2])+rad;k++){
                label[U(i,j,k)]=
                ((SQUARED(i-locations[0])+SQUARED(j-locations[1])+SQUARED(k-locations[2]))<=SQUARED(radius.back()))? bubble_number : 0;
                
                //prints out pts in bubbles
                //                if (label[U(i,j,k)]){
                //                    pts=pts+1;
                //                    std::cout<<i<<","<<j<<","<<k<<"\n";
                //                }
            }
        }          //if the point less than or equal to the radius assign it bubble number otherwise let it be zero
    }
    //  std::cout<<pts<<"\n";
}
//********************************************************************
//**************                              ************************
//**************      compute boundary label  ************************
//**************                              ************************
//********************************************************************
//sets the boundary label mask, any points in the tissue but at the bubble bounrary should have value 1 everything else should be zero.
void Bubble::compute_boundary_label(uint8* boundary_label,uint8* bub_label,int tt){
    int rad=ceil(radius.back())+1; //to make sure we get all boundary points
    for(int i=floor(locations[0])-rad;i<=ceil(locations[0])+rad;i++)for(int j=floor(locations[1])-rad;j<=ceil(locations[1])+rad;j++)for(int k=floor(locations[2])-rad;k<=ceil(locations[2])+rad;k++){
        //these are negative
        if(bub_label[U(i+1,j,k)]!=bub_label[U(i,j,k)] or bub_label[U(i-1,j,k)]!=bub_label[U(i,j,k)]
           or  bub_label[U(i,j+1,k)]!=bub_label[U(i,j,k)] or bub_label[U(i,j-1,k)]!=bub_label[U(i,j,k)]
           or bub_label[U(i,j,k+1)]!=bub_label[U(i,j,k)] or bub_label[U(i,j,k-1)]!=bub_label[U(i,j,k)])
        {boundary_label[U(i,j,k)]=bubble_number;}
        
        //prints out boundary pts
        // std::cout<<i<<","<<j<<","<<k<<"\n";
        else{ boundary_label[U(i,j,k)]=0;}
        
    }
}
//********************************************************************
//**************                              ************************
//**************     mass of gas in bubble    ************************
//**************                              ************************
//********************************************************************
//This calcultes the mass of gas in the bubble at a particular time point via updating the mass from the mass flux, inital mass is caluclated with Gent formulation for bubble pressure then using perfect gas formular, this mass is then used to set the concs for the bubble boundary conditions!!
double Bubble::mass_in_bubble(int tt,double timestep){
    double bubble_gas_massN2;
    double bubble_gas_massO2;
    // std::cout<<tt<<"time point\n";
    if (tt==0){
        double rad_initial=initial_radius;////1e-6/rbar;//min_radius;//
        
        double initial_pressure=P_amb[tt]+(2.0*sigma/radius.back())+(mu/2)*(5.0-4.0*(rad_initial/radius.back())-pow((rad_initial/radius.back()),4));
        
        bubble_gas_massN2=(pp_fraction_N2*initial_pressure*(4.0/3.0)*PI*pow(get_radius(tt),3))/L_N2;
        bubble_gas_massO2=(pp_fraction_O2*initial_pressure*(4.0/3.0)*PI*pow(get_radius(tt),3))*(MrO2/(MrN2*L_N2));
    }
    
    else{
        bubble_gas_massN2=(get_dcdrN2()*4*PI*pow(get_radius(tt-1),2)*(diffusivity2/(rbar*rbar))*timestep)+bubmassN2[tt-1];
        bubble_gas_massO2=(get_dcdrO2()*4*PI*pow(get_radius(tt-1),2)*(diffusivity2/(rbar*rbar))*timestep)+bubmassO2[tt-1];
    }
    
    
    bubmassN2.push_back(bubble_gas_massN2);
    bubmassO2.push_back(bubble_gas_massO2);
    bubbleconc(tt);
    return (bubble_gas_massN2 + bubble_gas_massO2);
}


//********************************************************************
//**************                              ************************
//**************      conc_in_bubble          ************************
//**************                              ************************
//********************************************************************
//concentration inside the bubble caluclated from mass and Henry's Law
void Bubble::bubbleconc(int tt){
    // double rad_initial=initial_radius;////1e-6/rbar;//min_radius;//
    bubconcN2.push_back(get_mmassN2()*L_N2/((4.0*PI/3.0)*pow(get_radius(tt),3)));
    bubconcO2.push_back(get_mmassO2()*L_O2/((4.0*PI/3.0)*pow(get_radius(tt),3)));
    
    //  bubpress.push_back(P_amb[tt]+(2.0*sigma/radius.back())+(mu/2)*(5.0-4.0*(rad_initial/radius.back())-pow((rad_initial/radius.back()),4));
}



//********************************************************************
//**************                              ************************
//**************      dc/dr                   ************************
//**************                              ************************
//********************************************************************
////calculate the dc/dr across the bubble boundary using JP's equation for second order finite difference on cartesian grid.
//void Bubble::dcdr(double *conc,double tissue_step){
//    double node_dist;
//    double node_step;
//
//    //the distance between the edge of the bubble and the node in the surrounding tissue that is closest is given by: (ceil applies only to the first bracket)
//    node_dist=ceil((radius.back()*rbar/tissue_step))-(radius.back()*rbar/tissue_step);
//   // std::cout<<"node dist is "<<ceil((radius.back()*rbar/tissue_step))<<"\n";
//    if(node_dist==0){
//        node_dist=1.0;
//
//    }
//    // If the node is too close then the 2nd order finite difference will become very inaccurate, sensitivty is set in parameters.h
//    if(node_dist>=DC_DR_SENS)node_step = 1;
//    else{
//        node_step = 2;   //if it is too close use the node beyond it to calculate FD
//        node_dist += 1.0;} //in that case node distance will incease one node step.
//
//    int rad_temp=(int)floor(radius.back()*rbar/tissue_step);
//
//       //this is dc/dr as from JP's code. The points in the bubconc matrix are one and two nodes in or out from the bubble radius at the 6 corners of the bubble. It is [x][y][z][0] because the bubconc matrix has already had old turned to new.
//    double  dcdr=
//    //in the xplane
//    (-(1+2*node_dist)*bubconc.back()+
//     pow(1+node_dist,2.0)*conc[locations[0]+rad_temp+node_step][locations[1]][locations[2]][0]-
//     node_dist*node_dist*conc[locations[0]+rad_temp+node_step+1.0][locations[1]][locations[2]][0]-
//
//     (1+2*node_dist)*bubconc.back()+
//     pow(1+node_dist,2.0)*conc[locations[0]-rad_temp-node_step][locations[1]][locations[2]][0]-
//     node_dist*node_dist*conc[locations[0]-rad_temp-node_step-1][locations[1]][locations[2]][0] -
//
//     //in the y-plane
//     (1+2*node_dist)*bubconc.back()+
//     pow(1+node_dist,2.0)*conc[locations[0]][locations[1]+rad_temp+node_step][locations[2]][0]-
//     node_dist*node_dist*conc[locations[0]][locations[1]+rad_temp+node_step+1][locations[2]][0] -
//
//     (1+2*node_dist)*bubconc.back()+
//     pow(1+node_dist,2.0)*conc[locations[0]][locations[1]-rad_temp-node_step][locations[2]][0]-
//     node_dist*node_dist*conc[locations[0]][locations[1]-rad_temp-node_step-1][locations[2]][0]-
//
//     //in the z-plane
//     (1+2*node_dist)*bubconc.back()+
//     pow(1+node_dist,2.0)*conc[locations[0]][locations[1]][locations[2]+rad_temp+node_step][0]-
//     node_dist*node_dist*conc[locations[0]][locations[1]][locations[2]+rad_temp+node_step+1][0] -
//
//     (1+2*node_dist)*bubconc.back()+
//     pow(1+node_dist,2.0)*conc[locations[0]][locations[1]][locations[2]-rad_temp-node_step][0]-
//     node_dist*node_dist*conc[locations[0]][locations[1]][locations[2]-rad_temp-node_step-1][0])/
//    6.0/(node_dist*(1+node_dist)*tissue_step/rbar);

//i'm not sure if there should be an rbar here or not??
//
//    dc_dr.push_back(dcdr);
//   // std::cout<<"dc/dr= "<<dc_dr.back()<<"\n";
//
//}
//********************************************************************
//************                                    ********************
//************   dist of point to bubble surface  ********************
//************                                    ********************
//********************************************************************
double Bubble::surf_to_point(int xx, int yy, int zz){
    //converts the point into spherical co-ordinates then moves along r direction till it meets the bubble surface, the surface point is then converted back to cartesian co-ordinates.
    std::vector<double>dist;
    double x=(xx-get_locations(0));
    double y=(yy-get_locations(1));
    double z=(zz-get_locations(2));
    double theta=acos(z/sqrt((x*x)+(y*y)+(z*z)));
    double phi=atan2(y,x);
    double rad_dist=sqrt(x*x+y*y+z*z);
    
    int face=which_face(x,y,z);
    
    //std::cout<<theta<<" theta\n"<<phi<<" phi\n";
    double  surf_x=get_radius()*sin(theta)*cos(phi);
    double surf_y=get_radius()*sin(theta)*sin(phi);
    double surf_z=get_radius()*cos(theta);
    
    
    double diffs[3]={(x-surf_x),(y-surf_y),(z-surf_z)};
    
    const double tol=0.5e-15;
    for (int i=0;i<3;i++){
        if(fabs(diffs[i])<tol){diffs[i]=0;}
        dist.push_back(sqrt(diffs[i]*diffs[i]));
        
    }
    double length=length_line_square(theta*(180/PI),phi*(180/PI),face);
    double radial_dist=fabs(rad_dist-get_radius());
    
    if ((length/2.0)<radial_dist){return 1.0;} // the actual boundary lies in the other cell.
    else{//std::cout<<((length/2.0)-radial_dist)/length<<"\n";
        return(((length/2.0)-radial_dist)/length); }
    
    
    
    
    
}

//********************************************************************
//************                                    ********************
//************   dist of point to bubble surface  ********************
//************                                    ********************
//********************************************************************
double Bubble::surf_to_point(int xx, int yy, int zz, int tt){
    //converts the point into spherical co-ordinates then moves along r direction till it meets the bubble surface, the surface point is then converted back to cartesian co-ordinates.
    std::vector<double>dist;
    double x=(xx-get_locations(0));
    double y=(yy-get_locations(1));
    double z=(zz-get_locations(2));
    double theta=acos(z/sqrt((x*x)+(y*y)+(z*z)));
    double phi=atan2(y,x);
    double rad_dist=sqrt(x*x+y*y+z*z);
    
    int face=which_face(x,y,z);
    
    //std::cout<<theta<<" theta\n"<<phi<<" phi\n";
    double  surf_x=get_radius(tt)*sin(theta)*cos(phi);
    double surf_y=get_radius(tt)*sin(theta)*sin(phi);
    double surf_z=get_radius(tt)*cos(theta);
    
    
    double diffs[3]={(x-surf_x),(y-surf_y),(z-surf_z)};
    
    const double tol=0.5e-15;
    for (int i=0;i<3;i++){
        if(fabs(diffs[i])<tol){diffs[i]=0;}
        dist.push_back(sqrt(diffs[i]*diffs[i]));
        
    }
    double length=length_line_square(theta*(180/PI),phi*(180/PI),face);
    double radial_dist=fabs(rad_dist-get_radius(tt));
    
    if ((length/2.0)<radial_dist){return 1.0;} // the actual boundary lies in the other cell.
    else{//std::cout<<((length/2.0)-radial_dist)/length<<"\n";
        return(((length/2.0)-radial_dist)/length); }
    
    
    
    
    
}





//********************************************************************
//**************                              ************************
//**************      length_line_square      ************************
//**************                              ************************
//********************************************************************
// Given the radial angles for a line and the face of the cube it passes through this caculates the length of the line it is needed in the volume fraction caluclation for the conservation of mass check. (NB phi is the azmuith (xy angle) theta is the incline (z axis)

double Bubble::length_line_square(double theta, double phi,int face){
    double    PI=4*atan(1);
    // double theta_lim=90-(acos(sqrt(2)/sqrt(3))*(180/PI));
    double incline=(90-theta)*(PI/180);
    double xy_ang= phi*(PI/180);
    // std::cout<<incline<<" incline\n";
    // std::cout<<xy_ang<<" xy ang\n";
    
    if (face==0 or face==3){//std::cout<<"x=1\n";
        double length=1/(cos(xy_ang)*cos(incline));
        //        std::cout<<"xy diagonal= "<<1/cos(xy_ang)<<"\n";
        //        std::cout<<"y = "<< length*(cos(incline)*sin(xy_ang))<<"\n";
        //        std::cout<<"z = "<< length*(sin(incline))<<"\n";
        return fabs(length);}
    
    else if(face==1 or face==4){//std::cout<<"y=1\n";
        double length=1/(sin(xy_ang)*cos(incline));
        //        std::cout<<"xy diagonal= "<<1/sin(xy_ang)<<"\n";
        //        std::cout<<"x = "<< length*(cos(incline)*cos(xy_ang))<<"\n";
        //        std::cout<<"z = "<< length*(sin(incline))<<"\n";
        return fabs(length);}
    
    else if(face==2 or face==5){//std::cout<<"z=1\n";
        double length=1/sin(incline);
        //        std::cout<<"xy diagonal= "<< length*(cos(incline))<<"\n";
        //        std::cout<<"y length = "<<length*(cos(incline))*sin(xy_ang)<<"\n";
        //        std::cout<<"x length = "<<length*(cos(incline))*cos(xy_ang)<<"\n";
        return fabs(length);}
    
    
    else {std::cout<<"angle not handelled\n"; return 5.0;}
}



//********************************************************************
//**************                              ************************
//**************      dist to point_x           ************************
//**************                              ************************
//********************************************************************
// function that works out the distance from a bulk point to the edge of the bubble. Done by the intercept of line x=... or y=... and the bubble equation.
double Bubble::xdist(int i, int j, int k){
    double boundary_dist;
    std::vector<double> surface_x_pt;
    surface_x_pt= LineSphereIntersections({static_cast<double>(i-1),static_cast<double>(j),static_cast<double> (k)}, {static_cast<double>(i+1),static_cast<double>(j),static_cast<double>(k)});
    if(surface_x_pt.size()<3){boundary_dist=1; return boundary_dist;}
    boundary_dist= i-surface_x_pt[0];
    //  std::cout<<"xdist= "<<boundary_dist<<"\n";
    return boundary_dist;
}
//********************************************************************
//**************                              ************************
//**************      dist to point_y           ************************
//**************                              ************************
//********************************************************************
// function that works out the distance from a bulk point to the edge of the bubble. Done by the intercept of line x=... or y=... and the bubble equation.
double Bubble::ydist(int i, int j, int k){
    double boundary_dist;
    std::vector<double> surface_y_pt;
    surface_y_pt= LineSphereIntersections({static_cast<double>(i),static_cast<double>(j-1),static_cast<double> (k)}, {static_cast<double>(i),static_cast<double>(j+1),static_cast<double>(k)});
    if(surface_y_pt.size()<3){boundary_dist=1; return boundary_dist;}
    
    boundary_dist= j-surface_y_pt[1];
    //std::cout<<"ydist= "<<boundary_dist<<"\n";
    return boundary_dist;}
//********************************************************************
//**************                              ************************
//**************      dist to point_z         ************************
//**************                              ************************
//********************************************************************
// function that works out the distance from a bulk point to the edge of the bubble. Done by the intercept of line x=... or y=... and the bubble equation.
double Bubble::zdist(int i, int j, int k){
    double boundary_dist;
    std::vector<double> surface_z_pt;
    surface_z_pt= LineSphereIntersections({static_cast<double>(i),static_cast<double>(j),static_cast<double> (k-1)}, {static_cast<double>(i),static_cast<double>(j),static_cast<double>(k+1)});
    if(surface_z_pt.size()<3){boundary_dist=1; return boundary_dist;}
    
    boundary_dist=k-surface_z_pt[2];
    // std::cout<<"zdist= "<<boundary_dist<<"\n";
    return boundary_dist;}

//********************************************************************
//**************                              ************************
//**************      which face              ************************
//**************                              ************************
//********************************************************************
//Given a node point (i,j,k) it works out which face of the node cube a line from the center of the bubble would cross. It is needed for the volume fraction calculation.
int Bubble::which_face(double i, double j, double k){
    double    PI=4*atan(1);
    double x=i;
    double y=j;
    double z=k;
    std::vector<double> ang_dists;
    double mag_vec=sqrt(x*x +y*y + z*z);
    ang_dists.push_back(acos(x/mag_vec));
    ang_dists.push_back(acos(y/mag_vec));
    ang_dists.push_back(acos(z/mag_vec));
    ang_dists.push_back(acos(-x/mag_vec));
    ang_dists.push_back(acos(-y/mag_vec));
    ang_dists.push_back(acos(-z/mag_vec));
    
    //        std::cout<<ang_dists[0]*(180/PI)<<","<<ang_dists[1]*(180/PI)<<","<<ang_dists[2]*(180/PI)<<","<<ang_dists[3]*(180/PI)<<","<<ang_dists[4]*(180/PI)<<","<<ang_dists[5]*(180/PI)<<"\n";
    std::vector<double>::iterator result= std::min_element(std::begin(ang_dists), std::end(ang_dists));
    // std::cout << "min element at: " << std::distance(std::begin(ang_dists), result)<<"\n";
    return std::distance(std::begin(ang_dists), result);
}


//********************************************************************
//**************                              ************************
//**************      dc/dr                   ************************
//**************                              ************************
//********************************************************************
//// The second order finite difference that I actually understand
//double Bubble::dcdr_2nd_order()
//dcdr_points
//

//********************************************************************
//**************                              ************************
//**************     spherical transform      ************************
//**************                              ************************
//********************************************************************
/*Part of the function for calcualting dc/dr based on a linear interpolation rather than JP's original formulation. In this method a unique set of points in in spherical polar co-ords are transformed into cartesian co-ords so that the conc can be found from the tissue concentration grid*/
std::vector<std::vector<double> >Bubble::spherical_tansform_r1(double deltar,double dangle){
    std::vector<std::vector<double> > dcdr_points_r1;
    int dr=1;
    // std::cout<<dangle*(180/PI)<<" angle \n";
    for (double phi=0; phi<2.0*PI; phi+=dangle){
        for (double theta=0; theta<=PI; theta+=dangle){
            // int dr=1;
            // std::cout<<phi*(180/PI)<<" phi "<<theta*(180/PI)<<" theta "<<cos(0)<<"\n";
            std::vector<double>  row;
            std::vector<double> angular_result;
            angular_result.push_back(sin(theta)*cos(phi));
            angular_result.push_back(sin(theta)*sin(phi));
            angular_result.push_back(cos(theta));
            
            const double tol=1e-14;
            for (int t=0;t<3;t++){
                if(fabs( angular_result[t])<tol){ angular_result[t]=0.0;}
            }
            
            double x=angular_result[0]*(deltar+get_radius())+get_locations(0);
            double y=angular_result[1]*(deltar+get_radius())+get_locations(1);
            double z=angular_result[2]*(deltar+get_radius())+get_locations(2);
            //  std::cout<<x<<","<<y<<","<<z<<" spherical pts\n";
            row.push_back(x); //this is here to deal with imprecision through using PI
            row.push_back(y);
            row.push_back(z);
            
            dcdr_points_r1.push_back(row);
            
        }
    }
    std::sort(dcdr_points_r1.begin(),dcdr_points_r1.end()); //sort the points
    
    dcdr_points_r1.erase(std::unique(dcdr_points_r1.begin(),dcdr_points_r1.end()),dcdr_points_r1.end()); //removing duplicate points
    assert(dcdr_points_r1.size()==((2.0*PI/dangle)*((PI/dangle)-1)+2)); //checking there are the correct number of points for thetheta and phi specified
    //   for (int i=0;i<dcdr_points_r1.size();i++){
    //  std::cout<<dcdr_points_r1[i][0]<<","<<dcdr_points_r1[i][1]<<","<<dcdr_points_r1[i][2]<<" r1 points \n";
    // }
    return dcdr_points_r1;
}
//*****************************************************************************************************************************************//
//std::vector<std::vector<double> >Bubble::spherical_tansform_r2(double deltar,double dangle){
//    std::vector<std::vector<double> > dcdr_points_r2;
//
//   // int dr=2;
//    for (double phi=0; phi<2.0*PI; phi+=dangle){
//        for (double theta=0; theta<=PI; theta+=dangle){                // int dr=1;
//                std::vector<double>  row;
//           // std::cout<<phi*(180/PI)<<" phi "<<theta*(180/PI)<<" theta "<<cos(0)<<"\n";
//            std::vector<double> angular_result;
//            angular_result.push_back(sin(theta)*cos(phi));
//            angular_result.push_back(sin(theta)*sin(phi));
//            angular_result.push_back(cos(theta));
//
//            const double tol=1e-14;
//            for (int t=0;t<3;t++){
//                if(fabs( angular_result[t])<tol){ angular_result[t]=0.0;}
//            }
//
//            double x=angular_result[0]*((deltar*2.0)+get_radius())+get_locations(0);
//            double y=angular_result[1]*((deltar*2.0)+get_radius())+get_locations(1);
//            double z=angular_result[2]*((deltar*2.0)+get_radius())+get_locations(2);
//          //  std::cout<<x<<","<<y<<","<<z<<" spherical pts\n";
//            row.push_back(x); //this is here to deal with imprecision through using PI
//            row.push_back(y);
//            row.push_back(z);
//                dcdr_points_r2.push_back(row);
//            }
//        }
//
//
//    std::sort(dcdr_points_r2.begin(),dcdr_points_r2.end()); //sort the points
//    dcdr_points_r2.erase(std::unique(dcdr_points_r2.begin(),dcdr_points_r2.end()),dcdr_points_r2.end()); //removing duplicate points
//     assert(dcdr_points_r2.size()==((2.0*PI/dangle)*((PI/dangle)-1)+2));
//    return dcdr_points_r2;
//    }
//
//********************************************************************
//**************                              ************************
//**************     dc_dr_spherical          ************************
//**************                              ************************
//********************************************************************
//calculates dc/dr based on simple difference method (dc/dr=(c_out - c_int)/diffusion region thickness). In this case dc/dr for each point on the bubble surface is calulated and averaged to give the final dc/dr.

double Bubble::dcdr_sphere(std::vector<double>newconcvec_r1_N2,std::vector<double>newconcvec_r1_O2,std::vector<double>newconcvec_r2_N2,std::vector<double>newconcvec_r2_O2,double dr,int tt){
    assert(newconcvec_r2_O2.size()==newconcvec_r1_O2.size());
    double dcdrN2_tot=0;
    double dcdrO2_tot=0;
    //  std::cout<<"conc outside bub is "<<newconcvec_r1_N2[1]<<"\n";
    for (int i =0; i<newconcvec_r1_N2.size(); i++){
        dcdrN2_tot+=(1.0/(2.0*dr))*((-newconcvec_r1_N2[i])+(4.0*newconcvec_r2_N2[i])-(3.0*get_concN2()));
        dcdrO2_tot+=(1.0/(2.0*dr))*((-newconcvec_r1_O2[i])+(4.0*newconcvec_r2_O2[i])-(3.0*get_concO2()));}
    
    
    double denom=newconcvec_r1_N2.size();
    double dcdrN2=dcdrN2_tot/denom;
    double dcdrO2=dcdrO2_tot/denom;
    dc_dr.push_back(dcdrN2+dcdrO2);
    dc_drN2.push_back(dcdrN2);
    dc_drO2.push_back(dcdrO2);
    
    // std::cout<<dcdr<<"\n";
    // }
    return dc_dr.back();
}

//********************************************************************
//**************                              ************************
//**************      variable perm           ************************
//**************                              ************************
//********************************************************************
//this function is always called by dr/dt but not necessarily used
void Bubble::var_perm(){
    // double perm=diffusivity*(mu/2)*(5-4*(radius[0]/radius.back())-pow((radius[0]/radius.back()),4)); //diffusion coeff changes as a function of the non-dimensionalised gent pressure
    double perm=diffusivity2*pow(1.0-(min_radius/radius.back()),2);
    //std::cout<<"new diffusivity is "<<perm<<"\n";
    double tauvar=(k_h_N2*perm*RT_N2*RT_O2)/(pow(rbar,2)*((pp_fraction_N2*RT_O2)+pp_fraction_O2*RT_N2)); // non dim param with variable bubble perm;
    set_tauvar(tauvar);
}

//********************************************************************
//**************                              ************************
//**************      dr/dt  for gent         ************************
//**************                              ************************
//********************************************************************
////this is the calcualtion of dr/dt using the formular derived in model derivation.tex it needs to know the current time step in order to get the previous P_amb and it needs to know the size of the time step. In tissue.cpp each bubble will be fed in to this function in each case the new radius is pushed back ready for the next iteration.
void Bubble::drdt(int t,double timestep){
    var_perm();
    double rad_initial=initial_radius;////1e-6/rbar;min_radius;//
    double dpdt = (P_amb[t] - P_amb[t-1])/timestep;
    // double dcdr_term=tau*get_dcdr();
    //to use varying permeability un-comment below
    double dcdr_term=tauvar*get_dcdr();
    double surface_tension_term=4.0*sigma/(3.0*radius.back());
    double elasticity_term=(5.0*mu/2.0)-(4.0*mu/3.0)*(rad_initial/radius.back()) + (mu/6.0)*pow((rad_initial/radius.back()),4);
    
    // std::cout<<"dc/dr is "<<dc_dr.back()<<"\n";
    double radnow=(timestep*(dcdr_term - dpdt*(radius.back()/3.0))/(P_amb[t]+surface_tension_term+elasticity_term))+radius.back();
    
    
    /* uncomment to see how each component contributes the the bubble radius */
    //          std::cout<<"denominator of dr/dt "<<P_amb[t]+(4*sigma/(3*radius.back())) +(5*mu/2) -(4*mu/3)*(radius[0]/radius.back()) + (mu/6)*pow((radius[0]/radius.back()),4)<<"\n";
    //     std::cout<<"dpdt P_amb  "<<dpdt<<"\n";
    //        std::cout<<"new elasticity "<<(4*sigma/(3*radius.back())) +(5*mu/2) -(4*mu/3)*(radius[0]/radius.back()) + (mu/6)*pow((radius[0]/radius.back()),4)<<"\n";
    //    std::cout<<"tau is"<<tau<<"\n";
    //      std::cout<<"boyles law term = "<<-timestep*((radius.back()/3.0)*dpdt)<<"\n";
    //     std::cout<<"dc/dr term = "<<timestep*tauvar*get_dcdr()<<"\n";
    //        // std::cout<<
    //    std::cout<<"radius1 is "<<radnow<<"\n";
    
    
    
    //if(radnow<radius[0]){radnow=radius[0];}
    if (radnow<min_radius)
    {radnow=min_radius;}
    radius.push_back(radnow); //put the new radius into the variable
    mass_in_bubble(t, timestep); ///put the new bubble conc into the variable
    
}

//********************************************************************
//**************                              ************************
//**************      dr/dt gern              ************************
//**************                              ************************
//********************************************************************
////this is the calcualtion of dr/dt using Gernhardts formulation of elasticity. It needs to know the current time in order to get the previous P_amb and it needs to know the size of the timestep. In tissue each bubble will be fed in to this function in each case the new radius is pushed back ready for the next iteration.
void Bubble::drdt_gern(int t,double timestep){
    double dpdt = (P_amb[t] - P_amb[t-1])/timestep;
    // std::cout<<"dc/dr is "<<dc_dr.back()<<"\n";
    double radnow=timestep*((tau*get_dcdr() - dpdt*radius.back()/3.0)/
                            (P_amb[t]+(4.0*sigma/(3.0*radius.back())) + (8.0*PI/3.0)*M*pow(radius.back(),3)))+radius.back();
    
    //Looking at each component separately
    //     std::cout<<"denominator of dr/dt "<<P_amb[t]+(4*sigma/(3*radius.back())) +(5*mu/2) -(4*mu/3)*(radius[0]/radius.back()) + (mu/6)*pow((radius[0]/radius.back()),4)<<"\n";
    //     std::cout<<"P_amb  "<<P_amb[t]<<"\n";
    //     std::cout<<"new elasticity "<<(4*sigma/(3*radius.back())) +(5*mu/2) -(4*mu/3)*(radius[0]/radius.back()) + (mu/6)*pow((radius[0]/radius.back()),4)<<"\n";
    //
    //    std::cout<<"boyles law term = "<<-timestep*((radius.back()/3.0)*dpdt)<<"\n";
    // std::cout<<"dc/dr term = "<<timestep*tau*get_dcdr()<<"\n";
    // std::cout<<
    // std::cout<<"radius now is "<<radnow<<"\n";
    //if(radnow<radius[0]){radnow=radius[0];}
    if (radnow<min_radius)
    {radnow=min_radius;}
    radius.push_back(radnow); //put the new radius into the variable
    mass_in_bubble(t, timestep); ///put the new bubble conc into the variable
    
}

//********************************************************************
//**************                              ************************
//**************      dr/dt laplace           ************************
//**************                              ************************
//********************************************************************
void Bubble::drdt_laplace(int t,double timestep){
    double dpdt = (P_amb[t] - P_amb[t-1])/timestep;
    //std::cout<<"dc/dr is "<<dc_dr.back()<<"\n";
    double radnow=+radius.back()+ timestep*((tau*get_dcdr() - dpdt*radius.back()/3.0)/
                                            (P_amb[t]+(4.0*sigma/(3.0*radius.back()))));
    
    
    // std::cout<<"denominator of dr/dt "<<(P_amb[t]+(4*sigma/(3*radius.back())))<<"\n";
    //std::cout<<"P_amb  "<<P_amb[t]<<"\n";
    // std::cout<<"tau is"<<tau<<"\n";
    //std::cout<<"boyles law term = "<<-timestep*((radius.back()/3.0)*dpdt)<<"\n";
    //std::cout<<"dc/dr term = "<<timestep*tau*dc_dr.back()<<"\n";
    if (radnow<min_radius)
    {radnow=min_radius;}
    radius.push_back(radnow); //put the new radius into the variable
    mass_in_bubble(t, timestep); ///put the new bubble conc into the variable
    
}


//********************************************************************
//**************                              ************************
//**************    sphere_line intersection  ************************
//**************                              ************************
//********************************************************************
std::vector<double> Bubble::LineSphereIntersections( std::vector<double> linePoint0,  std::vector<double> linePoint1){
    
    double px = linePoint0[0];
    double py = linePoint0[1];
    double pz = linePoint0[2];
    //   std::cout<<px<<","<<py<<","<<pz<<" in line intersect pt1\n";
    // std::cout<<linePoint1[0]<<","<<linePoint1[1]<<","<<linePoint1[2]<<" in line intersect pt2\n";
    double vx = linePoint1[0] - px;
    double vy = linePoint1[1] - py;
    double vz = linePoint1[2] - pz;
    
    
    double A = vx * vx + vy * vy + vz * vz;
    double B = 2.0 * ((px * vx) + (py * vy) + (pz * vz) - (vx * locations[0]) - (vy * locations[1]) - (vz * locations[2]));
    double C = px * px - 2.0 * px * locations[0] + locations[0] * locations[0] + py * py - 2.0 * py * locations[1] + locations[1] * locations[1] + pz * pz - 2.0 * pz * locations[2] + locations[2] * locations[2] - radius.back() * radius.back();
    
    
    // discriminant
    double D = B * B - 4 * A * C;
    // std::cout<<D<<" discriminant\n";
    
    if ( D < 0 )
    {
        return {NAN};
    }
    
    double t1 = ( -B - sqrt ( D ) ) / ( 2.0 * A );
    // std::cout<<t1<<" \n";
    std::vector<double> solution1;
    solution1= {( linePoint0[0] * ( 1 - t1 ) + t1 * linePoint1[0]),(linePoint0[1] *( 1 - t1 ) + t1 * linePoint1[1]),(linePoint0[2] * ( 1 - t1 ) + t1 * linePoint1[2] )};
    
    
    if ( D == 0 )
    {
        return {NAN}; //this is because we are not interested in tangents (they will have no volume in the cube of interest)
    }
    
    double t2 = ( -B + sqrt( D ) ) / ( 2.0 * A );
    //std::cout<<t2<<" \n";
    std::vector<double> solution2;
    solution2= {( linePoint0[0] * ( 1 - t2 ) + t2 * linePoint1[0]),
        (linePoint0[1] * ( 1 - t2 ) + t2 * linePoint1[1]),
        (linePoint0[2] * ( 1 - t2 ) + t2 * linePoint1[2] )};
    
    
    
    // prefer a solution that's on the line segment itself
    if ( abs( t1 - 0.5 ) < abs( t2 - 0.5 ) )
    {
        std::vector<double> solutions12;
        solutions12=solution1;
        //solution1.insert(std::end(solution1), std::begin(solution2), std::end(solution2));
        //std::cout<<solution1.insert(std::end(solution1), std::begin(solution2), std::end(solution2))[0]<<" test \n";
        //   std::cout<<solution1[0]<<solution1[1]<<solution1[2]<<"\n";
        
        return solutions12;
        
    }
    
    
    std::vector<double> solutions21;
    solutions21=solution2;
    // std::cout<<solutions21[0]<<"\n";
    //std::cout<<solution2.insert(std::end(solution2), std::begin(solution1), std::end(solution1))[0]<<" test \n";
    //  std::cout<<solution2[0]<<solution2[1]<<solution2[2]<<"\n";
    return solutions21;
    
}
