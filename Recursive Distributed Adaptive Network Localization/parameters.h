/*
 *  system.h
 *  Simulator
 *
 *  Created by Stefan Jorgensen on 10/8/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _PARAMETERS_H__
#define _PARAMETERS_H__


//this should really be its own class...

//Default parameter values: 
#define NETWORK_SIZE 50
#define ARENA_SIZE_DEFAULT 300

//Standard deviation index values
#define SDEV_R 0 //= 0.7; //Variance of R meaeasurements
#define SDEV_U 1 //= 0.025; //Variance of U_x, U_y measurements (= 1/2 * variance of theta)
#define SDEV_G 2 //= 10;   //Variance of G_x, G_y measurements

//Adaption parameters
//Step Size:
#define MU_PARAM 1.7

//Simulation type
#define BERNOULLI 0
#define DEPENDENT_BERNOULLI 1
#define SINR_FAIL 2

#include "coord.h"

//Simulator prototype:
class Simulator;

//Node prototype:
class Node;

//get_meas(j,k, r, u)
//
// if j=k, returns noisy g measurement in place of u
// else returns noisy r, u measurements
//bool get_meas(int j, int k, float* r, coord u);

#endif //_PARAMETERS_H__
