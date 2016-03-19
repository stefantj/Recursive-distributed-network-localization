/*
 *  system.h
 *
 *  Created by Stefan Jorgensen for Stanford EE373A project (Winter 2016).
 *  based on earlier work for EE359.
 *
 *  Code released under the MIT License (MIT)
 *
 */

#ifndef _PARAMETERS_H__
#define _PARAMETERS_H__

//Default parameter values: 
#define NETWORK_SIZE 20
#define ARENA_SIZE_DEFAULT 30

//Standard deviation index values
#define SDEV_R 0  //Variance of R meaeasurements
#define SDEV_U 1  //Variance of U_x, U_y measurements (= 1/2 * variance of theta)
#define SDEV_G 2  //Variance of G_x, G_y measurements

//Adaption parameters
//Step Size:
#define MU_PARAM 1.0

//Simulation type
#define BERNOULLI 0
#define DEPENDENT_BERNOULLI 1
#define SINR_FAIL 2

#include "coord.h"

//Simulator prototype:
class Simulator;

//Node prototype:
class Node;

#endif //_PARAMETERS_H__
