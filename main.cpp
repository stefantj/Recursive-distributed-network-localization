/*  main.cpp
 *  EE373A Project "Adaptive Recursive Distributed Filtering for Network Localization"
 *  Written by Stefan Jorgensen
 *  Stanford University, Winter 2016
 *  This work is funded in part by NSF grant DGE-114747
 *
 *  PROGRAM DESCRIPTION:
 *   Heavily based on earlier project for EE359.
 *   Uses distributed filtering methods developed by Professor Sayed at UCLA, but with output feedback.
 *
 *  LICENCE:
 *  This program is released under GPL 3.0:
 *  http://www.gnu.org/copyleft/gpl.html
 */

/* General libraries */
#include <thread>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

/* GSL helpers */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort_vector.h>

/* Simulator parameters */
#include "Globals.h"

/* Program header files */
#include "coord.h"
#include "Node.h"
#include "simulator.h"



//Node list:
Node* NODES[NETWORK_SIZE];

//Neighborhood container
//first index: k,
//second index:
// first element: |N_k|
// next |N_k|-1 elements: indices of neighbors
int NEIGHBORS[NETWORK_SIZE][NETWORK_SIZE+1];

//Global simulator object (needs to be global for nodes?
Simulator Sim;

int main (int argc, char * const argv[]) {
	
	//run simulations!
	float sdevs[3];
	sdevs[SDEV_R] = 0.5;
	sdevs[SDEV_U] = 2500;
	sdevs[SDEV_G] = 10;
	char prefix[3] = {'B','0','0'};
	Sim.setLoops(1, 1);
	//Mobile Trials:
	
	//Run with beta = 1:
	prefix[0] = 'M';
	Sim.tracking_trial(prefix, 0.5, sdevs);

	//Run with beta = 0.3
//	prefix[2] = '1';
//	Sim.tracking_trial(prefix, 0.3, sdevs);

	//Run with beta = 0
//	prefix[2] = '2';
//	Sim.tracking_trial(prefix, 0, sdevs);
	
	return 0;
}

