/*
 *  simulator.h
 *
 *  Created by Stefan Jorgensen for Stanford EE373A project (Winter 2016).
 *  based on earlier work for EE359.
 *
 *  Code released under the MIT License (MIT)
 *
 */

#ifndef _SIMULATOR_H__
#define _SIMULATOR_H__

#include "Globals.h"
#include "coord.h"
#include <iostream>
#include <fstream>
#include <math.h>

/* GSL random number libraries */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std;
class Node;

#define NUM_FILES		5
#define W_0_FILE		0
#define W_K_FILE		1
#define MSE_FILE		2
#define INST_MSE		3
#define PARAM_FILE		4

#define ITER_MAX 1000

class Simulator
{
public:
	
	//sim_init() replaced by constructor
	//initializes the simulation
	// set up file pointers
	// create nodes 
	Simulator();
	
	//Runs fixed position trial (returns final MSE)
	float fixed_pos_trial(char* prefix, float param, float* stdevs);
	void tracking_trial(char* prefix, float param, float* stdevs);
	//Randomizes node positions
	void newPositions();
	//resets mse container
	void reset_mse(); 
	void setArenaWidth(int width);
	int getArenaWidth();
	void setArenaHeight(int height);
	int getArenaHeight();
	int getNradius();
	void setMradius(int mrad);
	int getMradius();
	float getMbar();
	void setLoops(float loops, float iloops);
	bool get_meas(int j, int k, float* r, coord u);
	
private:

	//reset node information
	void reset_nodes();
	
	//run simulation
	// initialize nodes
	// for n = 1:ITER_MAX
	//  update_neighbors();
	//  for all j
	//   NODES[j]->initialize();
	//  for all j
	//   NODES[j]->measure();
	//  for all j
	//   NODES[j]->calc_psi();
	//  for all j
	//   NODES[j]->calc_w();
	//
	//  save_frame();
	//
	void run_sim(int simtype, const char* prefix);
	
	//Computes network MSE measurements.
	// mse = 0
	// for all j
	//   mse += |w_jj - w_j_0|^2
	// return mse/NETWORK_SIZE
	float get_NETMSE();

	//Saves position information for the network
	//kinda messy, but we'll use a global filename container.
	void save_scene(const char* prefix, int file_index);
    void save_julia_data(bool close);
	
//MEASUREMENT HELPERS
	//updates neighborhood 
	// varying degrees of complexity according to FAIL_MODE
	// for all j
	//   NEIGHBORS[j][0] = 1; //default with self in neighborhood
	//	 NEIGHBORS[j][1] = j; //include self in neighborhood
	//   for l > j //check everyone else:
	//     if l,j connected //<- secret sauce is here...
	//       NEIGHBORS[j][NEIGHBORS[j][0]] = l;
	//       NEIGHBORS[j][0]++;
	//       NEIGHBORS[l][NEIGHBORS[l][0]] = j;
	//       NEIGHBORS[l][0]++;
	//
	void update_neighbors();

	//calculate the distance between nodes j, l
	float dist_nodes(int j, int l);

//VARIABLES:
    //Random number generator:
    gsl_rng* RNG;
	//Network MSE container:
	float NET_MSE[ITER_MAX];	
	//Arena parameters:
	int ARENA_WIDTH;
	int ARENA_HEIGHT;
	int M_RADIUS;
	int N_RADIUS;
	//Standard deviation values
	// SDEV_R 0
	// SDEV_U 1
	// SDEV_G 2
	float STDEVS[3];	
	//Communications mode:
	int FAIL_MODE;	
	//Number of realizations to test
	float LOOPS;
	//Number of times to test each realization
	float ILOOPS; 
	//Bernoulli failure probability
	float BETA;
	//Parameters for path loss modelling
	int NUM_TIMESLOTS;
	float SINR_THRESH; //units mW
	
	//average number of neighbors
	float Nbar_tot; 
	int Nbar_cnt;

	float Mbar;
	
	//filenames:
	char filenames[NUM_FILES][12];
	FILE* files[NUM_FILES];
    bool firstWrite;
};

#endif //_SIMULATOR_H__