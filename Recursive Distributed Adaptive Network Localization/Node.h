/*
 *  Node.h
 *
 *  Created by Stefan Jorgensen for Stanford EE373A project (Winter 2016).
 *  based on earlier work for EE359.
 *
 *  Code released under the MIT License (MIT)
 *
 */

#ifndef _NODE_H__
#define _NODE_H__

#include "Globals.h" //contains system parameters
#include "coord.h" //contains coord class

#define USE_RECURSIVE_LMS

#ifdef USE_RECURSIVE_LMS
#define RECURSIVE_HIST 15
#endif
//extern Node* NODES[NETWORK_SIZE]; //global network list
class Simulator;

class Node
{
public: 
	//Constructor
	// k: index of the node.
	// sets position w_0 to fit within system boundaries
	Node(int n_index);
	
	//initialize()
	// fills in initial data
	//
	// for all j
	//  if(j in N_k)
	//   w_jk[j][X] = NODES[j]->get_w_jk(j, X);
	//   w_jk[j][Y] = NODES[j]->get_w_jk(j, Y);
	//   mu_j[j] = MU_PARAM; //default u value from parameters.h
	//  else
	//   w_jk[j] = {0,0};
	//   mu_j[j] = 1;
	//
	void initialize();

	//ATC functions: 
	
	//measure()
	// for j in Nk\{k}
	//   float r_j
	//   coord u_j
	//   if(get_meas(j,k,&r_j, &u_j))
	//     d_jk[j] = r_j - coord_mult(u_j, (coord)w_jk[k]);
	//     u_jk[j][X] = u_j[X];
	//	   u_jk[j][Y] = u_j[Y];
	//     M_k[j] = true;
	//   else 
	//     M_k[j] = false;
	// M_k[k] = true;
	// coord g_j;
	// get_meas(k,k, NULL, g_j);
	// u_j = coord_diff(g_j, w_jk[k]);
	// r_j = coord_norm(u_j);
	// coord_scale(1/r_j, u_j);
	// u_jk[j][X] = u_j[X];
	// u_jk[j][Y] = u_j[Y];
	// d_kk[j] = r_j-coord_mult(u_j, w_jk[k]);
	// 
	//Returns number of nodes in measurement space (for analysis)
	int measure(); 
	
	//calc_psi()
	// calculates preliminary psi
	// for all j 
	//   float c_j = 0;
	//   float d, tmp;
	//   coord u;
	//   coord psi = {0,0};
	//   for l in N_k
	//     if(NODES[l]->get_du(j, &d, u))
	//		c_j++;
	//      tmp = coord_mult(u, w_jk[j]);
	//      tmp = d - tmp;
	//      coord_scale(tmp, u);
	//      psi = coord_add(psi, u, psi);
	//   if(c_j != 0)
	//     NM_k[j] = true;
	//     coord_scale(mu_j[j]/c_j, psi);
	//     coord_add(w_jk[j], psi, w_jk[j]);
	//     p_jk[j][X] = psi[X];
	//     p_jk[j][Y] = psi[Y];
	//   else
	//     NM_k[j] = false;
	//returns number of measurements used to calc psi_kk
	int calc_psi();

	//calc_w()
	// calculates estimates w_jk
	// for all j
	//   float a_j = 0;
	//   coord p;
	//   coord w = {0,0};
	//   for l in N_k
	//     if(NODES[l]->get_psi(l,p))
	//       a_j++;
	//       coord_add(w, p, w);
	//   coord_scale(1/a_j, w);
	//   w_jk[j][X] = w[X];
	//   w_jk[j][Y] = w[Y];
	void calc_w();
	
	//calc unaided estimate (normal localization)
	
	//Data getters
	
	//check neighborhood for node j
	bool check_neighborhood(int j); 
	
	//measurement data getter:
	//returns false if l not in M_k
	bool get_du(int l, float* d, coord u);	
	
	//intermediate estimate data getter:
	//returns false if j not in NM_k
	bool get_psi(int l, coord p);
    
#ifdef USE_RECURSIVE_LMS
    //feedback intermediate estimate data getter:
    //returns false if j not in NM_k
    bool get_psif(int l, float* pf);
#endif
	
	//get_w_jk(j,dim)
	// return w_jk[j][dim];
	float get_w_jk(int j, int dim); 
	
	//get_w_0(dim)
	//return w_0[dim]
	float get_w_0(int dim);
	
	//get noncooperative values
	float get_nc_w_jk(int j, int dim);
	
	//For power adaptation:
	float get_power(void);
	void  set_power(float pwr);
	
    //MSE estimate
    float get_error_est(void);
    
	//For motion:
	//Given arena boundaries, update position.
	void update_position(float arena_x, float arena_y);
	void randomize_position();
	void set_velocity(float vx, float vy);
	
private:
	bool M_k[NETWORK_SIZE];  // if M_k[j] == true, then j is in M_k
	bool NM_k[NETWORK_SIZE]; // if NM_k[j] == true, then j is in U M_l for l in N_k 
	//Estimators:
	float d_jk[NETWORK_SIZE];    //range
	float u_jk[NETWORK_SIZE][2]; //bearing
	float p_jk[NETWORK_SIZE][2]; //intermediate weight
	float w_jk[NETWORK_SIZE][2]; //weight (position estimate)
	float nc_w_jk[NETWORK_SIZE][2];
	float mu_j[NETWORK_SIZE];
    float e_jk[NETWORK_SIZE]; //errors

    //Recursive estimators
#ifdef USE_RECURSIVE_LMS
    float ef_jk[NETWORK_SIZE];
	float y_jk[NETWORK_SIZE][NETWORK_SIZE][RECURSIVE_HIST];
	float wf_jk[NETWORK_SIZE][RECURSIVE_HIST];
    float pf_jk[NETWORK_SIZE][RECURSIVE_HIST];

#endif
	//true position
	coord w_0;
	//velocity (units m/s)
	coord vel;
	int self_index;
	float power;
};

#endif //_NODE_H__