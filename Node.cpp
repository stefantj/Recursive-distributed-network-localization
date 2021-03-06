/*
 *  Node.cpp
 *
 *  Created by Stefan Jorgensen for Stanford EE373A project (Winter 2016).
 *  based on earlier work for EE359.
 *
 *  Code released under the MIT License (MIT)
 *
 */

#include "Node.h"
#include "simulator.h"
#include <math.h>
#include <iostream>
using namespace std;
class Simulator;
extern Simulator Sim; //include global simulator object
//Node list:
extern Node* NODES[NETWORK_SIZE];
extern int NEIGHBORS[NETWORK_SIZE][NETWORK_SIZE+1]; //global neighbor list

Node::Node(int n_index) : self_index(n_index)
{	
	float pos = 0; 
	//resolution 0.1m
	pos = rand() % (Sim.getArenaWidth()*10);
	w_0[X] = pos/10;
	pos = rand() % (Sim.getArenaHeight()*10);
	w_0[Y] = pos/10;

    //Simple test: Start near center and pointing at each other.
//    w_0[X] = Sim.getArenaWidth()/2+(1-2*n_index);
//    w_0[Y] = Sim.getArenaHeight()/2+5*(1-2*n_index);
    
//    vel[X] = 0;
//    vel[Y] = (1-2*n_index);
    
	//assign random velocity in +-[5, 15] m/s
//	float sx = 1 - 2*(rand()%2);
//	float sy = 1 - 2*(rand()%2);
//	vel[X] = sx*(((float)rand())/(INT_MAX))*2+1;
//	vel[Y] = sy*(((float)rand())/(INT_MAX))*2+1;
    vel[X] = (float)rand()/INT_MAX-0.5;
    vel[Y] = (float)rand()/INT_MAX-0.5;
	
	printf("Node %d created @ (%f, %f) v = <%f, %f>\n", self_index, w_0[X], w_0[Y], vel[X], vel[Y]);
	
	power = 100; //power in mW
	//initialize all arrays to a sensible value
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		M_k[j] = false;
		NM_k[j] = false;
		d_jk[j] = 0;
		mu_j[j] = MU_PARAM;
		u_jk[j][X] = 0;
		u_jk[j][Y] = 0; 
		p_jk[j][X] = 0;
		p_jk[j][Y] = 0;
		w_jk[j][X] = 0;
		w_jk[j][Y] = 0;
        e_jk[j] = 0;
#ifdef USE_RECURSIVE_LMS
        ef_jk[j] = 0;
        for(int jj = 0; jj < RECURSIVE_HIST; jj++){
            wf_jk[j][jj] = 0;
            pf_jk[j][jj] = 0;
            for(int l = 0; l < NETWORK_SIZE; l++)
                y_jk[j][l][jj]  = 0;
        }
#endif
	}
}

void Node::initialize(void)
{
	//initialize all arrays to a sensible value
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		M_k[j] = false;
		NM_k[j] = false;
		d_jk[j] = 0;
		mu_j[j] = MU_PARAM;
		u_jk[j][X] = 0;
		u_jk[j][Y] = 0;
		p_jk[j][X] = 0;
		p_jk[j][Y] = 0;
		w_jk[j][X] = 0;
		w_jk[j][Y] = 0;
        e_jk[j] = 0;
#ifdef USE_RECURSIVE_LMS
        ef_jk[j] = 0;
        for(int jj = 0; jj < RECURSIVE_HIST; jj++){
            wf_jk[j][jj] = 0;
            pf_jk[j][jj] = 0;
            for(int l = 0; l < NETWORK_SIZE; l++)
                y_jk[j][l][jj]  = 0;
        }
#endif
    }

	if(!Sim.get_meas(self_index, self_index, NULL, w_jk[self_index]))
		printf("Error measuring self (NODE %d) \n", self_index);
	
}

int Node::measure()
{
	int M = 0;
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		float r_j;
		coord u_j,g_j;
		if(!Sim.get_meas(self_index, self_index, NULL, g_j))
			printf("Error measuring self (NODE %d) \n", self_index);

		if(j!=self_index)
		{
			if(Sim.get_meas(j,self_index,&r_j,u_j))
			{
				d_jk[j] = r_j + coord_mult(u_j, w_jk[self_index]);
				u_jk[j][X] = u_j[X];
				u_jk[j][Y] = u_j[Y];
				M_k[j] = true;
				M++;
				//compute non-cooperative values
				nc_w_jk[j][X] = g_j[X]+r_j*u_j[X];
				nc_w_jk[j][Y] = g_j[Y]+r_j*u_j[Y];
			}else {
				M_k[j] = false;
			}
		}else{			
			M_k[j] = true;
			M++;
			coord_diff(g_j, w_jk[j], g_j);
			r_j = sqrt(coord_mult(g_j, g_j));
			coord_scale(1/r_j, g_j);
			u_jk[j][X] = g_j[X];
			u_jk[j][Y] = g_j[Y];
			d_jk[j] = r_j + coord_mult(g_j, w_jk[j]);
			//Store non-cooperative values
			nc_w_jk[j][X] = g_j[X];
			nc_w_jk[j][Y] = g_j[Y];
		}
	}
	return M;
}

//returns number of measurements used to compute W_k
int Node::calc_psi()
{
	int M = 0;
    // Update estimator for each node in the network:
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		int c_j =0;
		float d = 0;
		coord u = {0,0};
		coord psi = {0,0};
		int N_s = NEIGHBORS[self_index][0];
        
#ifdef USE_RECURSIVE_LMS
        float psi_f[RECURSIVE_HIST] = {};
        
        // y_jk[l][j] is the vector of histories of node l's measurement of node j.
        //shift history - we always do this whether or not there are measurements available.
        for(int l = 1; l <= NETWORK_SIZE; l++){
            for(int jj = RECURSIVE_HIST-1; jj >0; jj--)
                y_jk[l][j][jj] = y_jk[l][j][jj-1];
            // Fill in current value as zero in case no measurements.
            y_jk[l][j][0] = 0;
        }
#endif
        
        //Query Neighbors for measurements:
		//if you can talk to them
		for(int l = 1; l <= N_s; l++)
		{
			int ll = NEIGHBORS[self_index][l];
			//if they have measurements of node j
			if(NODES[ll]->get_du(j, &d, u))
			{
				c_j++;
                
#ifdef USE_RECURSIVE_LMS
                
                // Adaptive recursive case: y = a'x + b'y (x = u, a = w, b = w_f)
                // error_l(j) = d_l(j) - u_l(j)'w(j) - y_l(j)'w_f(j)
                
                //compute current y
                y_jk[ll][j][0] = coord_mult(u, w_jk[j]);
                if(isnan(y_jk[ll][j][0]))
                    printf("Neighbor's data is nan");
                
                for(int jj = 0; jj < RECURSIVE_HIST-1; jj++)
                    y_jk[ll][j][0] += (y_jk[ll][j][jj+1])*(wf_jk[j][jj+1]);
                if(isnan(y_jk[ll][j][0]))
                    printf("y-estimate went nan");
                
                // Compute Y'W_f for error term
                float ysum = 0;
                for(int jj = 0; jj < RECURSIVE_HIST-1; jj++)
                    ysum += wf_jk[j][jj]*y_jk[ll][j][jj];
                if(isnan(ysum))
                    printf("ysum went nan");
                
                //Compute error
                 ef_jk[j] = d - ysum - coord_mult(u, w_jk[j]);
                
                //Compute intermediate feed-back weight vector
                for(int jj=0; jj < RECURSIVE_HIST; jj++){
                    psi_f[jj] += ef_jk[j]*y_jk[ll][j][jj];
                    if(isnan(psi_f[jj]))
                        printf("psi_f went nan");
                }
                
                //Compute intermediate feed-forward weight vector
				coord_scale(ef_jk[j], u);
				coord_add(psi, u, psi);

                // Likely need some sanity checking here...
#else
                // Sum for standard case:
                float e_i = d - coord_mult(u,w_jk[j]);
				coord_scale(e_i, u);
				coord_add(psi, u, psi);
#endif
			}
		}
		
        //Scale and add offset
		if(c_j)
		{
            // Optional: normalize step size? Might help with stability
            mu_j[j] = .5;
			if(j == self_index)
				M = c_j;
			NM_k[j] = true;
			coord_scale(mu_j[j]/(float)c_j, psi);
			coord_add(w_jk[j], psi, psi);
			p_jk[j][X] = psi[X];
			p_jk[j][Y] = psi[Y];
            
#ifdef USE_RECURSIVE_LMS
            for(int jj = 0; jj<RECURSIVE_HIST; jj++)
                pf_jk[j][jj] =  wf_jk[j][jj] + 0.1*(mu_j[j]/(float)c_j)*psi_f[jj];
#endif
		} else{
			NM_k[j] = false;
		}
	}
	return M;
}

void Node::calc_w()
{

	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		float a_j = 0;
		coord p; 
		coord w = {0,0};
#ifdef USE_RECURSIVE_LMS
        float pf[RECURSIVE_HIST]={};
        float wf[RECURSIVE_HIST]={};
#endif
		int N_s = NEIGHBORS[self_index][0];
		for(int l = 1; l <= N_s; l++)
		{
			int ll = NEIGHBORS[self_index][l];
			if(NODES[ll]->get_psi(j,p))
			{
				a_j++;
				coord_add(w, p, w);
			}
#ifdef USE_RECURSIVE_LMS
            if(NODES[ll]->get_psif(j,pf))
            {
                for(int jj= 0; jj < RECURSIVE_HIST; jj++)
                    wf[jj]+=pf[jj];
            }
            
#endif
		}
		if(a_j)
		{
			coord_scale(1/a_j, w);
			w_jk[j][X] = w[X];
			w_jk[j][Y] = w[Y];
#ifdef USE_RECURSIVE_LMS
            //Normalize wf for stabilitiy?
            float wsum = 0;
            for(int jj=0; jj < RECURSIVE_HIST; jj++)
                wsum+=wf[jj];
            
            for(int jj = 0; jj < RECURSIVE_HIST; jj++){
                if(isnan(wf_jk[j][jj])){
                   printf("Recursive filter for node %d diverged.\n",self_index);
                    exit(1);
                }else{
                    wf_jk[j][jj] = wf[jj]/a_j;
                    if(wsum > 1 || wsum < -1)
                        wf_jk[j][jj] /= wsum;
                }
            }
#endif
            
		}else{
            printf("No data!\n");
        }
	}
    
}

bool Node::check_neighborhood(int j)
{
	for(int jj = 1; jj <=NEIGHBORS[self_index][0]; jj++)
	{
		if(NEIGHBORS[self_index][jj] == j) return true;
	}
	return false;
}


bool Node::get_du(int l, float* d, coord u)
{
	//check that l in M_k:
	if(M_k[l])
	{
		d[0] = d_jk[l];
		u[X] = u_jk[l][X];
		u[Y] = u_jk[l][Y];
		return true;
	}else{return false;}
}

bool Node::get_psi(int l, coord p)
{
	if(NM_k[l])
	{
		p[X] = p_jk[l][X];
		p[Y] = p_jk[l][Y];
		return true;
	}else{return false;}
}


#ifdef USE_RECURSIVE_LMS
bool Node::get_psif(int l, float* pf)
{
	if(NM_k[l])
	{
        for(int jj = 0; jj < RECURSIVE_HIST; jj++)
            pf[jj] = pf_jk[l][jj];
		return true;
	}else{return false;}
}
#endif

float Node::get_w_jk(int j, int dim)
{
	return w_jk[j][dim];
}

float Node::get_nc_w_jk(int j, int dim)
{
	return nc_w_jk[j][dim];
}

float Node::get_w_0(int dim)
{
	return w_0[dim];
}

float Node::get_power(void)
{
	return power;
}



void Node::set_power(float pwr)
{
	power = pwr;
}

void Node::randomize_position()
{
	float pos = 0; 
	//resolution 0.1m
	pos = rand() % (Sim.getArenaWidth()*10);
	w_0[X] = pos/10;
	pos = rand() % (Sim.getArenaHeight()*10);
	w_0[Y] = pos/10;	
}

float Node::get_error_est(){
    float error = 0;
    int N_s = NEIGHBORS[self_index][0];
    for(int l = 1; l <= N_s; l++)
        error+=e_jk[NEIGHBORS[self_index][l]]*e_jk[NEIGHBORS[self_index][l]];
    return error/N_s;
}

void Node::update_position(float arena_x, float arena_y)
{

    
    /* Potential fields model - repulsed by neighbors */
    float dt = 0.1; //units seconds
    float mag = 0;
    float mass = 10; // "mass" in f = m*a
    coord force = {}; //"force" from potential field
    coord attraction_force = {};

    
    for(int l = 1; l < NETWORK_SIZE; l++)
    {
        if(l != self_index){
            coord dir = {0,0};
            Sim.get_meas(l, self_index, &mag, dir);
            if(mag < 5 && mag > 0){
                force[X] += (dir[X]/mag)/(mag*mag);
                force[Y] += (dir[Y]/mag)/(mag*mag);
            }
        }
    }
    
    float F_tot = force[X]*force[X] + force[Y]+force[Y];
    float force_mult = 2*dt/(mass*2);
    if(F_tot > force_mult){
        coord_scale((force_mult)/F_tot, force);
    }

    
    
    //force points to center:
    attraction_force[X] = Sim.getArenaWidth()/2;
    attraction_force[Y] = Sim.getArenaHeight()/2;
    coord_diff(w_0, attraction_force, attraction_force);
    coord_scale(-5*sqrtf(coord_mult(attraction_force, attraction_force)), attraction_force);
    
    coord_add(attraction_force, force, force);
    
    
    // Move - bounce if hits wall

    /* Marble model - go straight until collision with wall */
     //externalize if you want to mess with update rates
     coord tmp;
     tmp[X] = vel[X] + 0.5*dt*force[X]/mass;
     tmp[Y] = vel[Y] + 0.5*dt*force[Y]/mass; // dx = v + 1/2a dt
     vel[X] += dt*force[X]/mass;
     vel[Y] += dt*force[Y]/mass;
    float velocity = sqrt(vel[X]*vel[X] + vel[Y]+vel[Y]);
    if(velocity > 4){
        coord_scale(4/velocity, vel);
    }
    
     if(!tmp[X] && !tmp[Y]) return; //don't bother computing if vel == {0,0}
     coord_scale(dt, tmp); //tmp is now the change in position.
     if((w_0[X] + tmp[X]) > arena_x){ //coord system has (0,0) as lower left corner
         vel[X] = -vel[X]; //bounce on X
         w_0[X] += 2*(arena_x-w_0[X]) - tmp[X];
     } else if((w_0[X] + tmp[X]) < 0){
         vel[X] = -vel[X];
         w_0[X] -= 2*w_0[X]+tmp[X];
     } else {
         w_0[X] += tmp[X];
     }
     
     if((w_0[Y] + tmp[Y]) > arena_y){ //coord system has (0,0) as lower left corner
         vel[Y] = -vel[Y]; //bounce on Y
         w_0[Y] += 2*(arena_x-w_0[Y]) - tmp[Y];
     } else if((w_0[Y] + tmp[Y]) < 0){
         vel[Y] = -vel[Y];
         w_0[Y] -= 2*w_0[Y]+tmp[Y];
     } else {
         w_0[Y] += tmp[Y];
     }
}

void Node::set_velocity(float vx, float vy)
{
	vel[X] = vx;
	vel[Y] = vy;
}