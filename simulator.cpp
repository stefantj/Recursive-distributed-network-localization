/*
 *  simulator.cpp
 *
 *  Created by Stefan Jorgensen for Stanford EE373A project (Winter 2016).
 *  based on earlier work for EE359.
 *
 *  Code released under the MIT License (MIT)
 *
 */

#include "simulator.h"
#include "Node.h"



using namespace std;

//Node list:
extern Node* NODES[NETWORK_SIZE];

//Neighborhood container
//first index: k,
//second index:
// first element: |N_k|
// next |N_k|-1 elements: indices of neighbors
extern int NEIGHBORS[NETWORK_SIZE][NETWORK_SIZE+1];


Simulator::Simulator()
{
    // Initialize random number generator
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
//    gsl_rng_set(RNG, static_cast<unsigned int>(time(0)));
    gsl_rng_set(RNG, static_cast<unsigned int>(0));
    
	//Load Default simulator values:
	for(int j = 0; j < ITER_MAX; j++)
		NET_MSE[j] = 0;
	ARENA_WIDTH = ARENA_SIZE_DEFAULT;
	ARENA_HEIGHT = ARENA_WIDTH;
	M_RADIUS = 30;
	N_RADIUS = 100;
	STDEVS[SDEV_R] = 0.5;
	STDEVS[SDEV_U] = .025;
	STDEVS[SDEV_G] = 3; 
	FAIL_MODE = BERNOULLI;
	LOOPS = 200;
	BETA = 0.3;
	NUM_TIMESLOTS = 10;
	SINR_THRESH = 1000; //mW

	Nbar_cnt = 0;
	Nbar_tot = 0;	
	Mbar = 0;
	
	char fnames[NUM_FILES][12] = {"000_w_0.txt","000_w_k.txt","000_mse.txt","000imse.txt","000_par.txt"};
	memcpy(filenames, fnames, 12*NUM_FILES);

    //Reset files:
    firstWrite = true;
    
	//Initialize simulator:
//	srand(static_cast<unsigned int>(time(0)));
    srand(static_cast<unsigned int>(0));

	
	printf("Initializing simulation...\n");
	
	printf("Creating %d nodes in a %d by %d m arena...\n", NETWORK_SIZE, ARENA_WIDTH, ARENA_HEIGHT);
	
	for( int i = 0; i < NETWORK_SIZE; i++)
		NODES[i] = new Node(i);
	printf("Creating neighborhoods...\n");
	update_neighbors();
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		printf("[%d] |%d|",j, NEIGHBORS[j][0]);
		for(int jj = 1; jj <= NEIGHBORS[j][0]; jj++)
		{
			printf(", %d", NEIGHBORS[j][jj]);
		}
		printf("\n");
	}
}

float Simulator::fixed_pos_trial(char* prefix, float param, float* stdevs)
{
	bool save_files = true;
	if(prefix == NULL)
	{	
		FAIL_MODE = BERNOULLI;
		save_files = false;
	} else if(prefix[0] == 'B')
	{
		//Load simulation parameters:
		FAIL_MODE = BERNOULLI;
		BETA = param;
	} else if(prefix[0] == 'S')
	{
		FAIL_MODE = SINR_FAIL;
		SINR_THRESH = param;	
	} else{
		printf("Unknown prefix %c\n", prefix[0]); 
		return -1;
	}

	STDEVS[SDEV_R] = stdevs[SDEV_R];
	STDEVS[SDEV_U] = stdevs[SDEV_U];
	STDEVS[SDEV_G] = stdevs[SDEV_G];
	//Save node positions (doesn't change)
	if(save_files)
	{
		save_scene(prefix, W_0_FILE);
		save_scene(prefix, PARAM_FILE);
	}
	for(int j = 0; j < LOOPS; j++)
	{
		printf("%d\n",j);
		newPositions();
		for(int jj = 0; jj < ILOOPS; jj++)
			run_sim(0, prefix);
	}
//	FILE* averages = fopen("averages.txt", "a");
//	fprintf(averages, "%f, %f, %f, %f\n", param, Nbar_tot/(float)Nbar_cnt, Mbar, NET_MSE[ITER_MAX-1]);
//	fclose(averages);
	
	printf("%f, %f, %f, %f\n", param, Nbar_tot/(float)Nbar_cnt, Mbar, NET_MSE[ITER_MAX-1]);
	//save_scene(prefix, MSE_FILE);
	//average past 10 samples:
	float nmse = 0;
	for(int i = ITER_MAX-10; i < ITER_MAX; nmse+=NET_MSE[i]/10, i++);
	return NET_MSE[ITER_MAX-1];
	reset_nodes(); //may want to move this outside
	reset_mse();
}

void Simulator::tracking_trial(char* prefix, float param, float* stdevs)
{
	FAIL_MODE = BERNOULLI;
	BETA = param; //units mW
	STDEVS[SDEV_R] = stdevs[SDEV_R];
	STDEVS[SDEV_U] = stdevs[SDEV_U];
	STDEVS[SDEV_G] = stdevs[SDEV_G];
	
	run_sim(1,prefix);
	reset_nodes();
	reset_mse();
}

//Data getter/setters:
void Simulator::setArenaWidth(int width)
{
	ARENA_WIDTH = width;
}
int Simulator::getArenaWidth()
{
	return ARENA_WIDTH;
}
void Simulator::setArenaHeight(int height)
{
	ARENA_HEIGHT = height;
}
int Simulator::getArenaHeight()
{
	return ARENA_HEIGHT;
}
int Simulator::getNradius()
{
	return N_RADIUS;
}
void Simulator::setMradius(int mrad)
{
	M_RADIUS = mrad;
}
int Simulator::getMradius()
{
	return M_RADIUS;
}
float Simulator::getMbar()
{
	return Mbar;
}
//int Simulator::testNoise(float var)


void Simulator::setLoops(float loops, float iloops)
{
	LOOPS = loops;
	ILOOPS = iloops;
}

//returns measurements wrt node k
// u = nj-nk
bool Simulator::get_meas(int j, int k, float* r, coord u)
{
	if(j != k)
	{
		float r_j;
		coord u_j;
		if((r_j = dist_nodes(j, k)) > M_RADIUS) return false;
		//within measuring distance
		
		//Pass noisy distance
		r[0] = r_j + gsl_ran_gaussian(RNG, STDEVS[SDEV_R]); //add measurement noise

		//calculate direction vector
		u_j[X] = (NODES[j]->get_w_0(X) - NODES[k]->get_w_0(X));
		u_j[Y] = (NODES[j]->get_w_0(Y) - NODES[k]->get_w_0(Y));
		float n_theta = gsl_ran_gaussian(RNG, STDEVS[SDEV_U]);
		coord_rotate(n_theta, u_j);
		
		r_j = sqrt(coord_mult(u_j, u_j)); //holds magnitude sq of u_j
		if(r_j)
		{
			u[X] = u_j[X]/r_j;
			u[Y] = u_j[Y]/r_j;
		} else {printf("Error: u vec = 0\n"); return false;}
		return true;
	} else {
		//pass noisy version of self:
		u[X] = NODES[j]->get_w_0(X) + gsl_ran_gaussian(RNG, STDEVS[SDEV_G]);
		u[Y] = NODES[j]->get_w_0(Y) + gsl_ran_gaussian(RNG, STDEVS[SDEV_G]);
		return true;
	}
	
}


//mixes up positions
void Simulator::newPositions()
{
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		NODES[j]->randomize_position();
	}	
}

void Simulator::reset_mse()
{
	Nbar_cnt = 0;
	Nbar_tot = 0;
	Mbar = 0;

	for(int j = 0; j < ITER_MAX; j++)
		NET_MSE[j] = 0;
}

//*** Private functions ***//

void Simulator::reset_nodes()
{
	
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		NODES[j]->initialize();
	}
}	

void Simulator::run_sim(int simtype, const char* prefix)
{
	int M2bar = 0;
	if(simtype == 0)
	{
		for(int I = 0; I < ITER_MAX; I++)
		{
			for(int j = 0; j < NETWORK_SIZE; j++)
				M2bar += NODES[j]->measure();
			update_neighbors();
			for(int j = 0; j < NETWORK_SIZE; j++)
				Mbar += NODES[j]->calc_psi();
			update_neighbors();
			for(int j = 0; j < NETWORK_SIZE; j++)
				NODES[j]->calc_w();
			NET_MSE[I] += get_NETMSE();
		}
		Mbar /= (NETWORK_SIZE*ITER_MAX);
	} else if (simtype == 1) //mobile simulation for demonstration/tracking performance.
	{
		for(int I = 0; I < ITER_MAX; I++) 
		{
			//extra update_position step
			for(int j = 0; j < NETWORK_SIZE; j++)
				NODES[j]->update_position(ARENA_WIDTH, ARENA_HEIGHT);
			for(int j = 0; j < NETWORK_SIZE; j++)
				M2bar += NODES[j]->measure()/(float)NETWORK_SIZE;
			update_neighbors();
			for(int j = 0; j < NETWORK_SIZE; j++)
				Mbar += NODES[j]->calc_psi();
			update_neighbors();
			for(int j = 0; j < NETWORK_SIZE; j++)
				NODES[j]->calc_w();
            
			//Typically this is run with 1 realization, so OK to save here.
			//if too slow, save to array and then print array to file outside loop
			save_scene(prefix, W_0_FILE); //Save true locations
			save_scene(prefix, W_K_FILE); //Save node 0's estimates
            save_julia_data(I == (ITER_MAX - 1));
			printf("%f, (%f,%f)\n",get_NETMSE()*LOOPS*ILOOPS, NODES[0]->get_w_0(X), NODES[0]->get_w_0(Y));
			Mbar= 0;
			NET_MSE[I] += get_NETMSE();
            if(I != 0) firstWrite = false;
		}
		if(prefix != NULL)
			save_scene(prefix, MSE_FILE); //Save MSE data
	}	
}


void Simulator::update_neighbors()
{
	//initialize cardinality:
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		NEIGHBORS[j][1] = j;
		NEIGHBORS[j][0] = 1;
	}
	//always in its own neighborhood
	
	//compute connections
	
	switch (FAIL_MODE) {
		case BERNOULLI:
			for(int j = 0; j < NETWORK_SIZE; j++)
			{
				for(int l = 0; l < NETWORK_SIZE; l++)
				{
					if(j!= l && dist_nodes(j, l) < N_RADIUS)
					{	
						int b = rand()%10000;
						if(b > BETA*10000)
						{
							NEIGHBORS[j][0]++;
							NEIGHBORS[j][NEIGHBORS[j][0]] = l;
						}
					}
				}			
			}
			break;
		case DEPENDENT_BERNOULLI:
			for(int j = 0; j < NETWORK_SIZE; j++)
			{
				//compute number of nearest neighbors:
				int N = 0;
				for(int l = 0; l < NETWORK_SIZE; l++)
					if(dist_nodes(j, l)<N_RADIUS)
						N++;
				//compute failed connections:
				for(int l = 0; l < NETWORK_SIZE; l++)
				{
					if(j!= l && dist_nodes(j, l) < N_RADIUS)
					{	
						int b = rand()%10000;
						if(b > N*BETA*10000)
						{
							NEIGHBORS[j][0]++;
							NEIGHBORS[j][NEIGHBORS[j][0]] = l;
						}
					}
				}			
			}
			break;
			
		case SINR_FAIL:
			//Assign timeslot to each
			int timeslots[NETWORK_SIZE];
			for(int j = 0; j< NETWORK_SIZE; j++)
				timeslots[j] = rand()%NUM_TIMESLOTS; 
			//compute noise and interference:
			for(int t = 0; t < NUM_TIMESLOTS; t++)
			{	
				for(int j = 0; j < NETWORK_SIZE; j++)
				{
					int neighbs[NETWORK_SIZE] = {0};
					float neighbspwr[NETWORK_SIZE] = {0};
					float noisepwr = 0; 
					//can't listen when talking
					if(timeslots[j] != t)
					{
						for(int l = 0; l < NETWORK_SIZE; l++)
						{
							if(l!=j && timeslots[l] == t)
							{
								//note who is here
								neighbs[0]++;
								neighbs[neighbs[0]] = l;
								float tmp = dist_nodes(j, l);
								//super duper basic 1/r^2 model
								tmp = NODES[l]->get_power()/(tmp*tmp);
								neighbspwr[l] = tmp;
								noisepwr += tmp;
							}
						}
						for(int l = 1; l <= neighbs[0]; l++)
						{
							int k = neighbs[l];
							float tmp = neighbspwr[k]/noisepwr;
							if(tmp > SINR_THRESH)
							{
								NEIGHBORS[j][0]++;
								NEIGHBORS[j][NEIGHBORS[j][0]] = k;
							}
						}
					}
				}
			}
			
			break;
			
		default:
			break;
	}
	
	//statistics for analysis:
	float ntot = 0;
	for(int j = 0; j < NETWORK_SIZE; j++)
		ntot += NEIGHBORS[j][0]; //contains the number of neighbors
	ntot/=NETWORK_SIZE; //average size of the network at time time instant
	Nbar_tot += ntot; //sum
	Nbar_cnt++; //number	
}

float Simulator::dist_nodes(int j, int l)
{
	coord n1, n2;
	n1[X] = NODES[j]->get_w_0(X);
	n1[Y] = NODES[j]->get_w_0(Y);
	n2[X] = NODES[l]->get_w_0(X);
	n2[Y] = NODES[l]->get_w_0(Y);
	coord_diff(n1, n2, n1);
	return (float)sqrt(coord_mult(n1, n1));
}

float Simulator::get_NETMSE()
{
	float mse = 0;
	float x1, x2, y1, y2;
	float xerr, yerr;
	
	for(int j = 0; j < NETWORK_SIZE; j++)
	{
		x1 = NODES[j]->get_w_jk(j, X);
		x2 = NODES[j]->get_w_0(X);
		xerr = x1-x2;
		y1 = NODES[j]->get_w_jk(j, Y);
		y2 = NODES[j]->get_w_0(Y);
		yerr = y1-y2;
		mse += (xerr*xerr + yerr*yerr)/NETWORK_SIZE;
	}
	return mse/(LOOPS*ILOOPS);
}

//since there is so much data, save this in binary format
void Simulator::save_scene(const char* prefix, int file_index)
{
	char* tmp;
	float mse;
    char* writeType = "a";
    if(firstWrite)
      writeType = "w";
    
	switch (file_index) {
		case W_0_FILE:	// Save true positions:
			tmp = filenames[W_0_FILE];
			tmp[0] += prefix[0] - '0';
			tmp[1] += prefix[1] - '0';
			tmp[2] += prefix[2] - '0';
			
			files[W_0_FILE] = fopen(tmp, writeType);
			
			for( int i = 0; i < NETWORK_SIZE; i++)
			{
				coord pos; 
				pos[X] = NODES[i]->get_w_0(X);
				pos[Y] = NODES[i]->get_w_0(Y);
//				fwrite(pos,1, sizeof(pos), files[W_0_FILE]);
                fwrite(pos, sizeof(pos), 1, files[W_0_FILE]);
			}
			fclose(files[W_0_FILE]);
			tmp[0] -= prefix[0] - '0';
			tmp[1] -= prefix[1] - '0';
			tmp[2] -= prefix[2] - '0';
			break;
		case W_K_FILE:	//Save positions seen by node 0:
			tmp = filenames[W_K_FILE];
			tmp[0] += prefix[0] - '0';
			tmp[1] += prefix[1] - '0';
			tmp[2] += prefix[2] - '0';
			
			files[W_K_FILE] = fopen(tmp, writeType);
			for( int i = 0; i < NETWORK_SIZE; i++)
			{
				coord pos; 
				pos[X] = NODES[0]->get_w_jk(i,X); 
				pos[Y] = NODES[0]->get_w_jk(i,Y); 
				fwrite(pos, sizeof(float), 2, files[W_K_FILE]);
			}
			fclose(files[W_K_FILE]);
			tmp[0] -= prefix[0] - '0';
			tmp[1] -= prefix[1] - '0';
			tmp[2] -= prefix[2] - '0';			
			break;
		case MSE_FILE: //Record averaged Net MSE
			tmp = filenames[MSE_FILE];
			tmp[0] += prefix[0] - '0';
			tmp[1] += prefix[1] - '0';
			tmp[2] += prefix[2] - '0';
			
			files[MSE_FILE] = fopen(tmp, writeType);
			fwrite(NET_MSE, sizeof(float),ITER_MAX,files[MSE_FILE]);
			fclose(files[MSE_FILE]);
			tmp[0] -= prefix[0] - '0';
			tmp[1] -= prefix[1] - '0';
			tmp[2] -= prefix[2] - '0';			
			break;
						
		case INST_MSE: //record instantaneous Net MSE
			tmp = filenames[INST_MSE];
			tmp[0] += prefix[0] - '0';
			tmp[1] += prefix[1] - '0';
			tmp[2] += prefix[2] - '0';
			
			files[INST_MSE] = fopen(tmp, writeType);
			mse = get_NETMSE();
			fwrite(&mse, sizeof(float),1,files[INST_MSE]);
			fclose(files[INST_MSE]);
			tmp[0] -= prefix[0] - '0';
			tmp[1] -= prefix[1] - '0';
			tmp[2] -= prefix[2] - '0';			
			break;

		case PARAM_FILE://record parameter file (human readable??)
			printf("Saving parameter file\n");
			tmp = filenames[PARAM_FILE];
			tmp[0] += prefix[0] - '0';
			tmp[1] += prefix[1] - '0';
			tmp[2] += prefix[2] - '0';
			files[PARAM_FILE] = fopen(tmp, writeType);
			//prefix describes fail mode, 
			//Parameter data: beta/thresh, arena size, network size, M rad, N rad, sdevs, itermax
			float params[9]; 
			if(prefix[0] == 'B') //bernoulli
			{params[0] = BETA; }
			else if(prefix[0] == 'S') //SINR
			{params[0] = SINR_THRESH;}
			else{printf("Error saving parameters: tmp = %c\n", tmp[0]);}
			params[1] = (float)ARENA_WIDTH;
			params[2] = (float)NETWORK_SIZE;
			params[3] = (float)M_RADIUS;
			params[4] = (float)N_RADIUS;
			params[5] = STDEVS[SDEV_R];
			params[6] = STDEVS[SDEV_U];
			params[7] = STDEVS[SDEV_G];
			params[8] = (float)ITER_MAX;
			fwrite(params, sizeof(float), 9, files[PARAM_FILE]);
			fclose(files[PARAM_FILE]);
			printf("%f %f %f %f %f %f %f %f %f\n", params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8]);
			tmp[0] -= prefix[0] - '0';
			tmp[1] -= prefix[1] - '0';
			tmp[2] -= prefix[2] - '0';
			break;
		default:
			printf("Filename not found");
			break;
	}
}

// Function to save data in an ascii format for julia to parse
void Simulator::save_julia_data(bool close){
    char* writeType = "a";
    if(firstWrite)
    writeType = "w";
    FILE* msefile_handle = fopen("mse.jl", writeType);
    
    FILE* xfile_handle = fopen("data_x.jl",writeType);
    FILE* yfile_handle = fopen("data_y.jl",writeType);
    if(firstWrite){
        fprintf(xfile_handle, "x_pos = [");
        fprintf(yfile_handle, "y_pos = [");
        fprintf(msefile_handle, "mse = [");
    }
 
    fprintf(msefile_handle, "%f ", get_NETMSE()*LOOPS*ILOOPS);
    for( int i = 0; i < NETWORK_SIZE-1; i++){
        fprintf(xfile_handle, "%f ", NODES[i]->get_w_0(X));
        fprintf(yfile_handle, "%f ", NODES[i]->get_w_0(Y));
    }

    // print last line
    if(close){
        fprintf(msefile_handle, "]");
        fprintf(xfile_handle, "%f]", NODES[NETWORK_SIZE-1]->get_w_0(X));
        fprintf(yfile_handle, "%f]", NODES[NETWORK_SIZE-1]->get_w_0(Y));
    } else{
        fprintf(msefile_handle, ",");
        fprintf(xfile_handle, "%f;\n", NODES[NETWORK_SIZE-1]->get_w_0(X));
        fprintf(yfile_handle, "%f;\n", NODES[NETWORK_SIZE-1]->get_w_0(Y));
    }
    fclose(msefile_handle);
    fclose(xfile_handle);
    fclose(yfile_handle);
    
}