/*===================================================/
/  Comms routine for PX425 Assignment 5. Some basic  /
/  routines/ideas lifted from Assignment 4 code.     /
/  All interactions with MPI libraries are in here.  /
/                                                    /
/  Original code by C. Woodgate                      /
/  (based on code by N. Hine, D. Quigley)            /
/===================================================*/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "trajectories.h"
#include "comms.h"

int p;               /* Number of processor       */
int my_rank;         /* Rank of current processor */

MPI_Comm comm = MPI_COMM_WORLD;  /* Cartesian communicator    */

double t1,t2;		/* Time of initialisation and shutdown */

/*============/
/  FUNCTIONS  /
/============*/

void comms_wait()
{
	/*==================================================================/
	/ Function to allow me to call MPI_Barrier to sync threads.         /
	/==================================================================*/
	
	/* Call MPI_Barrier */
	MPI_Barrier(comm);
	return;
}

void comms_initialise(int argc, char **argv) 
{
	/*==================================================================/
	/ Function to initialise MPI, get the communicator size p and the   /
	/ rank my_rank of the current process within that communicator.     /
	/==================================================================*/

	// Initialise MPI 
	MPI_Init(&argc, &argv);

	// Get the rank of the current process
	MPI_Comm_rank(comm, &my_rank);

	// Get the size of the current communicator
	MPI_Comm_size (comm, &p);

	// Start the timer
	t1 = MPI_Wtime();

	return;
}

void comms_phase_space_local(int my_rank, int nFragments, int n_local)
{
	/*==================================================================/
	/ Function to set up the phase space for the system on each thread  /
	/ in 'local_' arrays from rank 0, where we called explodeTheMoon(). /
	/==================================================================*/
	
	//Number of elements to be sent when vectors
	int n_send = 3*n_local;
	
	// Broadcast positions array to all ranks
	MPI_Bcast(pos, nFragments*3, MPI_DOUBLE, 0, comm);

	//Boradcast mass information to all ranks
	MPI_Bcast(mass, nFragments, MPI_DOUBLE, 0, comm);
	
	//Give each rank information about whether fragment is a dummy or not
	MPI_Bcast(calc_active, nFragments, MPI_INT, 0, comm);
	
	// Send each rank the section of the velocity array on which they will work
	MPI_Scatter(global_vel, n_send, MPI_DOUBLE, local_vel, n_send, MPI_DOUBLE, 0, comm);

	// Send each rank their section of the acceleration array.
	// (Almost redundant but as we only call it once doesn't cost too much)
	MPI_Scatter(global_acc, n_send, MPI_DOUBLE, local_acc, n_send, MPI_DOUBLE, 0, comm);

	// Send each rank their relevant radii
	MPI_Scatter(global_radius, n_local, MPI_DOUBLE, local_radius, n_local, MPI_DOUBLE, 0, comm);

	// Send each rank the 'active' array
	MPI_Scatter(global_active, n_local, MPI_INT, local_active, n_local, MPI_INT, 0, comm);
	
	return;
}

void comms_positions_gather(int my_rank, int n_alloc, int n_local)
{
	/*==================================================================/
	/ Function to gather the positions array onto each thread after     /
	/ each timestep using MPI_Allgather.                                /
	/==================================================================*/

	// Number of elements to be sent
	int n_send = 3*n_local;

	// Wait until all threads have finished their calculations
	MPI_Barrier(comm);

	// Gather positions array on all threads
	MPI_Allgather(local_pos, n_send, MPI_DOUBLE, pos, n_send, MPI_DOUBLE, comm); 

	return;
}

void comms_impacts(int n_local)
{
	/*==================================================================/
	/ Function to communicate information needed about impacts on       /
	/ each thread.                                                      /
	/==================================================================*/

	// Work out how many fragments are active in each section	
	MPI_Reduce(&local_nActive, &nActive, 1, MPI_INT, MPI_SUM, 0, comm);

	// Work out which fragment is closest overall
	MPI_Reduce(&local_minRadsq, &minRadsq, 1, MPI_DOUBLE, MPI_MIN,0, comm);

	// Gather the local_active arrays onto thread zero
	MPI_Gather(local_active, n_local, MPI_INT, global_active, n_local, MPI_INT, 0, comm);
	return;
}

void comms_finalise() 
{
	/*===========================================================/
	/ Function to finalise MPI functionality and exit cleanly    /
	/===========================================================*/
  
	// Measure the time t2 using MPI_Wtime()
	t2 = MPI_Wtime();

	if (my_rank==0) 
	{
    	printf("Total time elapsed since MPI initialised :  %12.6f s\n",t2-t1);
	}

	// Shutdown MPI
	MPI_Finalize();

	return;
}

