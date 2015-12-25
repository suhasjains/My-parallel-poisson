#ifndef _PARALLEL_H
#define _PARALLEL_H
#include "mpi.h"

namespace poisson {

//Parallelization variables
int proc_rank;                                                      //rank of the current process
int P_grid_rank;                                                    //rank of the current proces in the virtual grid
int P_grid_coord[2];                                                //coordinates of the process in the virtual grid
int P_grid_top, P_grid_bottom, P_grid_left, P_grid_right;           //ranks of the neighbouring processes
int n_Procs;                                                        //total number of processes
int P_grid[2];                                                      //virtual grid dimensions
int offset[2];                                                      //offset for cell numbering for subdomains
MPI_Datatype exchange_buffer_type[2];                               //MPI_datatype for exchange of buffer cell data
MPI_Comm grid_comm;                                                 //grid COMMUNICATOR

MPI_Status status;


//Creates a virtual cartesian topology
void setup_proc_grid();

//Packs data for communication across processes
void setup_MPI_datatypes(int, int, int);

//exchanges data between processes
void exchange_buffers(Field *, int, int);





}
#endif
