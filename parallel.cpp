//#include "parallel_poisson.h"
#include "parallel.h"
namespace poisson {

//Creating a virtual cartesian topology
void setup_proc_grid() {

        int wrap_around[2];
        int reorder;

        //retrieve the number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &n_Procs);

        //number of processes per row and column
        P_grid[X_DIR] = 1;
        P_grid[Y_DIR] = 1;

        if(P_grid[X_DIR] * P_grid[Y_DIR] != n_Procs)    printf("Error: Mismatch of number of processes and process grid");

        //creating a virtual process grid topology
        wrap_around[X_DIR] = 0;
        wrap_around[Y_DIR] = 0;         //dont connect first and last processes
        reorder = 1;                    //reorder process ranks

        //creates a new communicator, grid_comm
        MPI_Cart_create(MPI_COMM_WORLD, 2, P_grid, wrap_around, reorder, &grid_comm);

        //retrieve new rank and cartesian coordinates of this process
        MPI_Comm_rank(grid_comm, &P_grid_rank);
        MPI_Cart_coords(grid_comm, P_grid_rank, 2, P_grid_coord);

        //calculate ranks of neighboring processes in the grid
        MPI_Cart_shift(grid_comm, 1, 1, &P_grid_left, &P_grid_right);
        MPI_Cart_shift(grid_comm, 0, 1, &P_grid_bottom, &P_grid_top);

}


//Packing data for communication across processes
void setup_MPI_datatypes(int N_int_x, int N_int_y, int N_x){


        //Datatype for horizontal data exchange
        MPI_Type_vector(N_int_y, 1, N_x, MPI_DOUBLE, &exchange_buffer_type[X_DIR]);
        MPI_Type_commit(&exchange_buffer_type[X_DIR]);

        //Datatype for vertical data exchange
        MPI_Type_vector(N_int_x, 1, 1, MPI_DOUBLE, &exchange_buffer_type[Y_DIR]);
        MPI_Type_commit(&exchange_buffer_type[Y_DIR]);
}


//communication between the processes
void exchange_buffers(Field *phi, int N_local_x, int N_local_y){


        //All traffic in direction "top"
         MPI_Sendrecv(&(phi->val[(N_local_y-2)*(N_local_x)+1]), 1, exchange_buffer_type[Y_DIR], P_grid_top, 0, &(phi->val[1]), 1, exchange_buffer_type[Y_DIR], P_grid_bottom, 0, grid_comm, &status);

        //All traffic in direction "bottom"
        MPI_Sendrecv(&(phi->val[N_local_x+1]), 1, exchange_buffer_type[Y_DIR], P_grid_bottom, 0, &(phi->val[(N_local_y-1)*(N_local_x)+1]), 1, exchange_buffer_type[Y_DIR], P_grid_top, 0, grid_comm, &status);

        //All traffic in direction "left"
        MPI_Sendrecv(&(phi->val[N_local_x+1]), 1, exchange_buffer_type[X_DIR], P_grid_left, 0, &(phi->val[(2*N_local_x)-1]), 1, exchange_buffer_type[X_DIR], P_grid_right, 0, grid_comm, &status);

        //All traffic in direction "right"
        MPI_Sendrecv(&(phi->val[2*(N_local_x-1)]), 1, exchange_buffer_type[X_DIR], P_grid_right, 0, &(phi->val[N_local_x]), 1, exchange_buffer_type[X_DIR], P_grid_left, 0, grid_comm, &status);

}



}
