#ifndef _PARALLEL_POISSON_H
#define _PARALLEL_POISSON_H
#include <stdio.h>
#include <cmath>
//#include "parallel.h"
#include<string.h>      //for memcpy
#include<string>        //for strings
//#include "output.cpp"
#include "mpi.h"


namespace poisson {

//Directions
#define EAST i+1
#define WEST i-1
#define NORTH i+Nx
#define SOUTH i-Nx
#define P i



//Boundary condition type
typedef enum{
        NONE,
        DIRICHLET,
        BUFFER,
} BC_type;

//Direction
typedef enum{
        X_DIR,
        Y_DIR,
} Direction;


//Field class
class Field {

        public:
        //Members
                int Nx,Ny;              //size
                int N;                  //size
                double *val;            //values
                BC_type *bc;            //BC type
                double bc_val[5];       //BC value

        //Member functions
        void set_field_value(double value) {

                for(int i=0;i<N;i++){
                        this->val[i] = value;
                }
        }


        //Constructors
        //allocates memory to the field variables equal to the number of cells in the domain
        Field( int N_x, int N_y ) : Nx(N_x), Ny(N_y) {

                N = N_x*N_y;
                val = new double [N];
                bc = new BC_type [N];
        }

        Field() {

                Nx = 0;
                Ny = 0;
                N = Nx*Ny;
                val = new double [N];
                bc = new BC_type [N];
        }


        //Copy constructor
        Field(const Field &obj) {


                Nx = obj.Nx;
                Ny = obj.Ny;
                N = obj.N;
                val = new double [N];
                memcpy(val,obj.val,sizeof(double)*N);
                bc = new BC_type [N];
                memcpy(bc,obj.bc,sizeof(BC_type)*N);
        }

	//Destructor
        ~Field() {

                delete val;
                delete  bc;
        }

        private:
        protected:

};


//Patch
typedef enum{
        INSIDE,
        XMIN,
        XMAX,
        YMIN,
        YMAX
} PATCH_type;
extern PATCH_type patch;


//Constant class
class Constant{

        public:
                double h;
                double f;

        private:
        protected:

};

//Domain class
class Domain{

        public:
        //Members
                Field *u;
                Constant *constant;

        //Constructors
        //adds the reference to the field and constant to the domain
        Domain( Field* field, Constant* con) : u(field), constant(con) {

        }

        private:
        protected:

};

//Parallelization variables
extern int proc_rank;                                                      //rank of the current process
extern int P_grid_rank;                                                    //rank of the current proces in the virtual grid
extern int P_grid_top, P_grid_bottom, P_grid_left, P_grid_right;           //ranks of the neighbouring processes
extern int n_Procs;                                                        //total number of processes
extern int P_grid[2];                                                      //virtual grid dimensions
extern MPI_Datatype exchange_buffer_type[2];                               //MPI_datatype for exchange of buffer cell data
extern MPI_Comm grid_comm;                                                 //grid COMMUNICATOR
extern int offset[2];                                                      //offset for cell numbering for subdomains
extern int P_grid_coord[2];                                                //coordinates of the process in the virtual grid
extern MPI_Status status;

//sets flags for ghost and buffer cells
void set_ghost_buffer_flag(Domain);

//Jacobi solver
void jacobi(Field, int, int, Constant);

//Gauss seidel solver 
void gauss_seidel(Field, int, int, Constant);

//Steepest descent solver
void solve_SD(Field*, Constant); 

//Conjugate gradient solver
void solve_CG(Field*, Constant); 

//AX maxtrix multiplication
void compute_AX(Field*, double*, Constant);

//Sets boundary condition values to boundary cells
void set_bc(Field);

//Creates a virtual cartesian topology
void setup_proc_grid();

//Packs data for communication across processes
void setup_MPI_datatypes(int, int, int);

//exchanges data between processes
void exchange_buffers(Field *, int, int);

void write_output(Domain, int);

}
#endif
