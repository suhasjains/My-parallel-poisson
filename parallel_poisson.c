#include <stdio.h>
#include <cmath>
#include "mpi.h"


/********************************************************
Solving a poisson equation in 2D serially
Finite volume discretized equation used is UE + UW + UN + US - 4*UP = f*h^2

 h is the spacing. 
 f is the forcing function.
 Nx and Ny are the number of CVs in x and y direction respectively_
***********************************************************/

#define EAST i+1
#define WEST i-1
#define NORTH i+Nx
#define SOUTH i-Nx
#define P i



int proc_rank; 							//rank of the current process
int P_grid_rank;							//rank of the current proces in the virtual grid
int P_grid_coord[2]; 						//coordinates of the process in the virtual grid
int P_grid_top, P_grid_bottom, P_grid_left, P_grid_right;	//ranks of the neighbouring processes
int n_Procs;							//total number of processes
int P_grid[2];							//virtual grid dimensions
int offset[2];							//offset for cell numbering for subdomains
MPI_Comm grid_comm;						//grid COMMUNICATOR

MPI_Status status;



//defining enums and structures
typedef enum{
	NONE,
	DIRICHLET,
	BUFFER,
} BC_type;


typedef enum{
	X_DIR,
	Y_DIR,
} Direction;


typedef struct _Field Field;
struct _Field{
	int Nx,Ny;
	int N;
	double *val;
	BC_type *bc;
	double bc_val[5];
};	


typedef enum{
	INSIDE,
	XMIN,
	XMAX,
	YMIN,
	YMAX
} PATCH_type;
PATCH_type patch;


typedef struct _Domain Domain;
struct _Domain{
	Field *u;
};

typedef struct _Constant Constant;
struct _Constant{
	double h;
	double f;
};




//declaring functions
//void set_ghosts(Domain);
//void set_bc(Field);
static Field *allocate_field(int,int);
//void jacobi1(Field , int, int, Constant);

void Setup_proc_grid(){

	int wrap_around[2];
	int reorder;
	
	//retrieve the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &n_Procs);

	//number of processes per row and column
	P_grid[X_DIR] = 2;
	P_grid[Y_DIR] = 2;

	if(P_grid[X_DIR] * P_grid[Y_DIR] != n_Procs) 	printf("Error: Mismatch of number of processes and process grid");

	//creating a virtual process grid topology
	wrap_around[X_DIR] = 0;
	wrap_around[Y_DIR] = 0;		//dont connect first and last processes
	reorder = 1;			//reorder process ranks
	
	//creates a new communicator, grid_comm
	MPI_Cart_create(MPI_COMM_WORLD, 2, P_grid, wrap_around, reorder, &grid_comm);

	//retrieve new rank and cartesian coordinates of this process 
	MPI_Comm_rank(grid_comm, &P_grid_rank); 
	MPI_Cart_coords(grid_comm, P_grid_rank, 2, P_grid_coord);

//	printf("I am %i and my location is %i, %i\n", P_grid_rank, P_grid_coord[0], P_grid_coord[1]);
   
	//calculate ranks of neighboring processes in the grid
	MPI_Cart_shift(grid_comm, 1, 1, &P_grid_left, &P_grid_right);
	MPI_Cart_shift(grid_comm, 0, 1, &P_grid_bottom, &P_grid_top);		

//	if(P_grid_top == MPI_PROC_NULL)		printf("My top neighbor doesnt exist");
	
//	if(P_grid_bottom == MPI_PROC_NULL)	printf("My bottom neighbor doesnt exist"); 

//	if(P_grid_right == MPI_PROC_NULL)	printf("My left neighbor doesnt exist"); 

//	if(P_grid_left == MPI_PROC_NULL)	printf("My right neighbor doesnt exist"); 
	

//	printf("My rank is %d and my neighbors are top %d, bottom %d, right %d, left %d\n", P_grid_rank, P_grid_top, P_grid_bottom, P_grid_right, P_grid_left);

}

void set_ghost_buffer_flag(Domain domain){
	
	int i, l, m;

	Field *u = domain.u;

	int N_Cells = u->Nx * u->Ny;
	int N_Cells_x = u->Nx;
	int N_Cells_y = u->Ny;

	for(i=0;i<N_Cells;i++){
		l = i%N_Cells_x;	//gives the position of the CV along x
		m = (int) i/N_Cells_x;	//gives the position of the CV along y
		if((l==0&&P_grid_left==MPI_PROC_NULL) || (l==N_Cells_x-1&&P_grid_right==MPI_PROC_NULL) || (m==0&&P_grid_bottom==MPI_PROC_NULL) || (m==N_Cells_y-1&&P_grid_top==MPI_PROC_NULL)){
			u->bc[i] = DIRICHLET;
		}else if(l==0||l==N_Cells_x-1||m==0||m==N_Cells_y-1){
			u->bc[i] = BUFFER;
		}else{
			u->bc[i] = NONE;
		}
	}return;
}

void jacobi(Field *phi, int Nx, int Ny, Constant constant){

	double res, e;	
	int i, l, m;
	int N_Cells_x = Nx;
	int N_Cells_y = Ny;
	int N = Nx * Ny;
	int N_Cells = N;
	
	long int loop = 0;

	//Allocating memory to the members of the temporary Field temp (temp is not part of the domain)
	Field *temp = allocate_field( N_Cells_x, N_Cells_y );

	//making res 1 so that it enters into the loop	
	res = 1.0;
	
	//Starting the iteration loop
	//for(t=0;t<10000;t++){
	while(res > pow(10,-10)){	
	
	//making res 0 so that any error greater than 0 can be equated to this
	res = 0.0;

		//Making the temp field zero after every iteration
		for(i=0;i<N_Cells;i++){
			temp->val[i] = 0.0;
		}


		double u_E, u_W, u_N, u_S, u_P;
	
		for(i=0;i<N_Cells;i++){
			if(phi->bc[i] == NONE){ 
				l=i%N_Cells_x;
				m=(int) i/N_Cells_x;

				u_E = phi->val[EAST];
				u_W = phi->val[WEST];
				u_N = phi->val[NORTH];
				u_S = phi->val[SOUTH];
				u_P = phi->val[P];
			
				if(l==1&&P_grid_left==MPI_PROC_NULL)		temp->val[P] += 2.0*u_W - u_P;
				else 	 					temp->val[P] += u_W;
				if(l==N_Cells_x-2&&P_grid_right==MPI_PROC_NULL) temp->val[P] += 2.0*u_E - u_P;
				else 						temp->val[P] += u_E;
				if(m==1&&P_grid_bottom==MPI_PROC_NULL)		temp->val[P] += 2.0*u_S - u_P;
				else 						temp->val[P] += u_S;
				if(m==N_Cells_y-2&&P_grid_top==MPI_PROC_NULL)	temp->val[P] += 2.0*u_N - u_P;
				else 						temp->val[P] += u_N;
 
				temp->val[P] -= constant.f*pow(constant.h,2);
				temp->val[P] = temp->val[P]/4.0;
	
				e = temp->val[P] - phi->val[P];
				if(e > res)	res = e; 

			}
		}
	
		//Transferring values from temp to u
		for(i=0;i<N_Cells;i++){
			if(phi->bc[i] == NONE){
				phi->val[i] = temp->val[i];
			}
		}
		
		loop++;
	}

	printf("Maximum residual is %e and number of iterations are %ld and I am process %d \n", res, loop, proc_rank);	
	return;
}

void gauss_seidel(Field *phi, int Nx, int Ny, Constant constant){

	double res, e;	
	int i, l, m;
	long int loop = 0;
	int N_Cells_x = Nx;
	int N_Cells_y = Ny;
	int N = Nx * Ny;
	int N_Cells = N;

	double relax = 1.95;

	//Allocating memory to the members of the temporary Field temp (temp is not part of the domain)
	Field *temp = allocate_field( N_Cells_x, N_Cells_y );

	//making res 1 so that it enters into the loop	
	res = 1.0;
	
	//Starting the iteration loop
	//for(t=0;t<10000;t++){
	while(res > pow(10,-10)){	
	
	//making res 0 so that any error greated than 0 can be equated to this
	res = 0.0;

		//Making the temp field zero after every iteration
		for(i=0;i<N_Cells;i++){
			temp->val[i] = 0.0;
		}


		double u_E, u_W, u_N, u_S, u_P;
	
		for(i=0;i<N_Cells;i++){
			if(phi->bc[i] == NONE){ 
				l=i%N_Cells_x;
				m=(int) i/N_Cells_x;

				u_E = phi->val[EAST];
				u_W = phi->val[WEST];
				u_N = phi->val[NORTH];
				u_S = phi->val[SOUTH];
				u_P = phi->val[P];
			
				if(l==1)		temp->val[P] += 2.0*u_W - u_P;
				else 	 		temp->val[P] += u_W;
				if(l==N_Cells_x-2) 	temp->val[P] += 2.0*u_E - u_P;
				else 			temp->val[P] += u_E;
				if(m==1)		temp->val[P] += 2.0*u_S - u_P;
				else 			temp->val[P] += u_S;
				if(m==N_Cells_y-2)	temp->val[P] += 2.0*u_N - u_P;
				else 			temp->val[P] += u_N;
 
				temp->val[P] -= constant.f*pow(constant.h,2);
				temp->val[P] = temp->val[P]/4.0;
	
				e = temp->val[P] - phi->val[P];
				if(e > res)	res = e; 
		
				phi->val[P] = (1.0 - relax) * phi->val[P] + relax * temp->val[P];

			}
		}
	
		//Transferring values from temp to u
		for(i=0;i<N_Cells;i++){
			if(phi->bc[i] == NONE){
				phi->val[i] = temp->val[i];
			}
		}
		
		loop++;
	}

	printf("Maximum residual is %e and number of iterations are %ld\n", res, loop);	
	return;
}

void set_bc(Field *phi){

	
	int i, l, m;
	int N_Cells_x = phi->Nx;
	int N_Cells_y = phi->Ny;
	int N_Cells = N_Cells_x * N_Cells_y;
	
	for(i=0;i<N_Cells;i++){
		l = i%N_Cells_x;
		m = (int) i/N_Cells_x;
		
		if(phi->bc[i] == DIRICHLET){
			if(l==0) 		phi->val[i] = phi->bc_val[XMIN];
			else if(l==N_Cells_x-1)	phi->val[i] = phi->bc_val[XMAX];
			else if(m==0)		phi->val[i] = phi->bc_val[YMIN];
			else if(m==N_Cells_y-1) phi->val[i] = phi->bc_val[YMAX];
		}
	} return;

} 
		

int main(int argc, char **argv){


	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
	
	//Setting up virtual process grid topology - to be done immediately after MPI Initialization  
	Setup_proc_grid();	
	

	FILE *fp;

	//Defining a domain of type Domain
	Domain domain;

	//Defining a constant of type Constant
	Constant constant;
	
	//Defining the number of total interior control volumes
	int N_int_x = 100;
	int N_int_y = 100;

	//With ghost cells
	int N_Cells_x = N_int_x + 2;
	int N_Cells_y = N_int_y + 2;
	int N_Cells = N_Cells_x * N_Cells_y;

	//Defining the number of interior local control volumes 
	int N_int_local_x = N_int_x / P_grid[X_DIR]; 
	int N_int_local_y = N_int_y / P_grid[Y_DIR];	

	//With buffer cells
	int N_local_x = N_int_local_x + 2;
	int N_local_y = N_int_local_y + 2;
	int N_local = N_local_x * N_local_y;	

	//offset of cell numbering for the subdomain
	offset[X_DIR] = P_grid_coord[X_DIR] * N_int_local_x;
	offset[Y_DIR] = P_grid_coord[Y_DIR] * N_int_local_y;

	char filename[10];

	//Defining the values of the members of constant	
	constant.h = 0.1;
	constant.f = 0.0;

	//Allocating local memory to the members of the Field u	
	Field *u = allocate_field( N_local_x, N_local_y );
	
	
	int i, l, m;
	//Initializing all boundary condition values to zero
	for(i=0;i<5;i++){
		u->bc_val[i] = 0.0;
	}
	
	//Setting up boundary condition values
	u->bc_val[XMAX] = 1.0;
	u->bc_val[XMIN] = 1.0;
	u->bc_val[YMAX] = 1.0;
	u->bc_val[YMIN] = 1.0;

	//Equating the address of the field u to the member u of the domain
	domain.u = u;

	//Assigning different flags to different boundary conditions
	set_ghost_buffer_flag(domain);

//	set_buffer();	

	//initialization of local CVs in the process
	for(i=0;i<N_local;i++){
	//	l=i%N_Cells_x;
	//	m=(int) i/N_Cells_x; 
		u->val[i] = 0.0;
	}

	//Assigning the ghost cells the respective boundary face values
	set_bc(u);

	//Calling jacobi solver
	jacobi(u, N_local_x, N_local_y, constant);

	//calling gauss seidel solver
//	gauss_seidel(u, N_Cells_x, N_Cells_y, constant);

	sprintf(filename, "data%i", proc_rank);

	fp = fopen(filename, "w");
	
	for(i=0;i<N_local;i++){
		l=i%N_local_x;
		m=(int) i/N_local_x;
		
		fprintf(fp,"%lf %lf %lf\n", l*(constant.h), m*(constant.h), u->val[i]);
	}

	fclose(fp);	

	MPI_Finalize();

}


static Field *allocate_field(int N_x, int N_y){
	
	Field *phi;
	phi = new Field;
	phi->Nx = N_x;
	phi->Ny = N_y;
	phi->N = N_x*N_y;
	phi->val = new double [phi->N];
	phi->bc = new BC_type [phi->N];
	return phi;
}

 
