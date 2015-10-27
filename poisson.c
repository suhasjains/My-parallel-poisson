#include <stdio.h>
#include <cmath>

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


//defining enums and structures
typedef enum{
	NONE,
	DIRICHLET,
} BC_type;


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
void set_ghosts(Domain);
//void set_bc(Field);
static Field *allocate_field(int,int);
//void jacobi1(Field , int, int, Constant);


void set_ghosts(Domain domain){
	
	int i, l, m;

	Field *u = domain.u;

	int N_Cells = u->Nx * u->Ny;
	int N_Cells_x = u->Nx;
	int N_Cells_y = u->Ny;

	for(i=0;i<N_Cells;i++){
		l = i%N_Cells_x;	//gives the position of the CV along x
		m = (int) i/N_Cells_x;	//gives the position of the CV along y
		if(l==0 || l==N_Cells_x-1 || m==0 || m==N_Cells_y-1){
			u->bc[i] = DIRICHLET;
		}else{
			u->bc[i] = NONE;
		}
	}return;
}
void jacobi1(Field *phi, int Nx, int Ny, Constant constant){

	double res;	
	int i, l, m, t;
	int N_Cells_x = Nx;
	int N_Cells_y = Ny;
	int N = Nx * Ny;
	int N_Cells = N;

	//Allocating memory to the members of the temporary Field temp (temp is not part of the domain)
	Field *temp = allocate_field( N_Cells_x, N_Cells_y );
	
	//Starting the iteration loop
	for(t=0;t<100000;t++){
		
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
			}
		}
	
		//Transferring values from temp to u
		for(i=0;i<N_Cells;i++){
			if(phi->bc[i] == NONE){
				phi->val[i] = temp->val[i];
			}
		}
	}return;


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
		

int main(){

	FILE *fp;

	//Defining a domain of type Domain
	Domain domain;

	//Defining a constant of type Constant
	Constant constant;
	
	//Defining the number of control volumes which includes the ghost cells
	int Nx,Ny,N;
	int N_Cells_x = Nx = 100+2;
	int N_Cells_y = Ny = 100+2;
	int N_Cells = N = N_Cells_x * N_Cells_y;

	//Defining the values of the members of constant	
	constant.h = 0.1;
	constant.f = 0.0;

	//Allocating memory to the members of the Field u	
	Field *u = allocate_field( N_Cells_x, N_Cells_y );
	
	
	int i, l, m;
	double res;
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
	set_ghosts(domain);


	//initialization
	for(i=0;i<N_Cells;i++){
	//	l=i%N_Cells_x;
	//	m=(int) i/N_Cells_x; 
		u->val[i] = 0.0;
	}

	//Assigning the ghost cells the respective boundary face values
	set_bc(u);

	//Calling jacobi solver
	jacobi1(u, N_Cells_x, N_Cells_y, constant);
	

	fp = fopen("data", "w");
	
	for(i=0;i<N_Cells;i++){
		l=i%N_Cells_x;
		m=(int) i/N_Cells_x;
		
		fprintf(fp,"%lf %lf %lf\n", l*(constant.h), m*(constant.h), u->val[i]);
	}

	fclose(fp);	


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

 
