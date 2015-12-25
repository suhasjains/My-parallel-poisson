#ifndef VARIABLES_H
#define VARIABLES_H
#include<string.h>	//for memcpy

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

	//Members
	public:
        	int Nx,Ny;
        	int N;
        	double *val;
        	BC_type *bc;
        	double bc_val[5];
	
	//Constructor
	Field(int N_x, int N_y);

	//Copy constructor
	Field(const Field &obj);
		
	//Destructor
	~Field();
};


//Constructor for Field class
Field :: Field( int N_x, int N_y ) : Nx(N_x), Ny(N_y) {

        N = N_x*N_y;
        val = new double [N];
        bc = new BC_type [N];
}

//Copy constructor for Field class
Field :: Field(const Field &obj) {

	Nx = obj.Nx;
	Ny = obj.Ny;
	N = obj.N;
        val = new double [N];
	memcpy(val,obj.val,sizeof(double)*N);
        bc = new BC_type [N];
	memcpy(bc,obj.bc,sizeof(double)*N);
}

//Destructor for Field class
Field :: ~Field() {

        delete val;
        delete  bc;
}

//Patch
typedef enum{
        INSIDE,
        XMIN,
        XMAX,
        YMIN,
        YMAX
} PATCH_type;
PATCH_type patch;

//Domain structure 
class Domain{

	public:
        	Field *u;
};


//Constants
class Constant{

	public:
        	double h;
        	double f;
};

}
#endif
