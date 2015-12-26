#ifndef VARIABLES_H
#define VARIABLES_H
#include<string.h>	//for memcpy
#include<string>	//for strings

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
        	int Nx,Ny;		//size
        	int N;			//size
        	double *val;		//values
        	BC_type *bc;		//BC type
        	double bc_val[5];	//BC value

	//Member functions
	void set_field_value(double value) {

		for(int i=0;i<N;i++){
                        this->val[i] = value;
                }
	}

	
	//Constructors
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
PATCH_type patch;


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
        	Field *u;
		Constant *constant;
	private:
	protected:
};



}
#endif
