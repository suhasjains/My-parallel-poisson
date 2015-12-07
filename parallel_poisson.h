#include <stdio.h>
#include <cmath>
//#include "parallel.h"

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


//Field structure
typedef struct _Field Field;
struct _Field{
        int Nx,Ny;
        int N;
        double *val;
        BC_type *bc;
        double bc_val[5];
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


//Domain structure 
typedef struct _Domain Domain;
struct _Domain{
        Field *u;
};


//Constants
typedef struct _Constant Constant;
struct _Constant{
        double h;
        double f;
};


//defines a new field variable
static Field *allocate_field(int,int);


