#ifndef _PARALLEL_POISSON_H
#define _PARALLEL_POISSON_H
#include <stdio.h>
#include <cmath>
#include "variables.h"
#include "parallel.cpp"
namespace poisson {

//defines a new field variable
//static Field *allocate_field(int,int);


//sets flags for ghost and buffer cells
void set_ghost_buffer_flag(Domain);

//Jacobi solver
void jacobi(Field, int, int, Constant);

//Gauss seidel solver 
void gauss_seidel(Field, int, int, Constant);

//Sets boundary condition values to boundary cells
void set_bc(Field);


}
#endif
