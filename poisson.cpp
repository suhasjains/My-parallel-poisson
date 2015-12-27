#include "poisson.h"
namespace poisson {

/********************************************************
Solving a poisson equation in 2D parallely
Finite volume discretized equation used is UE + UW + UN + US - 4*UP = f*h^2

 h is the spacing. 
 f is the forcing function.
 Nx and Ny are the number of CVs in x and y direction respectively_
***********************************************************/



//Setting flags for ghost and buffer cells
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

//Jacobi solver
void jacobi(Field *phi, int Nx, int Ny, Constant constant){

	double res, e;	
	int i, l, m;
	int N_Cells_x = Nx;
	int N_Cells_y = Ny;
	int N = Nx * Ny;
	int N_Cells = N;
	double global_res = 1.0;	

	long int loop = 0;

	//Defining a new temporary field (temp is not part of the domain)
	Field tempo( N_Cells_x, N_Cells_y); 
	Field *temp;
	temp = &tempo;

	//making res 1 so that it enters into the loop	
	res = 1.0;
	
	//Starting the iteration loop
	//for(t=0;t<10000;t++){
	while(global_res > pow(10,-4)){	
//	while(loop<45000){
	//making res 0 so that any error greater than 0 can be equated to this
	res = 0.0;

		//Making the temp field zero after every iteration
		temp->set_field_value(0.0);

		exchange_buffers(phi, Nx, Ny);	


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
	
//		Exchange_buffers(phi, Nx, Ny);	

//		printf("Done with iteration %ld\n",loop);

	//	printf("%ld iterations done and residual is %lf\n",loop, res);	

		if(loop%10==0)	MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MAX, grid_comm);
	
		loop++;
	}

	printf("Maximum residual is %e and number of iterations are %ld and I am process %d \n", res, loop, proc_rank);	
	
	return;
}


//Gauss seidel solver
void gauss_seidel(Field *phi, int Nx, int Ny, Constant constant){

	double res, e;	
	int i, l, m;
	long int loop = 0;
	int N_Cells_x = Nx;
	int N_Cells_y = Ny;
	int N = Nx * Ny;
	int N_Cells = N;

	double relax = 1.95;

	//Defining a new temporary field (temp is not part of the domain)
	Field tempo( N_Cells_x, N_Cells_y); 
	Field *temp;
	temp = &tempo;

	//making res 1 so that it enters into the loop	
	res = 1.0;
	
	//Starting the iteration loop
	//for(t=0;t<10000;t++){
	while(res > pow(10,-10)){	
	
	//making res 0 so that any error greated than 0 can be equated to this
	res = 0.0;

		//Making the temp field zero after every iteration
		temp->set_field_value(0.0);

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


//Setting boundary condition values to boundary cells
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

}



using namespace poisson;


int main(int argc, char **argv){

	double t1,t2;

        MPI_Init(&argc, &argv);

        t1 = MPI_Wtime();
	
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
	
	//Setting up virtual process grid topology - to be done immediately after MPI Initialization  
	setup_proc_grid();	
	
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


	//Defining the values of the members of constant	
	constant.h = 0.1;
	constant.f = 1.0;

	//Defining a new scalar field	
	Field u( N_local_x, N_local_y); 
	
	//Defining a domain of type Domain
	Domain domain(&u,&constant);
	
	
	//Equating the address of the field u to the member u of the domain
	//domain.u = &u;
	//domain.constant = &constant;
	
	int i, l, m;
	//Initializing all boundary condition values to zero
	for(i=0;i<5;i++){
		u.bc_val[i] = 0.0;
	}
	
	//Setting up boundary condition values
	u.bc_val[XMAX] = 1.0;
	u.bc_val[XMIN] = 1.0;
	u.bc_val[YMAX] = 1.0;
	u.bc_val[YMIN] = 1.0;


	//Assigning different flags to different boundary conditions
	set_ghost_buffer_flag(domain);

	//Setting up MPI datatypes for exchange of values in buffer in x and y direction.
	setup_MPI_datatypes(N_int_local_x, N_int_local_y, N_local_x);

//	set_buffer();	

	//initialization of local CVs in the process
	u.set_field_value(0.0);

	//Assigning the ghost cells the respective boundary face values
	set_bc(&u);

	//Calling jacobi solver
	jacobi(&u, N_local_x, N_local_y, constant);


	write_output(domain,proc_rank);


	t2 = MPI_Wtime();

        printf("Total time taken is %g\n",t2-t1);

	MPI_Finalize();

}


 
