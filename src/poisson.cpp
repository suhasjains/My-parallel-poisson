#include "poisson.h"
namespace poisson {

/********************************************************
Solving a poisson equation in 2D parallely
Finite volume discretized equation used is UE + UW + UN + US - 4*UP = f*h^2

 h is the spacing. 
 f is the forcing function.
 Nx and Ny are the number of CVs in x and y direction respectively_
***********************************************************/


//MPI_Datatype exchange_buffer_type[2];                               //MPI_datatype for exchange of buffer cell data
//int P_grid_coord[2];                                                //coordinates of the process in the virtual grid
//int P_grid[2];                                                      //virtual grid dimensions

//Parallelization variables
int proc_rank;                                                      //rank of the current process
int P_grid_rank;                                                    //rank of the current proces in the virtual grid
int P_grid_top, P_grid_bottom, P_grid_left, P_grid_right;           //ranks of the neighbouring processes
int n_Procs;                                                        //total number of processes
int P_grid[2];                                                      //virtual grid dimensions
MPI_Datatype exchange_buffer_type[2];                               //MPI_datatype for exchange of buffer cell data
MPI_Comm grid_comm;                                                 //grid COMMUNICATOR
int offset[2];                                                      //offset for cell numbering for subdomains
int P_grid_coord[2];                                                //coordinates of the process in the virtual grid
MPI_Status status;

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
	int N = Nx * Ny;
	double global_res = 1.0;	

	long int loop = 0;

	//Defining a new temporary field (temp is not part of the domain)
	Field temp( Nx, Ny); 

	//Starting the iteration loop
	while(global_res > pow(10,-7)){	
	
		//making res 0 so that any error greater than 0 can be equated to this
		res = 0.0;

		//Making the temp field zero after every iteration
		temp.set_field_value(0.0);

		//exchanges buffer cells
		exchange_buffers(phi, Nx, Ny);	


		double u_E, u_W, u_N, u_S, u_P;
	
		for(i=0;i<N;i++){
			if(phi->bc[i] == NONE){ 
				l=i%Nx;
				m=(int) i/Nx;

				u_E = phi->val[EAST];
				u_W = phi->val[WEST];
				u_N = phi->val[NORTH];
				u_S = phi->val[SOUTH];
				u_P = phi->val[P];
			
				if(l==1&&P_grid_left==MPI_PROC_NULL)		temp.val[P] += 2.0*u_W - u_P;
				else 	 					temp.val[P] += u_W;
				if(l==Nx-2&&P_grid_right==MPI_PROC_NULL) 	temp.val[P] += 2.0*u_E - u_P;
				else 						temp.val[P] += u_E;
				if(m==1&&P_grid_bottom==MPI_PROC_NULL)		temp.val[P] += 2.0*u_S - u_P;
				else 						temp.val[P] += u_S;
				if(m==Ny-2&&P_grid_top==MPI_PROC_NULL)		temp.val[P] += 2.0*u_N - u_P;
				else 						temp.val[P] += u_N;
 
				temp.val[P] -= constant.f*pow(constant.h,2);
				temp.val[P] = temp.val[P]/4.0;

	
				e = temp.val[P] - phi->val[P];
				if(e > res)	res = e;

	 

			}
		}
	
		//Transferring values from temp to u
		for(i=0;i<N;i++){
			if(phi->bc[i] == NONE){
				phi->val[i] = temp.val[i];
			}
		}


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
        int N = Nx * Ny;
        double global_res = 1.0;
	
	double relax = 1.5;
	double prev_val, curr_val;

        long int loop = 0;

        //Defining a new temporary field (temp is not part of the domain)
        Field temp( Nx, Ny);

        //Starting the iteration loop
        while(global_res > pow(10,-10)){

                //making res 0 so that any error greater than 0 can be equated to this
                res = 0.0;

                //Making the temp field zero after every iteration
                temp.set_field_value(0.0);

                //exchanges buffer cells
                exchange_buffers(phi, Nx, Ny);


                double u_E, u_W, u_N, u_S, u_P;

                for(i=0;i<N;i++){
                        if(phi->bc[i] == NONE && i%2 == 0 ){
                                l=i%Nx;
                                m=(int) i/Nx;

                                u_E = phi->val[EAST];
                                u_W = phi->val[WEST];
                                u_N = phi->val[NORTH];
                                u_S = phi->val[SOUTH];
                                u_P = phi->val[P];

                                if(l==1&&P_grid_left==MPI_PROC_NULL)            temp.val[P] += 2.0*u_W - u_P;
                                else                                            temp.val[P] += u_W;
                                if(l==Nx-2&&P_grid_right==MPI_PROC_NULL)        temp.val[P] += 2.0*u_E - u_P;
                                else                                            temp.val[P] += u_E;
                                if(m==1&&P_grid_bottom==MPI_PROC_NULL)          temp.val[P] += 2.0*u_S - u_P;
                                else                                            temp.val[P] += u_S;
                                if(m==Ny-2&&P_grid_top==MPI_PROC_NULL)          temp.val[P] += 2.0*u_N - u_P;
				else                                            temp.val[P] += u_N;

                                temp.val[P] -= constant.f*pow(constant.h,2);
                                temp.val[P] = temp.val[P]/4.0;

				prev_val = phi->val[P];
				
				//SOR
				phi->val[P] = (1.0 - relax) * phi->val[P] + relax * temp.val[P];

				curr_val = phi->val[P];

    	                   	e = curr_val - prev_val;

                                if(e > res)     res = e;

                        }
                }

                //exchanges buffer cells
                exchange_buffers(phi, Nx, Ny);

                for(i=0;i<N;i++){
                        if(phi->bc[i] == NONE && i%2 != 0 ){
                                l=i%Nx;
                                m=(int) i/Nx;

                                u_E = phi->val[EAST];
                                u_W = phi->val[WEST];
                                u_N = phi->val[NORTH];
                                u_S = phi->val[SOUTH];
                                u_P = phi->val[P];

                                if(l==1&&P_grid_left==MPI_PROC_NULL)            temp.val[P] += 2.0*u_W - u_P;
                                else                                            temp.val[P] += u_W;
                                if(l==Nx-2&&P_grid_right==MPI_PROC_NULL)        temp.val[P] += 2.0*u_E - u_P;
                                else                                            temp.val[P] += u_E;
                                if(m==1&&P_grid_bottom==MPI_PROC_NULL)          temp.val[P] += 2.0*u_S - u_P;
                                else                                            temp.val[P] += u_S;
                                if(m==Ny-2&&P_grid_top==MPI_PROC_NULL)          temp.val[P] += 2.0*u_N - u_P;
				else                                            temp.val[P] += u_N;

                                temp.val[P] -= constant.f*pow(constant.h,2);
                                temp.val[P] = temp.val[P]/4.0;

				prev_val = phi->val[P];
				
				//SOR
				phi->val[P] = (1.0 - relax) * phi->val[P] + relax * temp.val[P];

				curr_val = phi->val[P];

    	                   	e = curr_val - prev_val;

                                if(e > res)     res = e;

                        }
                }
				//printf("My residual = %g and I am %d\n",res,proc_rank); 

                //if(loop%10==0)  MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MAX, grid_comm);
		MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MAX, grid_comm);
                loop++;
        }

        printf("Maximum residual is %e and number of iterations are %ld and I am process %d \n", res, loop, proc_rank);

        return;
}


//Steepest gradient solver
void compute_AX(Field* phi, Field* temp, Constant constant) {

	int i, l, m;
	int Nx = phi->Nx;
	int Ny = phi->Ny; 
	int N = Nx*Ny;	
        double u_E, u_W, u_N, u_S, u_P;

//	for(i = 0 ; i < N ; i++)
//		temp[i] = 0.0;

	temp->set_field_value(0.0);
	
	//exchanges buffer cells
  	exchange_buffers(phi, Nx, Ny);
	
	for(i=0;i<N;i++){
		if(phi->bc[i] == NONE){ 
			l=i%Nx;
			m=(int) i/Nx;

			u_E = phi->val[EAST];
			u_W = phi->val[WEST];
			u_N = phi->val[NORTH];
			u_S = phi->val[SOUTH];
			u_P = phi->val[P];
		
			if(l==1&&P_grid_left==MPI_PROC_NULL)		temp->val[P] += 2.0*u_W - u_P;
			else 	 					temp->val[P] += u_W;
			if(l==Nx-2&&P_grid_right==MPI_PROC_NULL) 	temp->val[P] += 2.0*u_E - u_P;
			else 						temp->val[P] += u_E;
			if(m==1&&P_grid_bottom==MPI_PROC_NULL)		temp->val[P] += 2.0*u_S - u_P;
			else 						temp->val[P] += u_S;
			if(m==Ny-2&&P_grid_top==MPI_PROC_NULL)		temp->val[P] += 2.0*u_N - u_P;
			else 						temp->val[P] += u_N;

			temp->val[P] -= 4.0*u_P;
			temp->val[P] /= pow(constant.h,2);


		}
		
//		else temp->val[P] = phi->val[P];
		//else printf("Yes\n");

	}

	return;
}

//Steepest descent solver
void solve_SD(Field* phi, Constant constant) {

	long int loop = 0;
//	double *res, *temp, *x, *xo;
	double del, alpha, global_del;
	double rTAr, global_rTAr;
	int i, l, m;
	int Nx = phi->Nx;
	int Ny = phi->Ny; 
	int N = Nx*Ny;	

//	res = new double[N];
//	temp = new double[N];
//	x = new double[N];

	//Field temp(Nx, Ny);
	//Field res(Nx, Ny);
	//Field x(Nx, Ny);


	Field res = *phi;
	Field temp = *phi;
	Field x = *phi;
	
	res.set_field_value(0.0);
	x.set_field_value(0.0);
	temp.set_field_value(0.0);


	//calculates Ax
	compute_AX(phi,&temp,constant);


	//calculates residual and residual norm
	del = 0.0;
	for(i = 0 ; i < N ; i++) {
	//	if(res.bc[i]!=NONE) {
	//		printf("f = %e temp.val[i] = %e\n",constant.f,temp.val[i]);
	//	}
		if(res.bc[i]==NONE) {
			res.val[i] = constant.f - temp.val[i];
			del += res.val[i]*res.val[i]; 		
		}
	}
        //printf("del = %e\n", del);
	//printf("Maximum residual norm is %e and number of iterations are %ld and I am process %d \n", del, loop, proc_rank);

	//get the global summation of the norm
	MPI_Allreduce(&del, &global_del, 1, MPI_DOUBLE, MPI_SUM, grid_comm);
	

	//initial vector
	for(i = 0 ; i < N ; i++) {
		if(phi->bc[i]==NONE) {
			x.val[i] =  phi->val[i];  
		}
	}	

	while (global_del > 10E-10) {
	//while (loop < 3000) {

		//calculates Ar
		compute_AX(&res,&temp,constant);
			
		//calculates alpha
		rTAr = 0.0;	
		for(i = 0 ; i < N ; i++) {
			if(res.bc[i]==NONE) {
				rTAr += res.val[i]*temp.val[i];
			}
		}
	
		//get the global summation of the norm
		MPI_Allreduce(&rTAr, &global_rTAr, 1, MPI_DOUBLE, MPI_SUM, grid_comm);

		alpha = global_del/global_rTAr;

		//update vector
		for(i = 0 ; i < N ; i++) {
			if(res.bc[i]==NONE) {	
				x.val[i] = x.val[i] + alpha*res.val[i]; 
			}
		}


		if(loop%50==0) compute_AX(&x, &temp, constant);
	
		del = 0.0;
		//update residual and residual norm	
		for(i = 0 ; i < N ; i++) {
			if(res.bc[i]==NONE) {
				
				if(loop%50==0) 
					res.val[i] = constant.f - temp.val[i];

				else res.val[i] = res.val[i] - alpha*temp.val[i]; 
				
				del += res.val[i]*res.val[i]; 		
			}
		}

        	//printf("Maximum residual norm is %e and number of iterations are %ld and I am process %d \n", del, loop, proc_rank);
		
		//get the global summation of the norm
		MPI_Allreduce(&del, &global_del, 1, MPI_DOUBLE, MPI_SUM, grid_comm);
		
		
		loop++;

	}	
	
	//final vector
	for(i = 0 ; i < N ; i++) {
		if(phi->bc[i]==NONE) {
			phi->val[i] = x.val[i];  
		}
	}

        printf("Maximum residual norm is %e and number of iterations are %ld and I am process %d \n", global_del, loop, proc_rank);
		

	return;
}


//Conjugate gradient solver


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


        MPI_Init(&argc, &argv);


	double t1,t2;
	int i, l, m;
        t1 = MPI_Wtime();
	
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
	
	//Setting up virtual process grid topology - to be done immediately after MPI Initialization  
	setup_proc_grid();	
	
	
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

	//Defining a constant of type Constant
	Constant constant;

	//Defining the values of the members of constant	
	constant.h = 0.05;
	constant.f = 1.0;

	//Defining a new scalar field	
	Field u( N_local_x, N_local_y); 
	
	//Defining a domain of type Domain
	Domain domain(&u,&constant);
	
	
	//Initializing all boundary condition values to zero
	for(i=0;i<5;i++){
		u.bc_val[i] = 0.0;
	}
	
	//Setting up boundary condition values
	u.bc_val[XMAX] = 1.0;
	u.bc_val[XMIN] = 2.0;
	u.bc_val[YMAX] = 3.0;
	u.bc_val[YMIN] = 4.0;


	//Assigning different flags to different boundary conditions
	set_ghost_buffer_flag(domain);

	//Setting up MPI datatypes for exchange of values in buffer in x and y direction.
	setup_MPI_datatypes(N_int_local_x, N_int_local_y, N_local_x);

	//initialization of local CVs in the process
	u.set_field_value(0.0);

	//Assigning the ghost cells the respective boundary face values
	set_bc(&u);

	//Calling jacobi solver
	//jacobi(&u, N_local_x, N_local_y, constant);
	//gauss_seidel(&u, N_local_x, N_local_y, constant);
	solve_SD(&u, constant);

	//writes output to the file
	write_output(domain,proc_rank);


	t2 = MPI_Wtime();

        printf("Total time taken is %g\n",t2-t1);

	MPI_Finalize();

}


 
