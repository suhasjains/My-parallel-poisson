//#include "output.h"
#include "poisson.h"
namespace poisson {

void write_output(Domain domain, int proc_rank) {

        int i,l,m;
        int N_local = domain.u->N;
        int N_local_x = domain.u->Nx;
	int N_local_y = domain.u->Ny;
        FILE *fp;

        char filename[10];

	//offset of cell numbering for the subdomain
        offset[X_DIR] = P_grid_coord[X_DIR] * (N_local_x-2);
        offset[Y_DIR] = P_grid_coord[Y_DIR] * (N_local_y-2);


        sprintf(filename, "data%d", proc_rank);


        fp = fopen(filename, "w");


        for(i=0;i<N_local;i++){
               // if(domain.u->bc[i] == NONE){
                        l=i%N_local_x;
                        m=(int) i/N_local_x;

                        fprintf(fp,"%lf %lf %lf\n", (l+offset[X_DIR])*(domain.constant->h), (m+offset[Y_DIR])*(domain.constant->h), domain.u->val[i]);
                //}
        }

        fclose(fp);


}

}
