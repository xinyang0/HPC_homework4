/* MPI-parallel Jacobi smoothing to solve - Laplace u=f
 * Matrix size N*N. p=4^j processor. Each processor has Nl*Nl size chunk data (Nl=N/sqrt(p)).
 * Author: Xinyang Wang
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "util.h"
#include <string.h>


/* compuate global residual, assuming ghost values are updated */
double compute_residual(double **lu, int lN, double invhsq)
{
  int i, j;
  double tmp, gres = 0.0, lres = 0.0;
 for (i = 1; i <= lN; i++) {
  	for (j = 1; j <= lN; j++) {
    	tmp = ((4.0*lu[i][j] - lu[i-1][j] - lu[i][j-1]- lu[i+1][j] - lu[i][j+1]) * invhsq - 1);
    	lres += tmp * tmp;
  	}
  }
  
  /* use allreduce for convenience; a reduce would also be sufficient */
  MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sqrt(gres);
}


int main(int argc, char * argv[])
{
  int mpirank, i, j, p, q, N, lN, iter, max_iters;
  MPI_Status status, status1, status2, status3;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  
  /* get name of host running MPI process */
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);

  sscanf(argv[1], "%d", &N);
  sscanf(argv[2], "%d", &max_iters);

  /* compute number of unknowns handled by each process */
  q = (int)sqrt(p);
  lN = N / q;
  if ((N % q != 0) && mpirank == 0 ) {
    printf("N: %d, local N: %d\n", N, lN);
    printf("Exiting. N must be a multiple of sqrt(p)\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  
  
  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  timestamp_type time1, time2;
  get_timestamp(&time1);

  /* Allocation of vectors, including left/upper and right/lower ghost points */
  
  double ** lu    = (double **) calloc(sizeof(double *), lN+2);
  for (i = 0; i < lN+2; i++) {
  	lu[i] = (double *) calloc(sizeof(double), lN+2);
  	}
  
  double ** lunew = (double **) calloc(sizeof(double *), lN+2);
  for (i = 0; i < lN+2; i++) {
  	lunew[i] = (double *) calloc(sizeof(double), lN+2);
  	}

  double ** lutemp;
  double * send_left  = (double *) calloc(sizeof(double), lN); 
  double * send_right = (double *) calloc(sizeof(double), lN);  
  double * ghost_left  = (double *) calloc(sizeof(double), lN); 
  double * ghost_right = (double *) calloc(sizeof(double), lN);

  

  double h = 1.0 / (N + 1);
  double hsq = h * h;
  double invhsq = 1./hsq;
  double gres, gres0, tol = 1e-5;

  /* initial residual */
  gres0 = compute_residual(lu, lN, invhsq);
  gres = gres0;

  for (iter = 0; iter < max_iters && gres/gres0 > tol; iter++) {

    /* Jacobi step for local points */
    
    for (i = 1; i <= lN; i++) {
    	for (j = 1; j <= lN; j++){
      		lunew[i][j]  = 0.25 * (hsq + lu[i-1][j] + lu[i][j-1] + lu[i+1][j] + lu[i][j+1]);
      }
    }


    for (i = 1; i <= lN; i++){
      send_left[i-1]  = lunew[i][1];
    }
    for (i = 1; i <= lN; i++){
      send_right[i-1]  = lunew[i][lN];
    }

    
    /* communicate ghost values */
    /* If the chunk is not on the right side, send right boundary values to the right */
    /* and receive ghost right boundary values from the right.*/
    if (mpirank%q < q - 1) {
      MPI_Send(send_right, lN, MPI_DOUBLE, mpirank+1, 11, MPI_COMM_WORLD);
      MPI_Recv(ghost_right, lN, MPI_DOUBLE, mpirank+1, 22, MPI_COMM_WORLD, &status);
    }
    
    /* If the chunk is not on the left side, send left boundary values to the left */
    /* and receive ghost left boundary values from the left.*/
    if (mpirank%q > 0) {
      MPI_Recv(ghost_left, lN, MPI_DOUBLE, mpirank-1, 11, MPI_COMM_WORLD, &status1);
      MPI_Send(send_left, lN, MPI_DOUBLE, mpirank-1, 22, MPI_COMM_WORLD);
    }
    
    //In case
    //MPI_Barrier(MPI_COMM_WORLD);

    /* If the chunk is not on the bottom side, send top boundary values to the lower chunk */
    /* and receive ghost bottom boundary values from the lower chunk.*/
     if (mpirank >= q) {
      MPI_Send(lunew[1], lN+2, MPI_DOUBLE, mpirank-q, 33, MPI_COMM_WORLD);
      MPI_Recv(lunew[0], lN+2, MPI_DOUBLE, mpirank-q, 44, MPI_COMM_WORLD, &status2);
    }
    
    /* If the chunk is not on the top side, send top boundary values to the upper chunk */
    /* and receive ghost top boundary values from the upper chunk.*/
    if (mpirank < p - q) {
      MPI_Recv(lunew[lN+1], lN+2, MPI_DOUBLE, mpirank+q, 33, MPI_COMM_WORLD, &status3);
      MPI_Send(lunew[lN], lN+2, MPI_DOUBLE, mpirank+q, 44, MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (i = 1; i <= lN; i++){
      lunew[i][0] = ghost_left[i-1];
    }
    for (i = 1; i <= lN; i++){
     lunew[i][lN+1] = ghost_right[i-1];
    }

    /* copy newu to u using pointer flipping */
    lutemp = lu; lu = lunew; lunew = lutemp;
   if (0 == (iter % 10)) {
      gres = compute_residual(lu, lN, invhsq);
      if (0 == mpirank) {
	printf("Iter %d: Residual: %g\n", iter, gres);
      }
    }
  }
   

    /* Clean up */
  	for (i = 0; i < lN+2; i++) {
  		free(lu[i]);
  	}
  	free(lu);
    
    for (i = 0; i < lN+2; i++) {
  		free(lunew[i]);
  	}
  	free(lunew);
  	
  	free(send_left);
  	free(send_right);

  	
  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  if (0 == mpirank) {
    printf("Time elapsed is %f seconds.\n", elapsed);
  }
  MPI_Finalize();
  return 0;
}
