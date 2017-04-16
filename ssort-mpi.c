/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include "util.h"

static int compare(const void *a, const void *b)
{
  int *da = (int *)a;
  int *db = (int *)b;

  if (*da > *db)
    return 1;
  else if (*da < *db)
    return -1;
  else
    return 0;
}

int main( int argc, char *argv[])
{
  int rank;
  int i, j, N, P, S, SN, newN, d;
  double T1, T2;
  int *Samples = NULL;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &P);


  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */
   
  MPI_Barrier(MPI_COMM_WORLD);
  timestamp_type time1, time2;
  get_timestamp(&time1);

  sscanf(argv[1], "%d", &N);   // numbers on each tasks.

  int *vec = calloc(N, sizeof(int));
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  
  for (i = 0; i < N; ++i) {
    vec[i] = rand();
  }
  //printf("rank: %d, first entry: %d\n", rank, vec[0]);

  /* sort locally */
  qsort(vec, N, sizeof(int), compare);


  /* randomly sample s entries from vector or select local splitters,
   * i.e., every N/P-th entry of the sorted vector,  */
   S = P;
   d = (int) N/S;    /* every d-th entry of the sorted vector */
   int *Svec = calloc(S, sizeof(int));     //The sample vector has p entries.
   for (i = 0; i < S; i++) {
      Svec[i] = vec[i*d];             // For example, N=12, S=3, choose  vec[0],vec[4],vec[8].
  }

  /* every processor communicates the selected entries
   * to the root processor; use for instance an MPI_Gather */
  if(rank == 0){
    printf("N: %d, P: %d\n",N,P);
    SN = S*P;    // total sample number
    Samples = calloc(SN, sizeof(int));
  }
  

  MPI_Gather(Svec, S, MPI_INT, Samples, S, MPI_INT, 0, MPI_COMM_WORLD);   // Gather the samples to root.

  /* root processor does a sort, determinates splitters that
   * split the data into P buckets of approximately the same size,
   * need P-1 splitters*/
  
  int *splitter = calloc(P-1, sizeof(int));
  if (rank == 0) {
  qsort(Samples, SN, sizeof(int), compare);

  for (i = 1; i < P; i++) {
      splitter[i-1] = Samples[i*S+1];   // SN/P=S
      }
  }

  /* root process broadcasts splitters */
  MPI_Bcast(splitter, P-1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  /* every processor uses the obtained splitters to decide
   * which integers need to be sent to which other processor (local bins) 
   * number of entries and the displacement of the first entry in each bucket need to be sent later.*/

  int *sendcounts = calloc(P, sizeof(int));
  for (i = 0; i < N; i++) {                      // Compare vec[i] with each splitter.
  	if (vec[i] <= splitter[0]) {
  			sendcounts[0]++; }
  	else if (vec[i] > splitter[P-2]) {
  			sendcounts[P-1]++; }	
  	else{
  	for (j = 0; j < P-2; j++) {
		if (vec[i] <= splitter[j+1] && vec[i] > splitter[j]) {
			sendcounts[j+1]++; break; }
	        }
		}
  }
  
  int *sdispls = calloc(P, sizeof(int));
  for (i = 1; i < P; i++){
      sdispls[i] = sdispls[i-1] + sendcounts[i-1];
  }
  

  //pArray(sdispls,P,rank);
  /* send and receive: either you use MPI_AlltoallV, or
   * (and that might be easier), use an MPI_Alltoall to share
   * with every processor how many integers it should expect,
   * and then use MPI_Send and MPI_Recv to exchange the data */
   int *recvcounts = calloc(P, sizeof(int));
   MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);
   
  /* set receive displacements */
  int *rdispls = calloc(P, sizeof(int));
  for (i = 1; i < P; ++i){
      rdispls[i] = rdispls[i-1] + recvcounts[i-1];
  }

  newN = recvcounts[P-1] + rdispls[P-1];  //size of new vector
  int *newvec = calloc(newN, sizeof(int)); // receive buffer

  MPI_Alltoallv(vec, sendcounts, sdispls, MPI_INT, newvec, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);
  
  /* do a local sort */
  qsort(newvec, newN, sizeof(int), compare); 
  /* every processor writes its result to a file */
  /* timing */

  MPI_Barrier(MPI_COMM_WORLD);
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);

  if (0 == rank) {
    printf("Time elapsed is %f seconds.\n", elapsed);
  }

  {
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", rank);
    fd = fopen(filename,"w+");

    if(NULL == fd)
    {
      printf("Error opening file \n");
      return 1;
    }
    for(i = 0; i < newN; i++)
      fprintf(fd, "%d\n", newvec[i]);

    fclose(fd);
  }
  
  
  free(vec);
  free(Svec);
  free(splitter);
  free(Samples);
  free(sendcounts);
  free(sdispls);
  free(recvcounts);
  free(rdispls);
  free(newvec);
  
  MPI_Finalize();
  return 0;
}