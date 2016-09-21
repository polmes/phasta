#include<phasta.h>
#include<mpi.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include <stdio.h>


int 
main( int argc,   
      char* argv[] ) {

  struct rlimit lim;
  memset(&lim, 0, sizeof(struct rlimit));
  getrlimit(RLIMIT_STACK, &lim);
  printf("stack: soft: %d, hard: %d, inf: %d\n", lim.rlim_cur, lim.rlim_max, RLIM_INFINITY);


  MPI_Init(&argc,&argv);
  phasta ( argc, argv);  
  MPI_Finalize();
  return 0;
}
