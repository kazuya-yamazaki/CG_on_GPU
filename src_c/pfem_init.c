/**
 ** PFEM_INIT
**/
#include "pfem_util.h"

#ifdef OPENACC
#include <openacc.h>
#endif
/**
   INIT. PFEM-FEM process's
**/
void PFEM_INIT(int argc, char* argv[])
{
  int i, ngpu;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&PETOT);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  
#ifdef OPENACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  acc_set_device_num(my_rank % ngpu, acc_device_nvidia);
#endif
  
  for(i=0;i<100;i++)  pfemRarray[i]=0.0;
  for(i=0;i<100;i++)  pfemIarray[i]=0;
}
