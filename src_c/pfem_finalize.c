/***
 *** PFEM_FINALIZE
 ***/
#include <stdio.h>
#include <stdlib.h>
#include "pfem_util.h"
void  PFEM_FINALIZE()
{
     MPI_Finalize ();

     if( my_rank == 0 ){
       fprintf(stdout,"* normal termination\n");
       exit(0);
     }
}
