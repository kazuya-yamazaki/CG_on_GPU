#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "precision.h"
#include "allocate.h"
static MPI_Status  *sta;
static MPI_Request *req;
static KINT NFLAG=0;
extern FILE *fp_log;
void SOLVER_SEND_RECV ( int N0, int N, int NEIBPETOT,
			int NEIBPE[], int IMPORT_INDEX[], int IMPORT_ITEM[],
			int EXPORT_INDEX[], int EXPORT_ITEM[],
			KREAL WS[], KREAL WR[], KREAL X[], int my_rank)
{

  int  ii,k,neib,istart,inum;
/***
    INIT.
***/
  if( NFLAG == 0 ){
    sta=(MPI_Status*)allocate_vector(sizeof(MPI_Status),NEIBPETOT*2);
    req=(MPI_Request*)allocate_vector(sizeof(MPI_Request),NEIBPETOT*2);
    NFLAG=1;
  } 

#pragma acc data present(IMPORT_INDEX[:NEIBPETOT+1], IMPORT_ITEM[:IMPORT_INDEX[NEIBPETOT]]) \
                 present(EXPORT_INDEX[:NEIBPETOT+1], EXPORT_ITEM[:EXPORT_INDEX[NEIBPETOT]]) \
                 present(WS[:N], WR[:N], X[:N])
{

/***
    SEND
***/
  for( neib=1;neib<=NEIBPETOT;neib++){
    istart=EXPORT_INDEX[neib-1];
    inum  =EXPORT_INDEX[neib]-istart;
#pragma omp parallel for private (k,ii)
#pragma acc parallel loop
    for( k=istart;k<istart+inum;k++){
      ii= EXPORT_ITEM[k];
      WS[k]= X[ii-1];
    }
    
#pragma acc host_data use_device(WS)
    MPI_Isend(&WS[istart],inum,MPI_DOUBLE,
	      NEIBPE[neib-1],0,MPI_COMM_WORLD,&req[neib-1]);
  }
/***
    RECEIVE
***/
  for( neib=1;neib<=NEIBPETOT;neib++){
    istart = IMPORT_INDEX[neib-1] + N0;
    inum   = IMPORT_INDEX[neib] - IMPORT_INDEX[neib-1];
#pragma acc host_data use_device(X)
    MPI_Irecv(&X[istart],inum,MPI_DOUBLE,
	      NEIBPE[neib-1],0,MPI_COMM_WORLD,&req[NEIBPETOT+neib-1]);
  }
  
  MPI_Waitall (2*NEIBPETOT, req, sta);

} /* acc end data */

}
