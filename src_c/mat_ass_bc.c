/**
 ** MAT_ASS_BC
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pfem_util.h"
#include "allocate.h"
extern FILE *fp_log;
void MAT_ASS_BC()
{
  int i,j,k,in,ib,ib0,icel;
  int in1,in2,in3,in4,in5,in6,in7,in8;
  int iq1,iq2,iq3,iq4,iq5,iq6,iq7,iq8;
  int iS,iE;
  double STRESS,VAL;
  double t0;
  
  IWKX=(KINT**) allocate_matrix(sizeof(KINT),NP,2);
  
  t0 = MPI_Wtime();
#pragma acc data create(IWKX[:NP][:2]) copyin(NODGRP_ITEM[:NODGRP_INDEX[NODGRPtot]]) \
                 present(B[:NP], D[:NP], AMAT[:NPLU], itemLU[:NPLU], indexLU[:NP+1])
{
  printf("acc data time=%f\n", MPI_Wtime()-t0);

/**
   Z=Zmax
**/
#pragma omp parallel for private (i,j)
#pragma acc parallel loop
  for(i=0;i<NP;i++){
    for(j=0;j<2;j++){
      IWKX[i][j]=0;
    }
  }
  
  ib0=-1;
  for( ib0=0;ib0<NODGRPtot;ib0++){
    if( strcmp(NODGRP_NAME[ib0].name,"Zmax") == 0 ) break;
  }
#pragma omp parallel for private (ib,in)    
#pragma acc parallel loop
  for( ib=NODGRP_INDEX[ib0];ib<NODGRP_INDEX[ib0+1];ib++){
    in=NODGRP_ITEM[ib];
    IWKX[in-1][0]=1;
  }

#pragma omp parallel for private (in,k)      
#pragma acc kernels
  for(in=0;in<NP;in++){
    if( IWKX[in][0] == 1 ){
      B[in]= 0.e0;
      D[in]= 1.e0;
      for(k=indexLU[in];k<indexLU[in+1];k++){
	AMAT[k]= 0.e0;
      }
    }
  }
#pragma omp parallel for private (in,k)       
#pragma acc kernels   
  for(in=0;in<NP;in++){
    for(k=indexLU[in];k<indexLU[in+1];k++){
      if (IWKX[itemLU[k]][0] == 1 ) {
	AMAT[k]= 0.e0;
      }
    }
  }

} /* acc end data */

}
