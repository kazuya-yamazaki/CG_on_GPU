/**
 ** MAT_CON0
 **/
#include <stdio.h>
#include "pfem_util.h"
#include "allocate.h"
extern FILE *fp_log;
/*** external functions ***/
extern void mSORT(int*, int*, int);
/*** static functuons ***/
#pragma acc routine seq
static void FIND_TS_NODE (KINT*, KINT**, int,int,int,int);
void MAT_CON0()
{
  int i,j,k,icel,in;
  int i1, i2, in1, in2;
  int NN;
  
  NLU= 26;
  
  INLU=(KINT* )allocate_vector(sizeof(KINT),NP);
  IALU=(KINT**)allocate_matrix(sizeof(KINT),NP,NLU);
  
#pragma acc data copyout(INLU[:NP], IALU[:NP][:NLU]) copyin(ICELNOD[:ICELTOT][:8])
{
#pragma acc kernels loop independent
  for(i=0;i<NP;i++) INLU[i]=0;

#pragma acc kernels loop independent collapse(2)
  for(i=0;i<NP;i++) for(j=0;j<NLU;j++) IALU[i][j]=0;

#pragma acc kernels loop independent collapse(3) private(in1, in2)
  for( icel=0; icel<ICELTOT; icel++){
    for( i1=0; i1<8; i1++){
      for( i2=0; i2<8; i2++){
        if( i1 == i2 ) continue;
        in1 = ICELNOD[icel][i1];
        in2 = ICELNOD[icel][i2];
        FIND_TS_NODE (INLU,IALU,NP,NLU,in1,in2);
      }
    }
  }
}  /* acc end data */
  
  for(in=0;in<N;in++){
    NN=INLU[in];
    for (k=0;k<NN;k++){
      NCOL1[k]=IALU[in][k];
    }
    mSORT(NCOL1,NCOL2,NN);
    for(k=NN;k>0;k--){
      IALU[in][NN-k]= NCOL1[NCOL2[k-1]-1];
    }
  }
  
}
/***
 *** FIND_TS_NODE
**/
static void FIND_TS_NODE (KINT* INLU, KINT** IALU, int NP, int NLU, int ip1,int ip2)
{
  int inlu_ip1, kk, ks, icou, ind, retry;

  ks = 0;
  while(1){

#pragma acc atomic read
    inlu_ip1 = INLU[ip1-1];

    if (inlu_ip1 > NLU){
      /* If multiple threads increment INLU[ip1] at the same time, it may
         exceed the upper limit. In that case, start over. */
      continue;
    }
        
    retry = 0;
    for(kk=ks; kk<inlu_ip1; kk++){
#pragma acc atomic read
      ind = IALU[ip1-1][kk];
      
      if (ind == 0){
        /* It is possible that I might read the array value before other thread
           store the value in IALU[ip1][kk]. In that case, start over. */
        retry = 1;
        break;
      }else if (ip2 == ind){
        return;
      }
      ks = kk + 1;
    }
    if (retry) continue;

#pragma acc atomic capture
    icou = INLU[ip1-1]++;

    if (icou != inlu_ip1){
      /* Other threads may have incremented INLU(ip1).
         In that case, start over. */
#pragma acc atomic update
      INLU[ip1-1]--;
      
      continue;
    }
        
#pragma acc atomic write
    IALU[ip1-1][icou]= ip2;

    break;

  }
  
}

