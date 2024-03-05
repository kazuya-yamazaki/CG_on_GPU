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
static void FIND_TS_NODE (int,int);
void MAT_CON0()
{
  int i,j,k,icel,in;
  int in1,in2,in3,in4,in5,in6,in7,in8;
  int NN;
  
  NLU= 26;
  
  INLU=(KINT* )allocate_vector(sizeof(KINT),NP);
  IALU=(KINT**)allocate_matrix(sizeof(KINT),NP,NLU);
  
  for(i=0;i<NP;i++) INLU[i]=0;
  for(i=0;i<NP;i++) for(j=0;j<NLU;j++) IALU[i][j]=0;
  
  for( icel=0;icel< ICELTOT;icel++){
    in1=ICELNOD[icel][0];
    in2=ICELNOD[icel][1];
    in3=ICELNOD[icel][2];
    in4=ICELNOD[icel][3];
    in5=ICELNOD[icel][4];
    in6=ICELNOD[icel][5];
    in7=ICELNOD[icel][6];
    in8=ICELNOD[icel][7];
    
    FIND_TS_NODE (in1,in2);
    FIND_TS_NODE (in1,in3);
    FIND_TS_NODE (in1,in4);
    FIND_TS_NODE (in1,in5);
    FIND_TS_NODE (in1,in6);
    FIND_TS_NODE (in1,in7);
    FIND_TS_NODE (in1,in8);
    
    FIND_TS_NODE (in2,in1);
    FIND_TS_NODE (in2,in3);
    FIND_TS_NODE (in2,in4);
    FIND_TS_NODE (in2,in5);
    FIND_TS_NODE (in2,in6);
    FIND_TS_NODE (in2,in7);
    FIND_TS_NODE (in2,in8);
    
    FIND_TS_NODE (in3,in1);
    FIND_TS_NODE (in3,in2);
    FIND_TS_NODE (in3,in4);
    FIND_TS_NODE (in3,in5);
    FIND_TS_NODE (in3,in6);
    FIND_TS_NODE (in3,in7);
    FIND_TS_NODE (in3,in8);
    
    FIND_TS_NODE (in4,in1);
    FIND_TS_NODE (in4,in2);
    FIND_TS_NODE (in4,in3);
    FIND_TS_NODE (in4,in5);
    FIND_TS_NODE (in4,in6);
    FIND_TS_NODE (in4,in7);
    FIND_TS_NODE (in4,in8);

    FIND_TS_NODE (in5,in1);
    FIND_TS_NODE (in5,in2);
    FIND_TS_NODE (in5,in3);
    FIND_TS_NODE (in5,in4);
    FIND_TS_NODE (in5,in6);
    FIND_TS_NODE (in5,in7);
    FIND_TS_NODE (in5,in8);
    
    FIND_TS_NODE (in6,in1);
    FIND_TS_NODE (in6,in2);
    FIND_TS_NODE (in6,in3);
    FIND_TS_NODE (in6,in4);
    FIND_TS_NODE (in6,in5);
    FIND_TS_NODE (in6,in7);
    FIND_TS_NODE (in6,in8);
    
    FIND_TS_NODE (in7,in1);
    FIND_TS_NODE (in7,in2);
    FIND_TS_NODE (in7,in3);
    FIND_TS_NODE (in7,in4);
    FIND_TS_NODE (in7,in5);
    FIND_TS_NODE (in7,in6);
    FIND_TS_NODE (in7,in8);

    FIND_TS_NODE (in8,in1);
    FIND_TS_NODE (in8,in2);
    FIND_TS_NODE (in8,in3);
    FIND_TS_NODE (in8,in4);
    FIND_TS_NODE (in8,in5);
    FIND_TS_NODE (in8,in6);
    FIND_TS_NODE (in8,in7);
  }
  
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
static void FIND_TS_NODE (int ip1,int ip2)
{
  int kk,icou;
  for(kk=1;kk<=INLU[ip1-1];kk++){
    if(ip2 == IALU[ip1-1][kk-1]) return;
  }
  icou=INLU[ip1-1]+1;
  IALU[ip1-1][icou-1]=ip2;
  INLU[ip1-1]=icou;
  return;
}

