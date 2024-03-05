/***
 *** OUTPUT_UCD
 ***/
/** PARALLEL VERSION **/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "pfem_util.h"
#include "allocate.h"
extern FILE *fp_log;
/** external function **/
extern void JACOBI();
void OUTPUT_UCD()
{
  KREAL *XYZ_G;
  KREAL *VAL_G,*WSarray;
  int   *NODEflag;
  int   *ICELNOD_L,*ICELNOD_G;
  int   *PEnode,*PEelem;
  int   *rcountsN,*displsN;
  int   *rcountsE,*displsE;
  int   NNtotG,NEtotG;
  int   NNlocal,NElocal;
  
  char ETYPE[7];
  
  FILE *fp;
  int i,ie,ip,in,icou;
  int icel,icel0;
  int iS0;
  int N0,N1,N3,N4;
  int in1,in2,in3,in4,in5,in6,in7,in8;
  double ZERO;
  double XX,YY,ZZ;
  double igmax;
  double  VAL_G_max;
  double  VAL_G_node;
/***
    +--------------------+
    | INTERNAL ELEMENT's |
    +--------------------+
***/
  NODEflag =(int *)allocate_vector(sizeof(int),NP);
  ICELNOD_L=(int *)allocate_vector(sizeof(int),8*ICELTOT_INT);
  
  for(i=0;i<NP;i++) NODEflag[i]=0;
  for(i=0;i<8*ICELTOT_INT;i++) ICELNOD_L[i]=0;
  
  for( icel0=0;icel0<ICELTOT_INT;icel0++){
    icel=intELEM_list[icel0];
    
    in1= ICELNOD[icel-1][0];
    in2= ICELNOD[icel-1][1];
    in3= ICELNOD[icel-1][2];
    in4= ICELNOD[icel-1][3];
    in5= ICELNOD[icel-1][4];
    in6= ICELNOD[icel-1][5];
    in7= ICELNOD[icel-1][6];
    in8= ICELNOD[icel-1][7];
    
    NODEflag[in1-1]= 1;
    NODEflag[in2-1]= 1;
    NODEflag[in3-1]= 1;
    NODEflag[in4-1]= 1;
    NODEflag[in5-1]= 1;
    NODEflag[in6-1]= 1;
    NODEflag[in7-1]= 1;
    NODEflag[in8-1]= 1;
  }
  
  NElocal= ICELTOT_INT;
  
  icou= 0;
  for(i=0;i<NP;i++){
    if( NODEflag[i] != 0 ){
      icou++;
      NODEflag[i]=icou;
    }
  }
/*** 
     +---------------+
     | GLOBAL ARRAYs |
     +---------------+
***/
  rcountsN=(int *)allocate_vector(sizeof(int),PETOT);
  displsN =(int *)allocate_vector(sizeof(int),PETOT+1);
  rcountsE=(int *)allocate_vector(sizeof(int),PETOT);
  displsE =(int *)allocate_vector(sizeof(int),PETOT+1);
  
  for(i=0;i<PETOT;i++) rcountsN[i]=0;
  for(i=0;i<PETOT+1;i++) displsN[i]=0;
  for(i=0;i<PETOT;i++) rcountsE[i]=0;
  for(i=0;i<PETOT+1;i++) displsE[i]=0;
  
  NNlocal= icou;
  
  
  MPI_Allgather( &NNlocal, 1, MPI_INT,  rcountsN, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather( &NElocal, 1, MPI_INT,  rcountsE, 1, MPI_INT, MPI_COMM_WORLD);
  
  for(ip=1;ip<=PETOT;ip++){
    displsN[ip]=displsN[ip-1]+rcountsN[ip-1];
    displsE[ip]=displsE[ip-1]+rcountsE[ip-1];
  }
  
  iS0=displsN[my_rank];
  
  for( icel0=0;icel0< ICELTOT_INT;icel0++){
    icel=intELEM_list[icel0];
    
    in1= NODEflag[ICELNOD[icel-1][0]-1];
    in2= NODEflag[ICELNOD[icel-1][1]-1];
    in3= NODEflag[ICELNOD[icel-1][2]-1];
    in4= NODEflag[ICELNOD[icel-1][3]-1];
    in5= NODEflag[ICELNOD[icel-1][4]-1];
    in6= NODEflag[ICELNOD[icel-1][5]-1];
    in7= NODEflag[ICELNOD[icel-1][6]-1];
    in8= NODEflag[ICELNOD[icel-1][7]-1];
    
    ICELNOD_L[8*icel0  ]= in1 + iS0;
    ICELNOD_L[8*icel0+1]= in2 + iS0;
    ICELNOD_L[8*icel0+2]= in3 + iS0;
    ICELNOD_L[8*icel0+3]= in4 + iS0;
    ICELNOD_L[8*icel0+4]= in5 + iS0;
    ICELNOD_L[8*icel0+5]= in6 + iS0;
    ICELNOD_L[8*icel0+6]= in7 + iS0;
    ICELNOD_L[8*icel0+7]= in8 + iS0;
  }
  
  NNtotG= displsN[PETOT];
  NEtotG= displsE[PETOT];
  
/***
    +--------------+
    | ALL_GATHER_V |
    +--------------+
***/

  ICELNOD_G=(int *)allocate_vector(sizeof(int),8*NEtotG);
  VAL_G    =(KREAL *)allocate_vector(sizeof(KREAL),  NNtotG);
  XYZ_G    =(KREAL *)allocate_vector(sizeof(KREAL),3*NNtotG);
  WSarray  =(KREAL *)allocate_vector(sizeof(KREAL),3*NNtotG);
  
/***
    ELEMENT CONNECTIVITY
***/
  for(ip=1;ip<=PETOT;ip++){
    displsE [ip]  =8* displsE[ip];
    rcountsE[ip-1]=8* rcountsE[ip-1];
  }
  
  MPI_Allgatherv( ICELNOD_L,  NElocal*8, MPI_INT,
		  ICELNOD_G,  rcountsE , &displsE[0],MPI_INT,
		  MPI_COMM_WORLD);
/***
    NODE COORDINATE/VAL
***/
  for(ip=1;ip<=PETOT;ip++){
    displsN [ip]  = 3* displsN[ip];
    rcountsN[ip-1]= 3*rcountsN[ip-1];
  }
  
  for(i=0;i<NP;i++){
    if( NODEflag[i] != 0 ){
      in=NODEflag[i];
      WSarray[3*in-3]= XYZ[i][0];
      WSarray[3*in-2]= XYZ[i][1];
      WSarray[3*in-1]= XYZ[i][2];
    }
  }
  
  MPI_Allgatherv( WSarray,  3*NNlocal, MPI_DOUBLE,
		  XYZ_G  ,  rcountsN , &displsN[0], MPI_DOUBLE,
		  MPI_COMM_WORLD);

  for(ip=1;ip<=PETOT;ip++){
    displsN [ip]=   displsN[ip]/3;
    rcountsN[ip-1]= rcountsN[ip-1]/3;
  }
  
  for(i=0;i<NP;i++){
    if( NODEflag[i] != 0 ) {
      in= NODEflag[i];
      WSarray[in-1]= X[i];
    }
  }
  
  MPI_Allgatherv( WSarray,  NNlocal, MPI_DOUBLE,
		  VAL_G  ,  rcountsN , &displsN[0], MPI_DOUBLE,
		  MPI_COMM_WORLD);
  
/***
    +----------+
    | AVS file |
    +----------+
***/
  if (my_rank == 0 ){
    if( (fp=fopen("test.inp","w")) == NULL){
      fprintf(stdout,"output file cannot be opened!\n");
      exit(1);
    }
    
    N0= 0;
    N1= 1;
    N3= 3;
    N4= 4;
    ZERO= 0.e0;
    
    fprintf(fp,"%8d%8d%8d%8d%8d\n",NNtotG,NEtotG,N1,N0,N0);
    
    for(i=0;i<NNtotG;i++){
      XX=XYZ_G[3*i];
      YY=XYZ_G[3*i+1];
      ZZ=XYZ_G[3*i+2];
      fprintf(fp,"%8d%16.6e%16.6e%16.6e\n",i+1,XX,YY,ZZ);
    }
    
    for(ie=0;ie<NEtotG;ie++){
      strcpy(ETYPE," hex  ");
      in1= ICELNOD_G[8*ie  ];
      in2= ICELNOD_G[8*ie+1];
      in3= ICELNOD_G[8*ie+2];
      in4= ICELNOD_G[8*ie+3];
      in5= ICELNOD_G[8*ie+4];
      in6= ICELNOD_G[8*ie+5];
      in7= ICELNOD_G[8*ie+6];
      in8= ICELNOD_G[8*ie+7];
      fprintf(fp,"%8d%3d%-8s%8d%8d%8d%8d%8d%8d%8d%8d\n",
	      ie+1,N1,ETYPE,in1, in2, in3, in4, in5, in6, in7, in8);
    }
    
    fprintf(fp,"%3d%3d\n",N1, N1);
    fprintf(fp,"temp,temp\n");
    
    igmax    = 0;
    VAL_G_max= 0.e0;
    
    for(i=0;i<NNtotG;i++){
      fprintf(fp,"%8d%16.6e\n",i+1,VAL_G[i]);
    }
    
    fclose(fp);
    
  }
}

