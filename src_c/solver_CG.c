/**
 ** CG
 **/
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "precision.h"
#include "allocate.h"
extern FILE *fp_log;
extern void SOLVER_SEND_RECV (); 
/***
    CG solves the linear system Ax = b using the Conjugate Gradient 
    iterative method with the following preconditioners
***/
void  CG  (
	   KINT N,KINT NP,KINT NPLU, KREAL D[],
	   KREAL AMAT[],KINT indexLU[], KINT itemLU[],
	   KREAL B[],KREAL X[],KREAL RESID,KINT ITER, KINT *ERROR,int my_rank,
	   int NEIBPETOT,int NEIBPE[],
	   int IMPORT_INDEX[], int IMPORT_ITEM[],
	   int EXPORT_INDEX[], int EXPORT_ITEM[])
{
  int i,j,k;
  int ieL,isL,ieU,isU;
  double WVAL;
  double BNRM20,BNRM2,DNRM20,DNRM2;
  double S1_TIME,E1_TIME;
  double ALPHA,BETA;
  double C1,C10,RHO,RHO0,RHO1;
  int    iterPRE;
  
  KREAL *WS,*WR;
  KREAL **WW;
  
  KINT R=0,Z=1,Q=1,P=2,DD=3;
  KINT MAXIT;
  KREAL TOL;
  
  double COMPtime,COMMtime,R1;
  double START_TIME,END_TIME;
  
/**
   +-------+
   | INIT. |
   +-------+
**/
  ERROR= 0;
  
  COMPtime=0.0;
  COMMtime=0.0;
  
  WW=(KREAL**) allocate_matrix(sizeof(KREAL),4,NP);
  WS=(KREAL* ) allocate_vector(sizeof(KREAL),  NP);
  WR=(KREAL* ) allocate_vector(sizeof(KREAL),  NP);
  
#pragma acc data present(itemLU[:NPLU], indexLU[:NP+1]) \
	present(B[:NP], D[:NP], AMAT[:NPLU], X[:NP]) \
    create(WW[:4][:NP], WS[:NP], WR[:NP]) \
	copyin(IMPORT_INDEX[:NEIBPETOT+1], IMPORT_ITEM[:IMPORT_INDEX[NEIBPETOT]]) \
	copyin(EXPORT_INDEX[:NEIBPETOT+1], EXPORT_ITEM[:EXPORT_INDEX[NEIBPETOT]])
{
  
  MAXIT  = ITER;
  TOL   = RESID;          
  
#pragma omp parallel for private (i)
#pragma acc kernels
  for(i=0;i<NP;i++) {
    WW[0][i]=0.0;
    WW[1][i]=0.0;
    WW[2][i]=0.0;
    WW[3][i]=0.0;
    WS[i]= 0.0;
    WR[i]= 0.0;
    X [i]= 0.0;
  }     

/**
   +-----------------------+
   | {r0}= {b} - [A]{xini} |
   +-----------------------+
**/
/**
 ** INTERFACE data EXCHANGE
**/
  SOLVER_SEND_RECV
    ( N, NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
      EXPORT_INDEX, EXPORT_ITEM, WS, WR, X , my_rank);
  
/**
   BEGIN calculation
**/
#pragma omp parallel for private (j,k,i,WVAL)
#pragma acc parallel loop private(i,WVAL)
  for(j=0;j<N;j++){
    WW[DD][j]= 1.0/D[j];
    WVAL= B[j] - D[j]*X[j];
    
    for( k=indexLU[j];k<indexLU[j+1];k++){
      i= itemLU[k];
      WVAL+=  -AMAT[k]*X[i];
    }
    WW[R][j]= WVAL;
  }
  
  BNRM20= 0.e0;
#pragma omp parallel for private (i) reduction (+:BNRM20)
#pragma acc kernels loop             reduction (+:BNRM20)
  for(i=0;i<N;i++){
    BNRM20+= B[i]*B[i];
  }
  
  MPI_Allreduce (&BNRM20, &BNRM2, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
  
  if (BNRM2 == 0.e0) BNRM2= 1.e0;
  
  ITER = 0;
  
  S1_TIME= MPI_Wtime();
  for( ITER=1;ITER<= MAXIT;ITER++){
/**
************************************************* Conjugate Gradient Iteration
**/

/**
   +----------------+
   | {z}= [Minv]{r} |
   +----------------+
**/  
#pragma omp parallel for private (i) 
#pragma acc kernels
    for(i=0;i<N;i++){
      WW[Z][i]= WW[DD][i]*WW[R][i];
    }
/**
   +---------------+
   | {RHO}= {r}{z} |
   +---------------+
**/
    RHO0= 0.e0;
	
#pragma omp parallel for private (i) reduction (+:RHO0)
#pragma acc kernels loop             reduction (+:RHO0)
    for(i=0;i<N;i++){
      RHO0+= WW[R][i]*WW[Z][i];
    }

    MPI_Allreduce (&RHO0, &RHO, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
/**
   +-----------------------------+
   | {p} = {z} if      ITER=1    |
   | BETA= RHO / RHO1  otherwise |
   +-----------------------------+
**/
    if( ITER == 1 ){
#pragma omp parallel for private (i)
#pragma acc kernels
      for(i=0;i<N;i++){
	WW[P][i]=WW[Z][i];
      }
    }else{
      BETA= RHO / RHO1;
#pragma omp parallel for private (i)
#pragma acc kernels
      for(i=0;i<N;i++){
	WW[P][i]=WW[Z][i] + BETA*WW[P][i];
      }
    }
/**
   +-------------+
   | {q}= [A]{p} |
   +-------------+
**/      

    SOLVER_SEND_RECV
      ( N, NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
	EXPORT_INDEX, EXPORT_ITEM, WS, WR, WW[P], my_rank);

#pragma omp parallel for private (j,i,k,WVAL)
#pragma acc parallel loop
    for( j=0;j<N;j++){
      WVAL= D[j] * WW[P][j];
#pragma acc loop seq
      for(k=indexLU[j];k<indexLU[j+1];k++){
	i=itemLU[k];
	WVAL+= AMAT[k] * WW[P][i];
      }
      WW[Q][j]=WVAL;
    }

/**
   +---------------------+
   | ALPHA= RHO / {p}{q} |
   +---------------------+
**/
    C10= 0.e0;
#pragma omp parallel for private (i) reduction (+:C10)
#pragma acc kernels loop             reduction (+:C10)
    for(i=0;i<N;i++){
      C10+=WW[P][i]*WW[Q][i];
    }
    MPI_Allreduce (&C10,&C1, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    
    ALPHA= RHO / C1;
    
/**
   +----------------------+
   | {x}= {x} + ALPHA*{p} |
   | {r}= {r} - ALPHA*{q} |
   +----------------------+
**/
#pragma omp parallel for private (i)
#pragma acc parallel loop
    for(i=0;i<N;i++){
      X [i]   +=  ALPHA *WW[P][i];
      WW[R][i]+= -ALPHA *WW[Q][i];
    }

    DNRM20= 0.e0;
#pragma omp parallel for private (i) reduction (+:DNRM20)
#pragma acc kernels loop             reduction (+:DNRM20)
    for(i=0;i<N;i++){
      DNRM20+=WW[R][i]*WW[R][i];
    }
    MPI_Allreduce (&DNRM20,&DNRM2, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    RESID= sqrt(DNRM2/BNRM2);

/** ##### ITERATION HISTORY ***/
    if( my_rank == 0 ) fprintf(stdout,"%d %e\n",ITER,RESID);
    if( my_rank == 0 ) fprintf(fp_log,"%d %e\n",ITER,RESID);

    if ( RESID <= TOL   ) break;
    if ( ITER  == MAXIT ) *ERROR= -300;
    
    RHO1 = RHO ;                                                           
}
/** **/
/***
    INTERFACE data EXCHANGE
***/
  E1_TIME= MPI_Wtime();
  COMPtime= E1_TIME - S1_TIME;
  
  R1= 100.e0 * ( 1.e0 - COMMtime/COMPtime );

  if (my_rank == 0) {
    fprintf(stdout,"### comm.      :%e\n",COMMtime);
    fprintf(stdout,"### work ratio :%e\n",R1);
  }
  
  SOLVER_SEND_RECV
    ( N, NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
      EXPORT_INDEX, EXPORT_ITEM, WS, WR, X, my_rank);

} /* acc end data */

  free ( (KREAL**)WW);
  deallocate_vector ( (KREAL**)WR);
  deallocate_vector( (KREAL**)WS);
}
