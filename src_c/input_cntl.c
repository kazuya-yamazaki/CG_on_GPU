/***
 *** INPUT_CNTL
 ***/
/** PARALLEL VERSION **/
#include <stdio.h>
#include <stdlib.h>
#include "pfem_util.h"
/*** external functions **/
extern void ERROR_EXIT (int, int);
/** **/
void INPUT_CNTL()
{
  FILE *fp;
	
  if( my_rank == 0 ){
    if( (fp=fopen("INPUT.DAT","r")) == NULL){
      fprintf(stdout,"input file cannot be opened!\n");
      exit(1);
    }
    fscanf(fp,"%s",HEADER);
    fscanf(fp,"%d",&ITER);
    fscanf(fp, "%lf %lf", &COND, &QVOL);
    fscanf(fp, "%lf",     &RESID);
    fclose(fp);
  }
  
  MPI_Bcast(HEADER ,80,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&ITER   , 1,MPI_INTEGER,0,MPI_COMM_WORLD);
  MPI_Bcast(&COND   , 1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&QVOL   , 1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&RESID  , 1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  pfemRarray[0]= RESID;
  pfemIarray[0]= ITER;
}


