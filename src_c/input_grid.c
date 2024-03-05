/***
 *** INPUT_GRID
 ***/
/** PARALLEL VERSION **/
#include <stdio.h>
#include <stdlib.h>
#include "pfem_util.h"
#include "allocate.h"
/*** external functions **/
extern void ERROR_EXIT (int, int);
extern void DEFINE_FILE_NAME(char*,char*,int);
/** **/
void INPUT_GRID()
{
  FILE *fp;
  int i,j,k,ii,kk,kkk,nn,icel,iS,iE,ic0;
  int NTYPE,IMAT;
  
  DEFINE_FILE_NAME(HEADER, fname, my_rank);
  if( (fp=fopen(fname,"r")) == NULL){
    fprintf(stdout,"input file cannot be opened!\n");
    exit(1);
  }
/**
   NEIB-PE
**/
  fscanf(fp,"%d",&kkk);
  fscanf(fp,"%d",&NEIBPETOT);
  
  NEIBPE=(int*)allocate_vector(sizeof(int),NEIBPETOT);
  
  for(i=0;i<NEIBPETOT;i++) fscanf(fp,"%d",&NEIBPE[i]);
  
  for(i=0;i<NEIBPETOT;i++){
    if( NEIBPE[i] > PETOT-1 ){
      ERROR_EXIT (202,my_rank);
    }
  }
/**
   NODE
**/
  fscanf(fp,"%d %d",&NP,&N);
  
  XYZ    =(KREAL**)allocate_matrix(sizeof(KREAL),NP,3);
  NODE_ID=(KINT **)allocate_matrix(sizeof(KINT ),NP,2);
  
  for(i=0;i<NP;i++){
    for(j=0;j<3;j++){
      XYZ[i][j]=0.0;
    }
  }

  for(i=0;i<NP;i++){
    fscanf(fp,"%d %d %lf %lf %lf",&NODE_ID[i][0],&NODE_ID[i][1],&XYZ[i][0],&XYZ[i][1],&XYZ[i][2]);
  }
/**
   ELEMENT
**/
  fscanf(fp,"%d %d",&ICELTOT,&ICELTOT_INT);
  
  ICELNOD=(KINT**)allocate_matrix(sizeof(KINT),ICELTOT,8);
  intELEM_list=(KINT*)allocate_vector(sizeof(KINT),ICELTOT);
  ELEM_ID=(KINT**)allocate_matrix(sizeof(KINT),ICELTOT,2);

  for(i=0;i<ICELTOT;i++) fscanf(fp,"%d",&NTYPE);
  
  for(icel=0;icel<ICELTOT;icel++){
    fscanf(fp,"%d %d %d %d %d %d %d %d %d %d %d",
	   &ELEM_ID[icel][0],&ELEM_ID[icel][1],
	   &IMAT,
	   &ICELNOD[icel][0],&ICELNOD[icel][1],&ICELNOD[icel][2],&ICELNOD[icel][3],
	   &ICELNOD[icel][4],&ICELNOD[icel][5],&ICELNOD[icel][6],&ICELNOD[icel][7]);
  }

  for(ic0=0;ic0<ICELTOT_INT;ic0++) fscanf(fp,"%d",&intELEM_list[ic0]);

/**
   COMMUNICATION table
**/
  IMPORT_INDEX=(int*)allocate_vector(sizeof(int),NEIBPETOT+1);
  EXPORT_INDEX=(int*)allocate_vector(sizeof(int),NEIBPETOT+1);
  
  for(i=0;i<NEIBPETOT+1;++i) IMPORT_INDEX[i]=0;
  for(i=0;i<NEIBPETOT+1;++i) EXPORT_INDEX[i]=0;
  
  if( PETOT != 1 ) {
    for(i=1;i<=NEIBPETOT;i++) fscanf(fp,"%d",&IMPORT_INDEX[i]);
    nn=IMPORT_INDEX[NEIBPETOT];
    if( nn > 0 ) IMPORT_ITEM=(int*)allocate_vector(sizeof(int),nn);
    for(i=0;i<nn;i++) fscanf(fp,"%d ",&IMPORT_ITEM[i]);
    
    for(i=1;i<=NEIBPETOT;i++) fscanf(fp,"%d",&EXPORT_INDEX[i]);
    nn=EXPORT_INDEX[NEIBPETOT];
    if( nn > 0 ) EXPORT_ITEM=(int*)allocate_vector(sizeof(int),nn);
    for(i=0;i<nn;i++) fscanf(fp,"%d",&EXPORT_ITEM[i]);
  }
/**
   NODE grp. info.
**/
  fscanf(fp,"%d",&NODGRPtot);
  
  NODGRP_INDEX=(KINT*  )allocate_vector(sizeof(KINT),NODGRPtot+1);
  NODGRP_NAME =(CHAR80*)allocate_vector(sizeof(CHAR80),NODGRPtot);
  for(i=0;i<NODGRPtot+1;i++) NODGRP_INDEX[i]=0;
  
  for(i=0;i<NODGRPtot;i++) fscanf(fp,"%d",&NODGRP_INDEX[i+1]);
  nn=NODGRP_INDEX[NODGRPtot];
  NODGRP_ITEM=(KINT*)allocate_vector(sizeof(KINT),nn);
  
  for(k=0;k<NODGRPtot;k++){
    iS= NODGRP_INDEX[k];
    iE= NODGRP_INDEX[k+1];
    fscanf(fp,"%s",NODGRP_NAME[k].name);
    nn= iE - iS;
    if( nn != 0 ){
      for(kk=iS;kk<iE;kk++) fscanf(fp,"%d",&NODGRP_ITEM[kk]);
    }
  }
  fclose(fp);
}
