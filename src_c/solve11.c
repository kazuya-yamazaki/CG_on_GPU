#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pfem_util.h"
#include "allocate.h"
/** Paralle Version **/
extern FILE *fp_log;
extern void CG();
void SOLVE11()
{
  int  ERROR, ICFLAG=0;
  CHAR_LENGTH BUF;
/**
   +------------+
   | PARAMETERs |
   +------------+
**/
  ITER      = pfemIarray[0];
  RESID     = pfemRarray[0];

/**
   +------------------+
   | ITERATIVE solver |
   +------------------+
**/
  CG  ( N, NP, NPLU, D, AMAT, indexLU, itemLU, 
	B, X,  RESID, ITER, &ERROR, my_rank,
	NEIBPETOT,NEIBPE,IMPORT_INDEX,IMPORT_ITEM,
	EXPORT_INDEX,EXPORT_ITEM);

  ITERactual= ITER;
}

