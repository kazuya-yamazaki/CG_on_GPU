/***
	pfem_util.h
***/
#include "mpi.h"
#include "precision.h"
#ifdef GLOBAL_VALUE_DEFINE
#define GLOBAL
#else
#define GLOBAL extern
#endif
/***
	+--------------+
	| MPI settings |
	+--------------+
***/
	GLOBAL int PEsmpTOT;
	GLOBAL int PETOT, my_rank, errno;
	GLOBAL char  penum[4], penum_left[4];
	GLOBAL char  fname[80], HEADER[80];
/***
	+-------------------------+
	| DISTRIBUTED MESH FILE's |
	+-------------------------+
***/
/***
	COMM. TABLE
***/
	GLOBAL int NEIBPETOT;
	GLOBAL int *IMPORT_INDEX, *IMPORT_ITEM;
	GLOBAL int *EXPORT_INDEX, *EXPORT_ITEM;
	GLOBAL int *NEIBPE;
/***
	CONNECTIVITIES & BOUNDARY nodes
***/
        GLOBAL int ICELTOT, ICELTOT_INT, NODGRPtot;
	GLOBAL KREAL **XYZ;
	GLOBAL KINT  **ICELNOD;
	GLOBAL KINT  **NODE_ID;
	GLOBAL KINT  **ELEM_ID;
	GLOBAL KINT  *intELEM_list;
        GLOBAL KINT  *NODGRP_INDEX, *NODGRP_ITEM;
        GLOBAL CHAR80 *NODGRP_NAME;
/***
	+-----------------+
	| MATRIX & SOLVER |
	+-----------------+
***/
/***
	MATRIX SCALARs 
***/
GLOBAL KINT N, NP, N2, NLU, NPLU, ELMCOLORtot;
/***
	MATRIX arrays
***/
        GLOBAL KREAL *D, *B, *X;
        GLOBAL KREAL *AMAT;
        GLOBAL KINT *indexLU, *itemLU;
GLOBAL KINT *INLU, *ELMCOLORindex, *ELMCOLORitem;
        GLOBAL KINT **IALU;
	GLOBAL KINT **IWKX;
/***
	PARAMETER's for LINEAR SOLVER
***/
	GLOBAL KINT ITER, ITERactual;
	GLOBAL KREAL RESID, SIGMA_DIAG, SIGMA;
/***
	+-------------+
	| PARAMETER's |
	+-------------+
***/
/***
	GENERAL PARAMETER's
***/
	GLOBAL KINT  pfemIarray[100];
	GLOBAL KREAL pfemRarray[100];
#ifdef GLOBAL_VALUE_DEFINE
	GLOBAL KREAL O8th= 0.125e0;
#else
	GLOBAL KREAL O8th;
#endif
/***
	PARAMETER's for FEM
***/
	GLOBAL KREAL PNQ[2][2][8], PNE[2][2][8], PNT[2][2][8];
	GLOBAL KREAL WEI[2], POS[2];
	GLOBAL KINT  NCOL1[100], NCOL2[100];
	GLOBAL KREAL SHAPE[2][2][2][8];
	GLOBAL KREAL PNX[2][2][2][8],PNY[2][2][2][8],PNZ[2][2][2][8];
	GLOBAL KREAL DETJ[2][2][2];
/***
	PROBLEM PARAMETER's
***/
	GLOBAL KREAL COND;
	GLOBAL KREAL QVOL;

