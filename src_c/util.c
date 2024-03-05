#include <stdio.h>
#include <math.h>
#include "precision.h"
#include "allocate.h"
/**
 ** JACOBI
 **/
#pragma acc routine seq
void JACOBI(
	    KREAL DETJ[2][2][2],
	    KREAL PNQ[2][2][8],KREAL PNE[2][2][8],KREAL PNT[2][2][8],
	    KREAL PNX[2][2][2][8],KREAL PNY[2][2][2][8],KREAL PNZ[2][2][2][8],
	    KREAL X1,KREAL X2,KREAL X3,KREAL X4,KREAL X5,KREAL X6,KREAL X7,KREAL X8,
	    KREAL Y1,KREAL Y2,KREAL Y3,KREAL Y4,KREAL Y5,KREAL Y6,KREAL Y7,KREAL Y8,
	    KREAL Z1,KREAL Z2,KREAL Z3,KREAL Z4,KREAL Z5,KREAL Z6,KREAL Z7,KREAL Z8)
{
/**
	calculates JACOBIAN & INVERSE JACOBIAN
	             dNi/dx, dNi/dy & dNi/dz         
**/ 
  int ip,jp,kp;
  double dXdQ,dYdQ,dZdQ,dXdE,dYdE,dZdE,dXdT,dYdT,dZdT;
  double coef;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  
  for(ip=0;ip<2;ip++){
    for(jp=0;jp<2;jp++){
      for(kp=0;kp<2;kp++){
	PNX[ip][jp][kp][0]=0.0;
	PNX[ip][jp][kp][1]=0.0;
	PNX[ip][jp][kp][2]=0.0;
	PNX[ip][jp][kp][3]=0.0;
	PNX[ip][jp][kp][4]=0.0;
	PNX[ip][jp][kp][5]=0.0;
	PNX[ip][jp][kp][6]=0.0;
	PNX[ip][jp][kp][7]=0.0;
	
	PNY[ip][jp][kp][0]=0.0;
	PNY[ip][jp][kp][1]=0.0;
	PNY[ip][jp][kp][2]=0.0;
	PNY[ip][jp][kp][3]=0.0;
	PNY[ip][jp][kp][4]=0.0;
	PNY[ip][jp][kp][5]=0.0;
	PNY[ip][jp][kp][6]=0.0;
	PNY[ip][jp][kp][7]=0.0;
	
	PNZ[ip][jp][kp][0]=0.0;
	PNZ[ip][jp][kp][1]=0.0;
	PNZ[ip][jp][kp][2]=0.0;
	PNZ[ip][jp][kp][3]=0.0;
	PNZ[ip][jp][kp][4]=0.0;
	PNZ[ip][jp][kp][5]=0.0;
	PNZ[ip][jp][kp][6]=0.0;
	PNZ[ip][jp][kp][7]=0.0;
		
/**    
       DETERMINANT of the JACOBIAN
**/
	dXdQ = PNQ[jp][kp][0]*X1 + PNQ[jp][kp][1]*X2                               
	  + PNQ[jp][kp][2]*X3 + PNQ[jp][kp][3]*X4
	  + PNQ[jp][kp][4]*X5 + PNQ[jp][kp][5]*X6  
	  + PNQ[jp][kp][6]*X7 + PNQ[jp][kp][7]*X8;
	dYdQ = PNQ[jp][kp][0]*Y1 + PNQ[jp][kp][1]*Y2                               
	  + PNQ[jp][kp][2]*Y3 + PNQ[jp][kp][3]*Y4
	  + PNQ[jp][kp][4]*Y5 + PNQ[jp][kp][5]*Y6  
	  + PNQ[jp][kp][6]*Y7 + PNQ[jp][kp][7]*Y8;
	dZdQ = PNQ[jp][kp][0]*Z1 + PNQ[jp][kp][1]*Z2                               
	  + PNQ[jp][kp][2]*Z3 + PNQ[jp][kp][3]*Z4
	  + PNQ[jp][kp][4]*Z5 + PNQ[jp][kp][5]*Z6  
	  + PNQ[jp][kp][6]*Z7 + PNQ[jp][kp][7]*Z8;
	dXdE = PNE[ip][kp][0]*X1 + PNE[ip][kp][1]*X2                               
	  + PNE[ip][kp][2]*X3 + PNE[ip][kp][3]*X4
	  + PNE[ip][kp][4]*X5 + PNE[ip][kp][5]*X6  
	  + PNE[ip][kp][6]*X7 + PNE[ip][kp][7]*X8;			
	dYdE = PNE[ip][kp][0]*Y1 + PNE[ip][kp][1]*Y2                               
	  + PNE[ip][kp][2]*Y3 + PNE[ip][kp][3]*Y4
	  + PNE[ip][kp][4]*Y5 + PNE[ip][kp][5]*Y6  
	  + PNE[ip][kp][6]*Y7 + PNE[ip][kp][7]*Y8;
	dZdE = PNE[ip][kp][0]*Z1 + PNE[ip][kp][1]*Z2                               
	  + PNE[ip][kp][2]*Z3 + PNE[ip][kp][3]*Z4
	  + PNE[ip][kp][4]*Z5 + PNE[ip][kp][5]*Z6  
	  + PNE[ip][kp][6]*Z7 + PNE[ip][kp][7]*Z8;
	dXdT = PNT[ip][jp][0]*X1 + PNT[ip][jp][1]*X2                               
	  + PNT[ip][jp][2]*X3 + PNT[ip][jp][3]*X4
	  + PNT[ip][jp][4]*X5 + PNT[ip][jp][5]*X6  
	  + PNT[ip][jp][6]*X7 + PNT[ip][jp][7]*X8;			
	dYdT = PNT[ip][jp][0]*Y1 + PNT[ip][jp][1]*Y2                               
	  + PNT[ip][jp][2]*Y3 + PNT[ip][jp][3]*Y4
	  + PNT[ip][jp][4]*Y5 + PNT[ip][jp][5]*Y6  
	  + PNT[ip][jp][6]*Y7 + PNT[ip][jp][7]*Y8;
	dZdT = PNT[ip][jp][0]*Z1 + PNT[ip][jp][1]*Z2                               
	  + PNT[ip][jp][2]*Z3 + PNT[ip][jp][3]*Z4
	  + PNT[ip][jp][4]*Z5 + PNT[ip][jp][5]*Z6  
	  + PNT[ip][jp][6]*Z7 + PNT[ip][jp][7]*Z8;
	DETJ[ip][jp][kp]= dXdQ*(dYdE*dZdT-dZdE*dYdT) +
	  dYdQ*(dZdE*dXdT-dXdE*dZdT) + 
	  dZdQ*(dXdE*dYdT-dYdE*dXdT);
/**
   INVERSE JACOBIAN
**/
	
	coef=1.0 / DETJ[ip][jp][kp];
	
	a11= coef * ( dYdE*dZdT - dZdE*dYdT );
	a12= coef * ( dZdQ*dYdT - dYdQ*dZdT );
	a13= coef * ( dYdQ*dZdE - dZdQ*dYdE );
	
	a21= coef * ( dZdE*dXdT - dXdE*dZdT );
	a22= coef * ( dXdQ*dZdT - dZdQ*dXdT );
	a23= coef * ( dZdQ*dXdE - dXdQ*dZdE );
	
	a31= coef * ( dXdE*dYdT - dYdE*dXdT );
	a32= coef * ( dYdQ*dXdT - dXdQ*dYdT );
	a33= coef * ( dXdQ*dYdE - dYdQ*dXdE );
	
	DETJ[ip][jp][kp]=fabs(DETJ[ip][jp][kp]);
	
/**
	set the dNi/dX, dNi/dY & dNi/dZ components
**/
	PNX[ip][jp][kp][0]=a11*PNQ[jp][kp][0]+a12*PNE[ip][kp][0]+a13*PNT[ip][jp][0];
	PNX[ip][jp][kp][1]=a11*PNQ[jp][kp][1]+a12*PNE[ip][kp][1]+a13*PNT[ip][jp][1];
	PNX[ip][jp][kp][2]=a11*PNQ[jp][kp][2]+a12*PNE[ip][kp][2]+a13*PNT[ip][jp][2];
	PNX[ip][jp][kp][3]=a11*PNQ[jp][kp][3]+a12*PNE[ip][kp][3]+a13*PNT[ip][jp][3];
	PNX[ip][jp][kp][4]=a11*PNQ[jp][kp][4]+a12*PNE[ip][kp][4]+a13*PNT[ip][jp][4];
	PNX[ip][jp][kp][5]=a11*PNQ[jp][kp][5]+a12*PNE[ip][kp][5]+a13*PNT[ip][jp][5];
	PNX[ip][jp][kp][6]=a11*PNQ[jp][kp][6]+a12*PNE[ip][kp][6]+a13*PNT[ip][jp][6];
	PNX[ip][jp][kp][7]=a11*PNQ[jp][kp][7]+a12*PNE[ip][kp][7]+a13*PNT[ip][jp][7];
	
	PNY[ip][jp][kp][0]=a21*PNQ[jp][kp][0]+a22*PNE[ip][kp][0]+a23*PNT[ip][jp][0];
	PNY[ip][jp][kp][1]=a21*PNQ[jp][kp][1]+a22*PNE[ip][kp][1]+a23*PNT[ip][jp][1];
	PNY[ip][jp][kp][2]=a21*PNQ[jp][kp][2]+a22*PNE[ip][kp][2]+a23*PNT[ip][jp][2];
	PNY[ip][jp][kp][3]=a21*PNQ[jp][kp][3]+a22*PNE[ip][kp][3]+a23*PNT[ip][jp][3];
	PNY[ip][jp][kp][4]=a21*PNQ[jp][kp][4]+a22*PNE[ip][kp][4]+a23*PNT[ip][jp][4];
	PNY[ip][jp][kp][5]=a21*PNQ[jp][kp][5]+a22*PNE[ip][kp][5]+a23*PNT[ip][jp][5];
	PNY[ip][jp][kp][6]=a21*PNQ[jp][kp][6]+a22*PNE[ip][kp][6]+a23*PNT[ip][jp][6];
	PNY[ip][jp][kp][7]=a21*PNQ[jp][kp][7]+a22*PNE[ip][kp][7]+a23*PNT[ip][jp][7];
	
	PNZ[ip][jp][kp][0]=a31*PNQ[jp][kp][0]+a32*PNE[ip][kp][0]+a33*PNT[ip][jp][0];
	PNZ[ip][jp][kp][1]=a31*PNQ[jp][kp][1]+a32*PNE[ip][kp][1]+a33*PNT[ip][jp][1];
	PNZ[ip][jp][kp][2]=a31*PNQ[jp][kp][2]+a32*PNE[ip][kp][2]+a33*PNT[ip][jp][2];
	PNZ[ip][jp][kp][3]=a31*PNQ[jp][kp][3]+a32*PNE[ip][kp][3]+a33*PNT[ip][jp][3];
	PNZ[ip][jp][kp][4]=a31*PNQ[jp][kp][4]+a32*PNE[ip][kp][4]+a33*PNT[ip][jp][4];
	PNZ[ip][jp][kp][5]=a31*PNQ[jp][kp][5]+a32*PNE[ip][kp][5]+a33*PNT[ip][jp][5];
	PNZ[ip][jp][kp][6]=a31*PNQ[jp][kp][6]+a32*PNE[ip][kp][6]+a33*PNT[ip][jp][6];
	PNZ[ip][jp][kp][7]=a31*PNQ[jp][kp][7]+a32*PNE[ip][kp][7]+a33*PNT[ip][jp][7];
      }
    }
  }
}
/**
 ** mSORT
 **/
void mSORT(KINT STEM[], KINT INUM[], int NN)
{
  int ii,jj;
  int ITEM;
  
  for(ii=1;ii<=NN;ii++){
    INUM[ii-1]=ii;
  }
  
  for( ii=1;ii<=NN-1;ii++){
    for( jj=1;jj<=NN-ii;jj++){
      if( STEM[INUM[jj-1]-1] <  STEM[INUM[jj]-1] ){
	ITEM=INUM[jj];
	INUM[jj  ]=INUM[jj-1];
	INUM[jj-1]=ITEM;
      }
    }
  }
}
/**
 ** matconSORT
 **/
void  matconSORT( KINT STEM[], KINT INUM[], int N, int NN)
{
  int i,k,ii,ik1,ik2,icon;
  int **ISTACK;
  int ICONmax;
  
  ISTACK=(int**)allocate_matrix(sizeof(int),NN+2,2);
  
  ISTACK[0][0]=0;
  ISTACK[0][1]=0;
  
  for(i=1;i<=N;i++){
    INUM[i-1]=i;
    STEM[i-1]++;
  }
  
  for(i=1;i<=N+1;i++){
    ISTACK[i-1][0]=0;
  }
  
  ICONmax= -N;
  
  for(i=0;i<N;i++){
    ii=STEM[i];
    if( ii > ICONmax ){
      ICONmax=ii;
    }
    ISTACK[ii-1][0]++;
  }
  
  for(k=1;k<=ICONmax;k++){
    ISTACK[k][0]+=ISTACK[k-1][0];
    ISTACK[k][1] =ISTACK[k  ][0];
  }
  
  ISTACK[0][1]=ISTACK[1][1];
  
  for(k=1;k<=ICONmax;k++){
    ik1=ICONmax - k;
    ik2=ik1     + 1;
    ISTACK[k][0]=ISTACK[ik2][1]-ISTACK[ik1][1]+ISTACK[k-1][0];
  }
  
  for(k=1;k<=ICONmax;k++){
    ISTACK[k][1]= 0;
  }
  
  for(i=0;i<N;i++){
    ii=STEM[i];
    icon=ISTACK[ii][1]+1;
    ISTACK[ii][1]=icon;
    INUM[ISTACK[ICONmax-ii+1-1][0]+icon-1]= i;
  }
  
  for(i=0;i<N;i++){
    STEM[i]+=-1;
  }
  
  deallocate_vector(ISTACK);
}

