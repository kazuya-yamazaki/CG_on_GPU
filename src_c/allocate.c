/******************************************************************************
*         Allocate function                                                   *
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
/**************************
  allocate vector
**************************/
void* allocate_vector(long size,int m)
{
  void *a;
  if ( ( a=(void * )malloc( m * size ) ) == NULL ) {
     fprintf(stdout,"Error:Memory does not enough! in vector \n");
     exit(1);
  }
  return a;
}
/**************************
  deallocate vector
**************************/
void deallocate_vector(void *a)
{
  free( a );
}
/**************************
  allocate matrix 
**************************/
void** allocate_matrix(int size,int m,int n)
{
  void **aa;
  int i;
  if ( ( aa=(void ** )malloc( m * sizeof(void* ) ) ) == NULL ) {
     fprintf(stdout,"Error:Memory does not enough! aa in matrix \n");
     exit(1);
  }
  if ( ( aa[0]=(void * )malloc( m * n * size ) ) == NULL ) {
     fprintf(stdout,"Error:Memory does not enough! in matrix \n");
     exit(1);
  }
	for(i=1;i<m;i++) aa[i]=(char*)aa[i-1]+size*n;
  return aa;
}
/**************************
  deallocate matrix 
**************************/
void deallocate_matrix(void **aa)
{
	free( (void *) aa[0] );
	free( (void **) aa   );
}

