#include <stdio.h>
#include <string.h>
void DEFINE_FILE_NAME (char HEADERo[], char filename[],int my_rank)
{
  char string[80];
  sprintf(string,".%-d",my_rank);
  strcpy(filename,HEADERo);
  strcat(filename,string);
}
