
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"
/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */



#include "readfrd.h"


/* liest einen Record bis '\n'; uebergibt Anzahl gelesene Zeichen */
int frecord( FILE *handle1,  char *string)
{
  register int i, n, c;

  for (i=0; i<MAX_LINE_LENGTH-1; i++)
  {
    string[i] = getc(handle1);
    if (string[i] == '\n')
      {
      for (n=i+1; n<MAX_LINE_LENGTH; ++n) string[n] = '\0';
      return(i);
      }
    if (string[i] == '\r')
      {
      c = getc(handle1);
      if (c != '\n')
        ungetc(c, handle1);

      for (n=i+1; n<MAX_LINE_LENGTH; ++n) string[n] = '\0';
      return(i);
     }
    else if (string[i] == (char)EOF)
      {
      for (n=i+1; n<MAX_LINE_LENGTH; ++n) string[n] = '\0';
      return(i);
      }
  }
  string[MAX_LINE_LENGTH-1] = '\0';
  return(MAX_LINE_LENGTH-1);
}


