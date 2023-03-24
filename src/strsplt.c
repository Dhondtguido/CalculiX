
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

#include "readfrd.h"

/* split rec_str at each breakchar into dat[] */
/* dat should be unused before */
int strsplt( char *rec_str, char breakchar, char ***ptr)
{
  int i,j;
  int nextarg=0, letter=0, skip_breakchar=0;
  char **dat;

  if( (dat= (char **)malloc( 1*sizeof(char *))) == NULL )
    printf(" ERROR: malloc failed\n");
  if( (dat[0]= (char *)malloc(MAX_LINE_LENGTH *sizeof(char))) == NULL )
    printf(" ERROR: malloc failed\n");

  /* scan all args divided by breakchar */
  nextarg=0;letter=0;
  for(j=0; j<MAX_LINE_LENGTH; j++) dat[nextarg][j]='\0'; 
  for(i=0; i<MAX_LINE_LENGTH; i++)
  {
    if(rec_str[i]==(char)EOF) {  break; } 
    if(rec_str[i]=='\n') {  break; } 
    if(rec_str[i]==0) {  break; } 
    if((rec_str[i]==breakchar)&&(!skip_breakchar))
    {
      /* check if the former dat has chars */
      if(strlen(dat[nextarg]))
      {
        nextarg++;
        letter=0;
        if( (dat= (char **)realloc((char **)dat, (nextarg+1)*sizeof(char *))) == NULL )
          printf(" ERROR: realloc failed\n");
        if( (dat[nextarg]= (char *)malloc(MAX_LINE_LENGTH *sizeof(char))) == NULL )
          printf(" ERROR: malloc failed\n");
        for(j=0; j<MAX_LINE_LENGTH; j++) dat[nextarg][j]='\0'; 
      }
    }
    else
    {
      if(rec_str[i]=='"') skip_breakchar=!skip_breakchar;
      else if(skip_breakchar)
      {
        dat[nextarg][letter]=rec_str[i];
        letter++;
      }
      else if(rec_str[i]!=' ')
      {
        dat[nextarg][letter]=rec_str[i];
        letter++;
      }
    }
  }
  *ptr=dat;

  /* check if the former dat has chars */
  if(strlen(dat[nextarg])) nextarg++;

  /* for(i=0; i<nextarg; i++) printf("dat:%s|\n",dat[i]); */
  return(nextarg);
}

