
/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2014 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <stdlib.h>
#include "CalculiX.h"
extern int log_realloc;

//int log_realloc=-1;

/*
 Diehl program
*/

void *u_malloc(size_t size,const char *file,const int line, const char* ptr_name){

    /* allocating size bytes  */

  void *a;
  char *env;

  if(size==0){
    a=NULL;
    return(a);
  }
      
  a=malloc(size);
  if(a==NULL){
    printf("*ERROR in u_malloc: error allocating memory\n");
    printf("variable=%s, file=%s, line=%d, size=%ld\n",ptr_name,file,line,size);
    if(size<0){
	printf("\n It looks like you may need the i8 (integer*8) version of CalculiX\n");
    }
    exit(16);
  }
  else {
    if(log_realloc==-1) {
      log_realloc=0;
      env=getenv("CCX_LOG_ALLOC");
      if(env) {log_realloc=atoi(env);}
    }      
    if(log_realloc==1) {
	printf("ALLOCATION of variable %s, file %s, line=%d, size=%ld, address= %ld\n",ptr_name,file,line,size,(long int)a);
    }      
    return(a);
  }
}
