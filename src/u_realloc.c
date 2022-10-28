
/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

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

/*
 Diehl program
*/

void *u_realloc(void* ptr,size_t size,const char *file,const int line, const char* ptr_name){

    /* reallocating a field with pointer ptr to size bytes */

  void *a;
  char *env;

  a=realloc(ptr,size);

  if(a==NULL && ptr!=NULL && size!=0){
    printf(" *ERROR in u_realloc: error allocating memory\n");
    printf("variable=%s, file=%s, line=%d, size(bytes)=%ld, oldaddress=%ld\n",ptr_name,file,line,size,(long int)ptr);
    exit(16);
  }
  else {
    if(log_realloc==-1) {
      log_realloc=0;
      env=getenv("CCX_LOG_ALLOC");
      if(env) {log_realloc=atoi(env);}
    }      
    if(log_realloc==1) {
      printf("REALLOCATION of variable %s, file %s, line=%d: size(bytes)=%ld, oldaddress= %ld,address= %ld\n",ptr_name,file,line,size,(long int)ptr,(long int)a);
    }      
    return(a);
  }
}
