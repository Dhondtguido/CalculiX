
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

void *u_free(void* ptr,const char *file,const int line, const char* ptr_name){

    /* freeing a field with pointer ptr  */

  char *env;

  free(ptr);

  if(log_realloc==-1) {
      log_realloc=0;
      env=getenv("CCX_LOG_ALLOC");
      if(env) {log_realloc=atoi(env);}
  }      
  if(log_realloc==1) {
      printf("FREEING of variable %s, file %s, line=%d: oldaddress= %ld\n",ptr_name,file,line,(long int)ptr);
  }      
  return NULL;
}
