/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

void sensitivity_out(char *jobnamec,double *dgdxglob,ITG *neq,ITG *nobject,
       double *g0){
            
  char sensitivities[132]="",nominal[132]="";
  
  ITG i=0,iobject=0;
      
  FILE *f1;
           		 
  /* writing the sensitivities in the sen-file for optimizer */
        	
  strcpy2(sensitivities,jobnamec,132);
  strcat(sensitivities,".sen");
  
  if((f1=fopen(sensitivities,"w"))==NULL){
      printf(" *ERROR in sensitivity_out: cannot open sensitivity vector file for writing...");
      
      exit(0);
  }
  
  /* storing the sensitivity vectors */

  fprintf(f1,"---------------------------------- \n");
  fprintf(f1,"Objective \n");
  fprintf(f1,"---------------------------------- \n");
  
  for(i=0;i<neq[1];i++){
     for(iobject=0;iobject<*nobject;iobject++){
        fprintf(f1,"%12.5E",(double)dgdxglob[3+5*i+(5*neq[1]+2)*iobject]);
	fprintf(f1,";  ");    
     }
     fprintf(f1,"\n"); 
  }
  
  fclose(f1);

  /* writing the nominal values in the nom-file for optimizer */
        	
  strcpy2(nominal,jobnamec,132);
  strcat(nominal,".nom");
  
  if((f1=fopen(nominal,"w"))==NULL){
      printf(" *ERROR in sensitivity_out: cannot open sensitivity vector file for writing...");
      
      exit(0);
  }
  
  /* storing the sensitivity vectors */

  fprintf(f1,"---------------------------------- \n");
  fprintf(f1,"Objective \n");
  fprintf(f1,"---------------------------------- \n");
  
  for(iobject=0;iobject<*nobject;iobject++){
     fprintf(f1,"%12.5E",(double)g0[iobject]);
     fprintf(f1,";  ");    
  }
  fprintf(f1,"\n"); 
  
  fclose(f1);
  
  return;

}
