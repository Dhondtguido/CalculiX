/*     CalculiX - A 3-Dimensional finite element program                   */
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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

void writeheading(char *jobnamec,char *heading,ITG *nheading_){
    
    /* writes the headers in the frd-file */
    
    FILE *f1;

    char p1[6]="    1",fneig[132]="",c[2]="C",
         text[67]="                                                                  ";

    ITG i;
    
    strcpy2(fneig,jobnamec,132);
    strcat(fneig,".frd");
    
    if((f1=fopen(fneig,"ab"))==NULL){
	printf(" *ERROR in frd: cannot open frd file for writing...");
	exit(0);
    }

    /* first line */

    fprintf(f1,"%5s%1s\n",p1,c);

    /* header lines */

    for(i=0;i<*nheading_;i++){
	strcpy1(text,&heading[66*i],66);
	fprintf(f1,"%5sU%66s\n",p1,text);
    }

    fclose(f1);
    
    return;
}
