/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

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
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>

#include "CalculiX.h"
#include "readfrd.h"

void utempread (double *t1,ITG *istep,char *jobnamec)
{
 
    char  datin[MAX_LINE_LENGTH],text[13]="            ";
    Summen    anz[1]; 
    Nodes     *node=NULL;
    Elements  *elem=NULL;
    Datasets *lcase=NULL;

    ITG i,j,read_mode=0,loadcase,istep_global,nodenr;
    
    /* reading the global coordinates and the topology from file
       (if any, else return) */
    
    if(strcmp1(&jobnamec[660]," ")==0)return;
    strcpy1(datin,&jobnamec[660],132); 
    for(i=0;i<MAX_LINE_LENGTH;i++){
	if(strcmp1(&datin[i]," ")==0){
	    datin[i]='\0';
	    break;
	}
    }
    
    /* initialization of the size of fields used in readfrd.c */

    anz->orign=0;
    anz->n=0;
    anz->e=0;
    anz->f=0;
    anz->g=0;
    anz->t=0;
    anz->l=0;
    anz->olc=0;
    anz->orignmax=0;
    anz->nmax=0;
    anz->nmin=MAX_INTEGER;
    anz->emax=0;
    anz->emin=MAX_INTEGER;
    anz->sets=0;
    anz->mats=0;
    anz->amps=0;
    anz->nnext=0;
    anz->enext=0;
    
    readfrd( datin, anz, &node, &elem, &lcase, read_mode);

    /* check for the existence of nodes and/or elements */

    if((anz[0].n==0)||(anz[0].e==0)){
	printf(" *ERROR in utempread: there are either no nodes or\n no elements or neither nodes nor elements in the temperature frd-file\n");
	FORTRAN(stop,());
    }
    
    /* loading the step data : NDTEMP (1 variable) if present */
    
    /* reading the temperatures */
    /* 1. determining the appropriate temperature loadcase in the step */
    
    loadcase=-1;
    for(i=0;i<anz[0].l;i++){
	for(j=0;j<lcase[i].npheader;j++){
	    if(strcmp1(&lcase[i].pheader[j][5],"PSTEP")==0){
		strcpy1(text,&lcase[i].pheader[j][48],12);
		istep_global=atoi(text);
		break;
	    }
	}
	if((istep_global==*istep)&&
	   (strcmp1(lcase[i].name,"NDTEMP")==0)){
	    loadcase=i;
	}else if(istep_global>*istep){
	    break;
	}
    }
    
    /* 2. reading the data */
    
    if(loadcase>-1){
	if(!read_mode && readfrdblock(loadcase, anz, node, lcase )==-1) 
	{
	    printf(" *ERROR in utempread: Could not read data for Dataset:%" ITGFORMAT "\n", i+1); 
	    FORTRAN(stop,());
	}
	
    /* 3. storing the data */
    
	for(i=0;i<anz[0].n;i++){
	    nodenr=node[i].nr;
	    t1[nodenr-1]=lcase[loadcase].dat[0][nodenr];
	}
    }else{
	printf("INFO in utempread: no temperature data\n was found for step %d in the temperature frd-file\n\n",*istep);
    }
    
    for(j=0;j<anz->l;j++){
      freeDatasets(lcase,j);
    }
    SFREE(lcase);lcase=NULL;
    
    return;

}
