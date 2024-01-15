/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2023 Guido Dhondt                          */

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
#include <pthread.h>
#include "CalculiX.h"

void mastructread(ITG *ipompc,ITG *nodempc,ITG *nmpc,ITG *nactdof,
		  ITG *jq,ITG **mast1p,ITG *neq,ITG *ipointer, ITG *nzs_, 
		  ITG *nmethod,ITG *iperturb,ITG *mi,ITG **nextp,
		  ITG *ifree,ITG *i,ITG *ielmat,char *matname){

  /* reading stiffness and/or mass matrix of substructure
     (superelement) */
  
  FILE *f1;

  char filestiff[81]=" ";
  
  ITG id,index,jdof1,jdof2,idof1,idof2,mpc1,mpc2,id1,id2,ist1,ist2,
    index1,index2,ist,*mast1=NULL,icolumn,mt=mi[1]+1,*next=NULL,
    imat,node1,node2,j,k,m;

  double val;

  mast1=*mast1p;next=*nextp;
  imat=ielmat[*i*mi[2]];
  
  strcpy1(filestiff,&matname[80*(imat-1)],80);
  filestiff[80]=' ';
  for(j=0;j<81;j++){
    if(strcmp1(&filestiff[j]," ")==0){
      filestiff[j]='\0';
      break;
    }
  }
  
  if((f1=fopen(filestiff,"r"))==NULL){
    printf(" *ERROR in frd: cannot open substructure stiffness file for reading...");
    exit(0);
  }
  do{
    if(fscanf(f1,"%d,%d,%d,%d,%lf\n",&node1,&k,&node2,&m,&val)==5){
      mastructmatrix(ipompc,nodempc,nmpc,nactdof,jq,&mast1,neq,ipointer,nzs_,
		     nmethod,iperturb,mi,&next,&node1,&k,&node2,&m,ifree);
      //      printf("%d,%d,%d,%d,%e\n",node1,k,node2,m,val);
    }else{
      break;
    }
  }while(1);
  
  *mast1p=mast1;*nextp=next;
  
  return;
}
