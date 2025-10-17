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
#include <locale.h>
#include "CalculiX.h"

void frdgeneralvector(double *v,ITG *iset,ITG *ntrans,char * filabl,
               ITG *nkcoords,
               ITG *inum,char *m1,ITG *inotr,double *trab,double *co,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *mi,ITG *ngraph,
	       FILE *f1,char *output,char *m3,ITG *ioutall){

    ITG i,j,k,l,nksegment,*inodeset=NULL,nnodeset,nnodeset_,
    *iy=NULL,kflag;
      
  int iw;

  float fl;

  setlocale(LC_NUMERIC, "C");

  if(*iset==0){
    if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)){
      for(i=0;i<*nkcoords;i++){
	if(*ioutall==0){
	  if(inum[i]<=0) continue;
	}else{
	  if(inum[i]<0) continue;
	}
	if(strcmp1(output,"asc")==0){
	    if(mi[1]==4){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4]);
	    }else if(mi[1]==5){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			(float)v[(mi[1]+1)*i+5]);
	    }else if(mi[1]==6){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			(float)v[(mi[1]+1)*i+5],(float)v[(mi[1]+1)*i+6]);
	    }
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  if(strcmp1(output,"bin")==0){
	      for(j=1;j<=mi[1];j++){
		  fl=(float)v[(mi[1]+1)*i+j];fwrite(&fl,sizeof(float),1,f1);
	      }
	  }else{
	      for(j=1;j<=mi[1];j++){
		  fwrite(&v[(mi[1]+1)*i+j],sizeof(double),1,f1);
	      }
	  }
	}
      }
    }else{
      for(i=0;i<*nkcoords;i++){
	if(*ioutall==0){
	  if(inum[i]<=0) continue;
	}else{
	  if(inum[i]<0) continue;
	}
	if(inotr[2*i]==0){
	  if(strcmp1(output,"asc")==0){
	    if(mi[1]==4){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4]);
	    }else if(mi[1]==5){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			(float)v[(mi[1]+1)*i+5]);
	    }else if(mi[1]==6){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			(float)v[(mi[1]+1)*i+5],(float)v[(mi[1]+1)*i+6]);
	    }
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    if(strcmp1(output,"bin")==0){
		for(j=1;j<=mi[1];j++){
		    fl=(float)v[(mi[1]+1)*i+j];fwrite(&fl,sizeof(float),1,f1);
		}
	    }else{
		for(j=1;j<=mi[1];j++){
		    fwrite(&v[(mi[1]+1)*i+j],sizeof(double),1,f1);
		}
	    }
	  }
	}else{
	    printf(" *WARNING in frdgeneralvector:\n");
	    printf("          no output in local coordinates allowed\n");
	    printf("          for generalized vectors\n");
	    printf("          output request ist not performed;\n");
	}
      }
    }
  }else{

    /* output for a node set */

    nnodeset_=100;
    nnodeset=0;
    NNEW(inodeset,ITG,nnodeset_);

    /* collecting all nodes from the set */
    
    for(k=istartset[*iset-1]-1;k<iendset[*iset-1];k++){
      if(ialset[k]>0){
	nnodeset++;
	if(nnodeset>nnodeset_){
	  nnodeset_=(ITG)(1.1*nnodeset_);
	  RENEW(inodeset,ITG,nnodeset_);
	}
	inodeset[nnodeset-1]=ialset[k];
      }else{
	l=ialset[k-2];
	do{
	  l-=ialset[k];
	  if(l>=ialset[k-1]) break;
	  nnodeset++;
	  if(nnodeset>nnodeset_){
	    nnodeset_=(ITG)(1.1*nnodeset_);
	    RENEW(inodeset,ITG,nnodeset_);
	  }
	  inodeset[nnodeset-1]=l;
	}while(1);
      }
    }
	
    /* sorting the nodes in ascending order (required by the
       frd format */

    kflag=1;
    FORTRAN(isortii,(inodeset,iy,&nnodeset,&kflag));

    /* storing the results */
    
    nksegment=(*nkcoords)/(*ngraph);
    for(l=0;l<*ngraph;l++){
      for(k=0;k<nnodeset;k++){
	i=inodeset[k]+l*nksegment-1;
	if(*ioutall==0){
	  if(inum[i]<=0) continue;
	}else{
	  if(inum[i]<0) continue;
	}
	if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)||(inotr[2*i]==0)){
	  if(strcmp1(output,"asc")==0){
	    if(mi[1]==4){
	      fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
		      (float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
		      (float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4]);
	    }else if(mi[1]==5){
	      fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
		      (float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
		      (float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
		      (float)v[(mi[1]+1)*i+5]);
	    }else if(mi[1]==6){
	      fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
		      (float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
		      (float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
		      (float)v[(mi[1]+1)*i+5],(float)v[(mi[1]+1)*i+6]);
	    }
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    if(strcmp1(output,"bin")==0){
	      for(j=1;j<=mi[1];j++){
		fl=(float)v[(mi[1]+1)*i+j];fwrite(&fl,sizeof(float),1,f1);
	      }
	    }else{
	      for(j=1;j<=mi[1];j++){
		fwrite(&v[(mi[1]+1)*i+j],sizeof(double),1,f1);
	      }
	    }
	  }
	}else{
	  printf(" *WARNING in frdgeneralvector:\n");
	  printf("          no output in local coordinates allowed\n");
	  printf("          for generalized vectors\n");
	  printf("          output request ist not performed;\n");
	}
      }
    }
  }
      
  if(strcmp1(output,"asc")==0)fprintf(f1,"%3s\n",m3);

  SFREE(inodeset);

  return;

}

