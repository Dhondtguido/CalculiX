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

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void frdselect(double *field1,double *field2,ITG *iset,ITG *nkcoords,ITG *inum,
     char *m1,ITG *istartset,ITG *iendset,ITG *ialset,ITG *ngraph,ITG *ncomp,
     ITG *ifield,ITG *icomp,ITG *nfield,ITG *iselect,char *m2,FILE *f1,
     char *output, char*m3){

  /* storing scalars, components of vectors and tensors without additional
     transformations */

  /* number of components in field1: nfield[0]
     number of components in field2: nfield[1]

     number of entities to store: ncomp
     for each entity i, 0<=i<ncomp:
         - ifield[i]: 1=field1,2=field2
         - icomp[i]: component: 0...,(nfield[0]-1 or nfield[1]-1) */
 
  ITG i,j,k,l,m,n,nksegment,ioutall=0;
      
  int iw;

  float fl;

  if(strcmp1(&output[3],"a")==0) ioutall=1;
  
  if(*iset==0){
    for(i=0;i<*nkcoords;i++){

      /* check whether output is requested for solid nodes or
         network nodes */

      if(ioutall==0){
	  if(*iselect==1){
	      if(inum[i]<=0) continue;
	  }else if(*iselect==-1){
	      if(inum[i]>=0) continue;
	  }else{
	      if(inum[i]==0) continue;
	  }
      }else{
	  if(*iselect==1){
	      if(inum[i]<0) continue;
	  }else if(*iselect==-1){
	      if(inum[i]>0) continue;
	  }
      }

      /* storing the entities */

	for(n=1;n<=(ITG)((*ncomp+5)/6);n++){
	  if(n==1){
	    if(strcmp1(output,"asc")==0){
	      fprintf(f1,"%3s%10" ITGFORMAT "",m1,i+1);
	    }else{
	      iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    }
	    for(j=0;j<min(6,*ncomp);j++){
	      if(ifield[j]==1){
		if(strcmp1(output,"asc")==0){
		  fprintf(f1,"%12.5E",(float)field1[i*nfield[0]+icomp[j]]);
		}else if(strcmp1(output,"bin")==0){
		  fl=(float)field1[i*nfield[0]+icomp[j]];
		  fwrite(&fl,sizeof(float),1,f1);
		}else{
		  fwrite(&field1[i*nfield[0]+icomp[j]],sizeof(double),1,f1);
		}
	      }else{
		if(strcmp1(output,"asc")==0){
		  fprintf(f1,"%12.5E",(float)field2[i*nfield[1]+icomp[j]]);
		}else if(strcmp1(output,"bin")==0){
		  fl=(float)field2[i*nfield[1]+icomp[j]];
		  fwrite(&fl,sizeof(float),1,f1);
		}else{
		  fwrite(&field2[i*nfield[1]+icomp[j]],sizeof(double),1,f1);
		}
	      }
	    }
	    if(strcmp1(output,"asc")==0)fprintf(f1,"\n");
	  }else{
	    if(strcmp1(output,"asc")==0)fprintf(f1,"%3s          ",m2);
	    for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
	      if(ifield[j]==1){
		if(strcmp1(output,"asc")==0){
		  fprintf(f1,"%12.5E",(float)field1[i*nfield[0]+icomp[j]]);
		}else if(strcmp1(output,"bin")==0){
		  fl=(float)field1[i*nfield[0]+icomp[j]];
		  fwrite(&fl,sizeof(float),1,f1);
		}else{
		  fwrite(&field1[i*nfield[0]+icomp[j]],sizeof(double),1,f1);
		}
	      }else{
		if(strcmp1(output,"asc")==0){
		  fprintf(f1,"%12.5E",(float)field2[i*nfield[1]+icomp[j]]);
		}else if(strcmp1(output,"bin")==0){
		  fl=(float)field2[i*nfield[1]+icomp[j]];
		  fwrite(&fl,sizeof(float),1,f1);
		}else{
		  fwrite(&field2[i*nfield[1]+icomp[j]],sizeof(double),1,f1);
		}
	      }
	    }
	    if(strcmp1(output,"asc")==0)fprintf(f1,"\n");
	  }
	}

    }
  }else{
    nksegment=(*nkcoords)/(*ngraph);
    for(k=istartset[*iset-1]-1;k<iendset[*iset-1];k++){
      if(ialset[k]>0){
	for(l=0;l<*ngraph;l++){
	  i=ialset[k]+l*nksegment-1;

	  /* check whether output is requested for solid nodes or
	     network nodes */

	  if(*iselect==1){
	    if(inum[i]<=0) continue;
	  }else if(*iselect==-1){
	    if(inum[i]>=0) continue;
	  }else{
	    if(inum[i]==0) continue;
	  }
	  
	  /* storing the entities */

	  for(n=1;n<=(ITG)((*ncomp+5)/6);n++){
	    if(n==1){
	      if(strcmp1(output,"asc")==0){
		fprintf(f1,"%3s%10" ITGFORMAT "",m1,i+1);
	      }else{
		iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	      }
	      for(j=0;j<min(6,*ncomp);j++){
		if(ifield[j]==1){
		  if(strcmp1(output,"asc")==0){
		    fprintf(f1,"%12.5E",(float)field1[i*nfield[0]+icomp[j]]);
		  }else if(strcmp1(output,"bin")==0){
		    fl=(float)field1[i*nfield[0]+icomp[j]];
		    fwrite(&fl,sizeof(float),1,f1);
		  }else{
		    fwrite(&field1[i*nfield[0]+icomp[j]],sizeof(double),1,f1);
		  }
		}else{
		  if(strcmp1(output,"asc")==0){
		    fprintf(f1,"%12.5E",(float)field2[i*nfield[1]+icomp[j]]);
		  }else if(strcmp1(output,"bin")==0){
		    fl=(float)field2[i*nfield[1]+icomp[j]];
		    fwrite(&fl,sizeof(float),1,f1);
		  }else{
		    fwrite(&field2[i*nfield[1]+icomp[j]],sizeof(double),1,f1);
		  }
		}
	      }
	      if(strcmp1(output,"asc")==0)fprintf(f1,"\n");
	    }else{
	      if(strcmp1(output,"asc")==0)fprintf(f1,"%3s          ",m2);
	      for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
		if(ifield[j]==1){
		  if(strcmp1(output,"asc")==0){
		    fprintf(f1,"%12.5E",(float)field1[i*nfield[0]+icomp[j]]);
		  }else if(strcmp1(output,"bin")==0){
		    fl=(float)field1[i*nfield[0]+icomp[j]];
		    fwrite(&fl,sizeof(float),1,f1);
		  }else{
		    fwrite(&field1[i*nfield[0]+icomp[j]],sizeof(double),1,f1);
		  }
		}else{
		  if(strcmp1(output,"asc")==0){
		    fprintf(f1,"%12.5E",(float)field2[i*nfield[1]+icomp[j]]);
		  }else if(strcmp1(output,"bin")==0){
		    fl=(float)field2[i*nfield[1]+icomp[j]];
		    fwrite(&fl,sizeof(float),1,f1);
		  }else{
		    fwrite(&field2[i*nfield[1]+icomp[j]],sizeof(double),1,f1);
		  }
		}
	      }
	      if(strcmp1(output,"asc")==0)fprintf(f1,"\n");
	    }
	  }
	  
	}
      }else{
	l=ialset[k-2];
	do{
	  l-=ialset[k];
	  if(l>=ialset[k-1]) break;
	  for(m=0;m<*ngraph;m++){
	    i=l+m*nksegment-1;
	    
	    /* check whether output is requested for solid nodes or
	       network nodes */
	    
	    if(*iselect==1){
	      if(inum[i]<=0) continue;
	    }else if(*iselect==-1){
	      if(inum[i]>=0) continue;
	    }else{
	      if(inum[i]==0) continue;
	    }
	    
	    /* storing the entities */
	    
	    for(n=1;n<=(ITG)((*ncomp+5)/6);n++){
	      if(n==1){
		if(strcmp1(output,"asc")==0){
		  fprintf(f1,"%3s%10" ITGFORMAT "",m1,i+1);
		}else{
		  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
		}
		for(j=0;j<min(6,*ncomp);j++){
		  if(ifield[j]==1){
		    if(strcmp1(output,"asc")==0){
		      fprintf(f1,"%12.5E",(float)field1[i*nfield[0]+icomp[j]]);
		    }else if(strcmp1(output,"bin")==0){
		      fl=(float)field1[i*nfield[0]+icomp[j]];
		      fwrite(&fl,sizeof(float),1,f1);
		    }else{
		      fwrite(&field1[i*nfield[0]+icomp[j]],sizeof(double),1,f1);
		    }
		  }else{
		    if(strcmp1(output,"asc")==0){
		      fprintf(f1,"%12.5E",(float)field2[i*nfield[1]+icomp[j]]);
		    }else if(strcmp1(output,"bin")==0){
		      fl=(float)field2[i*nfield[1]+icomp[j]];
		      fwrite(&fl,sizeof(float),1,f1);
		    }else{
		      fwrite(&field2[i*nfield[1]+icomp[j]],sizeof(double),1,f1);
		    }
		  }
		}
		if(strcmp1(output,"asc")==0)fprintf(f1,"\n");
	      }else{
		if(strcmp1(output,"asc")==0)fprintf(f1,"%3s          ",m2);
		for(j=(n-1)*6;j<min(n*6,*ncomp);j++){
		  if(ifield[j]==1){
		    if(strcmp1(output,"asc")==0){
		      fprintf(f1,"%12.5E",(float)field1[i*nfield[0]+icomp[j]]);
		    }else if(strcmp1(output,"bin")==0){
		      fl=(float)field1[i*nfield[0]+icomp[j]];
		      fwrite(&fl,sizeof(float),1,f1);
		    }else{
		      fwrite(&field1[i*nfield[0]+icomp[j]],sizeof(double),1,f1);
		    }
		  }else{
		    if(strcmp1(output,"asc")==0){
		      fprintf(f1,"%12.5E",(float)field2[i*nfield[1]+icomp[j]]);
		    }else if(strcmp1(output,"bin")==0){
		      fl=(float)field2[i*nfield[1]+icomp[j]];
		      fwrite(&fl,sizeof(float),1,f1);
		    }else{
		      fwrite(&field2[i*nfield[1]+icomp[j]],sizeof(double),1,f1);
		    }
		  }
		}
		if(strcmp1(output,"asc")==0)fprintf(f1,"\n");
	      }
	    }
	    
	  }
	}while(1);
      }
    }
  }
  
  if(strcmp1(output,"asc")==0)fprintf(f1,"%3s\n",m3);

  return;

}

