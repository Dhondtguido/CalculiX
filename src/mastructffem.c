/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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

void mastructffem(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
		  ITG *nodeboun,ITG *ndirboun,ITG *nboun,ITG *ipompc,
		  ITG *nodempc,ITG *nmpc,ITG *nactdoh,
		  ITG *icolv,ITG *icolp,ITG *jqv,ITG *jqp,
		  ITG **mast1p,ITG **irowvp,ITG **irowpp,
		  ITG *neqp,ITG *ipointer,ITG *nzsv,ITG *nzsp,
		  ITG *nzs,ITG *compressible,ITG *inomat){

  ITG i,j,k,l,jj,ll,id,index,jdof1,jdof2,idof1,idof2,mpc1,mpc2,id1,id2,
    ist1,ist2,node1,node2,isubtract,nmast,ifree,istart,istartold,idir,
    index1,index2,node,nzs_,ist,kflag,indexe,nope,isize,*mast1=NULL,
    *irowv=NULL,*irowp=NULL,imaterial;

  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irowv=*irowvp;irowp=*irowpp;

  kflag=2;
  
  /* initialisation of nactdoh */
  
  for(i=0;i<*nk;++i){nactdoh[i]=0;}
  
  /* determining the active degrees of freedom due to elements */
  
  for(i=0;i<*ne;++i){
      
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"F")!=0) continue;
    indexe=ipkon[i];
    if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
    else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
    else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
    else continue;
      
    for(j=0;j<nope;++j){
      nactdoh[kon[indexe+j]-1]=1;
    }
  }
  
  /* determining the active degrees of freedom due to mpc's */
  
  for(i=0;i<*nmpc;++i){

    index=ipompc[i]-1;
    do{
      node=nodempc[3*index]-1;

      /* if node belongs to an element: mark the material number */

      if(nactdoh[node]!=0){
	imaterial=inomat[node];
	break;
      }
      index=nodempc[3*index+2];
      if(index==0) break;
      index--;
    }while(1);

    /* assign this material to all nodes of the MPC */
    
    index=ipompc[i]-1;
    do{
      node=nodempc[3*index]-1;
      idir=nodempc[3*index+1];
      inomat[node]=imaterial;
      index=nodempc[3*index+2];
      if(index==0) break;
      index--;
    }while(1);
  }
  
  /* subtracting pressure SPCs (only for incompressible fluids) */

  if(*compressible==0){
    for(i=0;i<*nboun;++i){
      if(ndirboun[i]!=4){continue;}
      nactdoh[nodeboun[i]-1]=-2*(i+1);
    }
  }      
  
  /* subtracting pressure MPCs (only for incompressible fluids) */
      
  if(*compressible==0){
    for(i=0;i<*nmpc;++i){
      index=ipompc[i]-1;
      if(nodempc[3*index+1]!=4) continue;
      nactdoh[nodempc[3*index]-1]=-2*i-1;
    }
  }      
  
  /* numbering the active degrees of freedom */
  
  *neqp=0;
  for(i=0;i<*nk;++i){
    if(nactdoh[i]>0){
      ++(*neqp);
      nactdoh[i]=*neqp;
    }
  }
  printf("neqppp=%d\n",*neqp);
  
  /* generic mass (temperature/velocity/compressible pressure/turbulent) entries*/

    ifree=0;
    nzs_=*nzs;
      
    /* determining the position of each nonzero matrix element
	 
       mast1(ipointer(i)) = first nonzero row in column i
       irow(ipointer(i))  points to further nonzero elements in 
       column i */
      
    for(i=0;i<*nk;++i){ipointer[i]=0;}
      
    for(i=0;i<*ne;++i){
	  
      indexe=ipkon[i];
      
      if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
      else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
      else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
      else continue;
	  
      for(j=0;j<nope;++j){
	jdof1=kon[indexe+j];
	for(l=j;l<nope;++l){
	  jdof2=kon[indexe+l];
	  insertcbs(ipointer,&mast1,&irowv,&jdof1,&jdof2,&ifree,&nzs_);
	}
      }
    }
      
    for(i=0;i<*nk;++i){
      if(ipointer[i]==0){
	printf(" *ERROR in mastructffem: zero column\n");
	FORTRAN(stop,());
      }
      istart=ipointer[i];
      while(1){
	istartold=istart;
	istart=irowv[istart-1];
	irowv[istartold-1]=i+1;
	if(istart==0) break;
      }
    }
      
    /* defining icolv and jqv */
      
    if(*nk==0){
      printf("\n *WARNING in mastructffem: no degrees of freedom in the generic mass matrix\n\n");
    }
      
    nmast=ifree;
      
    /* summary */
      
    printf(" number of generic mass equations\n");
    printf(" %d\n",*nk);
    printf(" number of nonzero generic mass matrix elements\n");
    printf(" %d\n",nmast);
    printf("\n");
      
    /* changing the meaning of icolv,jqv,mast1,irowv:
	 
       - irowv is going to contain the row numbers of the SUBdiagonal
       nonzero's, column per column
       - mast1 contains the column numbers
       - icolv(i)=# SUBdiagonal nonzero's in column i
       - jqv(i)= location in field irow of the first SUBdiagonal
       nonzero in column i
	 
    */
      
    /* switching from a SUPERdiagonal inventory to a SUBdiagonal one */
      
    FORTRAN(isortii,(mast1,irowv,&nmast,&kflag));
      
    /* filtering out the diagonal elements and generating icolv and jqv */
      
    isubtract=0;
    for(i=0;i<*nk;++i){icolv[i]=0;}
    k=0;
    for(i=0;i<nmast;++i){
      if(mast1[i]==irowv[i]){++isubtract;}
      else{
	mast1[i-isubtract]=mast1[i];
	irowv[i-isubtract]=irowv[i];
	if(k!=mast1[i]){
	  for(l=k;l<mast1[i];++l){jqv[l]=i+1-isubtract;}
	  k=mast1[i];
	}
	++icolv[k-1];
      }
    }
    nmast=nmast-isubtract;
    for(l=k;l<*nk+1;++l){jqv[l]=nmast+1;}
      
    for(i=0;i<*nk;++i){
      if(jqv[i+1]-jqv[i]>0){
	isize=jqv[i+1]-jqv[i];
	FORTRAN(isortii,(&irowv[jqv[i]-1],&mast1[jqv[i]-1],&isize,&kflag));
      }
    }
      
    if(*nk==0){*nzsv=0;}
    else{*nzsv=jqv[*nk]-1;}
  
  /* pressure entries */
  
  ifree=0;
  nzs_=*nzs;
  RENEW(mast1,ITG,nzs_);
  for(i=0;i<nzs_;i++){mast1[i]=0;}
  
  /* determining the position of each nonzero matrix element
     
     mast1(ipointer(i)) = first nonzero row in column i
     irowp(ipointer(i))  points to further nonzero elements in 
     column i */
  
  for(i=0;i<3**nk;++i){ipointer[i]=0;}
  
  for(i=0;i<*ne;++i){
      
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"F")!=0) continue;
    indexe=ipkon[i];
    if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
    else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
    else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
    else continue;
      
    for(jj=0;jj<nope;++jj){
	  
      j=jj;
	  
      node1=kon[indexe+j];
      jdof1=nactdoh[node1-1];
	  
      for(ll=jj;ll<nope;++ll){
	      
	l=ll;
	      
	node2=kon[indexe+l];
	jdof2=nactdoh[node2-1];
	      
	/* check whether one of the DOF belongs to a SPC or MPC */
	      
	if((jdof1>0)&&(jdof2>0)){
	  insertcbs(ipointer,&mast1,&irowp,&jdof1,&jdof2,&ifree,&nzs_);
	}
	else if((jdof1>0)||(jdof2>0)){
		  
	  /* idof1: genuine DOF
	     idof2: nominal DOF of the SPC/MPC */
		  
	  if(jdof1<=0){
	    idof1=jdof2;
	    idof2=jdof1;}
	  else{
	    idof1=jdof1;
	    idof2=jdof2;}
		  
	  if(*nmpc>0){


	    if(idof2!=2*(idof2/2)){
			  
	      /* regular DOF / MPC */
			  
	      id=(-idof2+1)/2;
	      ist=ipompc[id-1];
	      index=nodempc[3*ist-1];
	      if(index==0) continue;
	      while(1){
		idof2=nactdoh[nodempc[3*index-3]-1];
		if(idof2>0){
		  insertcbs(ipointer,&mast1,&irowp,&idof1,&idof2,&ifree,&nzs_);
		}
		index=nodempc[3*index-1];
		if(index==0) break;
	      }
	      continue;
	    }
	  }
	}
	      
	else{
	  idof1=jdof1;
	  idof2=jdof2;
	  mpc1=0;
	  mpc2=0;
	  if(*nmpc>0){
	    if(idof1!=2*(idof1/2)) mpc1=1;
	    if(idof2!=2*(idof2/2)) mpc2=1;
	  }
	  if((mpc1==1)&&(mpc2==1)){
	    id1=(-idof1+1)/2;
	    id2=(-idof2+1)/2;
	    if(id1==id2){
			  
	      /* MPC id1 / MPC id1 */
			  
	      ist=ipompc[id1-1];
	      index1=nodempc[3*ist-1];
	      if(index1==0) continue;
	      while(1){
		idof1=nactdoh[nodempc[3*index1-3]-1];
		index2=index1;
		while(1){
		  idof2=nactdoh[nodempc[3*index2-3]-1];
		  if((idof1>0)&&(idof2>0)){
		    insertcbs(ipointer,&mast1,&irowp,&idof1,&idof2,&ifree,&nzs_);}
		  index2=nodempc[3*index2-1];
		  if(index2==0) break;
		}
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	      }
	    }
		      
	    else{
			  
	      /* MPC id1 /MPC id2 */
			  
	      ist1=ipompc[id1-1];
	      index1=nodempc[3*ist1-1];
	      if(index1==0) continue;
	      while(1){
		idof1=nactdoh[nodempc[3*index1-3]-1];
		ist2=ipompc[id2-1];
		index2=nodempc[3*ist2-1];
		if(index2==0){
		  index1=nodempc[3*index1-1];
		  if(index1==0){break;}
		  else{continue;}
		}
		while(1){
		  idof2=nactdoh[nodempc[3*index2-3]-1];
		  if((idof1>0)&&(idof2>0)){
		    insertcbs(ipointer,&mast1,&irowp,&idof1,&idof2,&ifree,&nzs_);}
		  index2=nodempc[3*index2-1];
		  if(index2==0) break;
		}
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	      }
	    }
	  }
	}
      }
    }
  }
    
  for(i=0;i<*neqp;++i){
    if(ipointer[i]==0){
      if(i>=*neqp) continue;
      printf(" *ERROR in mastructffem: zero column\n");
      FORTRAN(stop,());
    }
    istart=ipointer[i];
    while(1){
      istartold=istart;
      istart=irowp[istart-1];
      irowp[istartold-1]=i+1;
      if(istart==0) break;
    }
  }
  
  /* defining icolp and jqp */
  
  if(*neqp==0){
    printf("\n*WARNING in matructf: no degrees of freedom in the pressure matrix\n\n");
  }
  
  nmast=ifree;
  
  /* summary */
  
  printf(" number of pressure equations\n");
  printf(" %d\n",*neqp);
  printf(" number of nonzero pressure matrix elements\n");
  printf(" %d\n",nmast);
  printf("\n");
  
  /* changing the meaning of icolp,jqp,mast1,irowp:
     
     - irowp is going to contain the row numbers of the SUBdiagonal
     nonzero's, column per column
     - mast1 contains the column numbers
     - icolp(i)=# SUBdiagonal nonzero's in column i
     - jqp(i)= location in field irow of the first SUBdiagonal
     nonzero in column i
  
  */
  
  /* switching from a SUPERdiagonal inventory to a SUBdiagonal one */
  
  FORTRAN(isortii,(mast1,irowp,&nmast,&kflag));
  
  /* filtering out the diagonal elements and generating icolp and jqp */
  
  isubtract=0;
  for(i=0;i<*neqp;++i){icolp[i]=0;}
  k=0;
  for(i=0;i<nmast;++i){
    if(mast1[i]==irowp[i]){++isubtract;}
    else{
      mast1[i-isubtract]=mast1[i];
      irowp[i-isubtract]=irowp[i];
      if(k!=mast1[i]){
	for(l=k;l<mast1[i];++l){jqp[l]=i+1-isubtract;}
	k=mast1[i];
      }
      ++icolp[k-1];
    }
  }
  nmast=nmast-isubtract;
  for(l=k;l<*neqp+1;++l){jqp[l]=nmast+1;}
  
  for(i=0;i<*neqp;++i){
    if(jqp[i+1]-jqp[i]>0){
      isize=jqp[i+1]-jqp[i];
      FORTRAN(isortii,(&irowp[jqp[i]-1],&mast1[jqp[i]-1],&isize,&kflag));
    }
  }
  
  if(*neqp==0){*nzsp=0;}
  else{*nzsp=jqp[*neqp]-1;}

  *mast1p=mast1;
  *irowvp=irowv;*irowpp=irowp;

  return;

}
