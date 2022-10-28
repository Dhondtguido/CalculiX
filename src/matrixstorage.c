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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"
#include "matrixstorage.h"

void matrixstorage(double *ad, double **aup, double *adb, double *aub, 
		   double *sigma,ITG *icol, ITG **irowp, 
		   ITG *neq, ITG *nzs, ITG *ntrans, ITG *inotr,
		   double *trab, double *co, ITG *nk, ITG *nactdof,
		   char *jobnamec, ITG *mi, ITG *ipkon, char *lakon,
		   ITG *kon, ITG *ne, ITG *mei, ITG *nboun, ITG *nmpc,
		   double *cs, ITG *mcs,ITG *ithermal,ITG *nmethod){

  char fsti[132]="",fmas[132]="",fdof[132]="";
  ITG i,j,l,*irow=NULL,*ai=NULL,*aj=NULL,kflag=2,ndim,jref,kstart,klen,
    npoint_,npoint,neq3,index,i3l,i3c,i3lo,i3co,idof,n,il,
    ic,id,itrans,ndim2,*ipoindex=NULL,mt=mi[1]+1,*nactdofinv=NULL,
    *nodorig=NULL,inode,idir;
  long long *ipoint=NULL,k;
  double *au=NULL,*aa=NULL,*trans=NULL,*aa3=NULL,a[9];
  FILE *f2,*f3,*f4;

  /* stiffness or heat conduction matrix */

  if(*ithermal!=2){
    strcpy2(fsti,jobnamec,132);
    strcat(fsti,".sti");
  }else{
    strcpy2(fsti,jobnamec,132);
    strcat(fsti,".con");
  }

  printf(" Storing the stiffness matrix in file %s \n\n",fsti);

  if(*mcs!=0){
    if(cs[1]>=0){
      printf(" For cyclic symmetric calculations the complex\n");
      printf(" Hermitian matrix is stored as a symmetric real\n");
      printf(" matrix double its size; if R stands for the real\n");
      printf(" part of the matrix and I for the imaginary part,\n");
      printf(" the resulting matrix takes the form:\n");
      printf("  _        _\n");
      printf(" |          |\n");
      printf(" |  R   -I  |\n");
      printf(" |  I    R  |\n");
      printf(" |_        _|\n\n");
      printf(" This applies to the stiffness and the mas matrix\n\n");
    }
  }

  if((f2=fopen(fsti,"wb"))==NULL){
    printf(" *ERROR in matrixstorage: cannot open %s for writing...\n",fsti);
    FORTRAN(stop,());
  }

  au=*aup;
  irow=*irowp;

  ndim=*neq+*nzs;

  itrans=0;
  if(*ntrans!=0){
    for(i=0;i<*nk;i++){
      if(inotr[2*i]!=0){
        itrans=1;
        break;
      }
    }
  }

  /* stiffness matrix */

  if((itrans==0)||(mei[0]==1)){
    
    /* no transformation */

    NNEW(aa,double,ndim);
    NNEW(ai,ITG,ndim);
    NNEW(aj,ITG,ndim);
    
    k=0;
    for(i=0;i<*neq;i++){
      ai[k]=i+1;
      aj[k]=i+1;
      aa[k]=ad[i];
      k++;
    }
    l=0;
    for(i=0;i<*neq;i++){
      for(j=0;j<icol[i];j++){
        ai[k]=i+1;
        aj[k]=irow[l];
        aa[k]=au[l];
        k++;l++;
      }
    }
  }
  else{
    
    /* transformation: storing the matrices in transformed coordinates 
       (needed by turboreduce (not part of CalculiX)) */

    if((*nboun!=0)||(*nmpc!=0)){
      printf(" *ERROR in matrixstorage: matrix storage in local\n");
      printf("        coordinates is only possible in the absence\n");
      printf("        of SPC's and MPC's\n\n");
      FORTRAN(stop,());
    }

    ndim2=*neq+2**nzs;

    /* aa contains the linear storage of the individual matrix elements */

    NNEW(aa,double,ndim2);

    /* aa3 contains the linear storage of the 3x3 matrices */

    NNEW(aa3,double,ndim2);
    NNEW(ai,ITG,ndim2);
    NNEW(aj,ITG,ndim2);
    
    k=0;
    for(i=0;i<*neq;i++){
      ai[k]=i+1;
      aj[k]=i+1;
      aa[k]=ad[i];
      k++;
    }
    l=0;
    for(i=0;i<*neq;i++){
      for(j=0;j<icol[i];j++){
        ai[k]=i+1;
        aj[k]=irow[l];
        aa[k]=au[l];
        k++;
        ai[k]=irow[l];
        aj[k]=i+1;
        aa[k]=au[l];
        k++;
        l++;
      }
    }

    FORTRAN(isortiid,(aj,ai,aa,&ndim2,&kflag));

    k=0;
    for(i=0;i<*neq;i++){
      jref=aj[k];
      kstart=k;
      do{
        k++;
	if(k==ndim2) break;
        if(aj[k]!=jref) break;
      }while(1);
      klen=k-kstart;
      FORTRAN(isortiid,(&ai[kstart],&aj[kstart],&aa[kstart],&klen,&kflag));
    }

    /* npoint is the number of 3x3 matrices
       ipoint contains the two-dimensional location of the matrices 
       (row and column stored as column x neq3 + row)
       ipoindex contains the one-dimensional location of the matrices in
       fields aa3 */

    npoint_=*neq;
    npoint=0;
    NNEW(ipoint,long long,npoint_);
    NNEW(ipoindex,ITG,npoint_);

    neq3=*neq/3;
    index=0;

    /* storing the matrix as a matrix of 3x3 submatrices */

    for(i=0;i<ndim2;i++){
      il=ai[i];
      ic=aj[i];
      i3co=(ic-1)%3;
      i3c=((ic-1)-i3co)/3;
      i3lo=(il-1)%3;
      i3l=((il-1)-i3lo)/3;
      k=(long long)i3c*neq3+i3l;
      FORTRAN(nidentll,(ipoint,&k,&npoint,&id));
      if(npoint==0){
        npoint++;
        ipoint[npoint-1]=k;
        ipoindex[npoint-1]=index;
      }
      else if(ipoint[id-1]!=k){
        npoint++;
        if(npoint>npoint_){
          npoint_=(ITG)(1.1*npoint_);
          RENEW(ipoint,long long,npoint_);
          RENEW(ipoindex,ITG,npoint_);
        }
        index+=9;
        ipoint[npoint-1]=k;
        ipoindex[npoint-1]=index;
      }
      else{
        index=ipoindex[id-1];
      }
      aa3[index+3*i3co+i3lo]=aa[i];
    }

    /* defining the transformation matrix (diagonal matrix of
       3x3 submatrices */

    NNEW(trans,double,9*neq3);
    for (i=0;i<*nk;i++){
      idof=nactdof[mt*i+1];
      if(idof<=0) continue;
      itrans=inotr[2*i];
      if(itrans==0){
        for(j=0;j<9;j++){
          trans[3*(idof-1)+j]=0.;
        }
        trans[3*(idof-1)]=1.;
        trans[3*(idof-1)+4]=1.;
        trans[3*(idof-1)+8]=1.;
      }
      else{
        FORTRAN(transformatrix,(&trab[7*itrans-7],&co[3*i],a));
        for(j=0;j<9;j++){
          trans[3*(idof-1)+j]=a[j];
        }
      }
    }

    /* postmultiplying the matrix with the transpose of the
       transformation matrix */

    for(i=0;i<npoint;i++){
      i3l=ipoint[i]%neq3;
      i3c=(ipoint[i]-i3l)/neq3;
      n=2;
      FORTRAN(mult,(&aa3[9*i],&trans[9*i3c],&n));
    }

    /* premultiplying the matrix with the transformation matrix */

    for(i=0;i<npoint;i++){
      i3l=ipoint[i]%neq3;
      n=1;
      FORTRAN(mult,(&aa3[9*i],&trans[9*i3l],&n));
    }

    /* storing the new matrix in a linear format */

    k=-1;
    for(i=0;i<npoint;i++){
      i3l=ipoint[i]%neq3;
      i3c=(ipoint[i]-i3l)/neq3;
      for(j=0;j<9;j++){
        i3lo=j%3;
        i3co=(j-i3lo)/3;
        ic=3*i3c+i3co+1;
        il=3*i3l+i3lo+1;
        if(il>ic) continue;
        k++;
        ai[k]=il;
        aj[k]=ic;
        aa[k]=aa3[9*i+j];
      }
    }
    SFREE(aa3);SFREE(ipoint);SFREE(ipoindex);SFREE(trans);
  }

  FORTRAN(isortiid,(aj,ai,aa,&ndim,&kflag));

  k=0;
  for(i=0;i<*neq;i++){
    jref=aj[k];
    kstart=k;
    do{
      k++;
      if(k==ndim) break;
      if(aj[k]!=jref) break;
    }while(1);
    klen=k-kstart;
    FORTRAN(isortiid,(&ai[kstart],&aj[kstart],&aa[kstart],&klen,&kflag));
  }

  for(i=0;i<ndim;i++){
    fprintf(f2,"%" ITGFORMAT " %" ITGFORMAT " %20.13e\n",ai[i],aj[i],aa[i]);
  }

  fclose(f2);

  SFREE(ai);SFREE(aj);SFREE(aa);

  /* mass or capacity (specific heat) matrix or stress stiffness (buckling) */

  /* exclude steady state heat calculations */

  if(*nmethod!=1){

    if(*nmethod!=3){
      if(*ithermal!=2){
	strcpy2(fmas,jobnamec,132);
	strcat(fmas,".mas");
      }else{
	strcpy2(fmas,jobnamec,132);
	strcat(fmas,".sph");
      }
    }else{
      strcpy2(fmas,jobnamec,132);
      strcat(fmas,".str");
    }
      
    printf(" Storing the mass matrix in file %s \n\n",fmas);
      
    if((f3=fopen(fmas,"wb"))==NULL){
      printf(" *ERROR in matrixstorage: cannot open %s for writing...\n",fmas);
      FORTRAN(stop,());
    }
      
    if((itrans==0)||(mei[0]==1)){
	  
      /* no transformation */
	  
      NNEW(aa,double,ndim);
      NNEW(ai,ITG,ndim);
      NNEW(aj,ITG,ndim);
	  
      k=0;
      for(i=0;i<*neq;i++){
	ai[k]=i+1;
	aj[k]=i+1;
	aa[k]=adb[i];
	k++;
      }
      l=0;
      for(i=0;i<*neq;i++){
	for(j=0;j<icol[i];j++){
	  ai[k]=i+1;
	  aj[k]=irow[l];
	  aa[k]=aub[l];
	  k++;l++;
	}
      }
    }
    else{
	  
	  
      /* transformation: storing the matrices in transformed coordinates 
	 (needed by turboreduce (not part of CalculiX)) */
	  
      ndim2=*neq+2**nzs;
	  
      /* aa contains the linear storage of the individual matrix elements */
	  
      NNEW(aa,double,ndim2);
	  
      /* aa3 contains the linear storage of the 3x3 matrices */
	  
      NNEW(aa3,double,ndim2);
      NNEW(ai,ITG,ndim2);
      NNEW(aj,ITG,ndim2);
	  
      k=0;
      for(i=0;i<*neq;i++){
	ai[k]=i+1;
	aj[k]=i+1;
	aa[k]=adb[i];
	k++;
      }
      l=0;
      for(i=0;i<*neq;i++){
	for(j=0;j<icol[i];j++){
	  ai[k]=i+1;
	  aj[k]=irow[l];
	  aa[k]=aub[l];
	  k++;
	  ai[k]=irow[l];
	  aj[k]=i+1;
	  aa[k]=aub[l];
	  k++;
	  l++;
	}
      }
	  
      FORTRAN(isortiid,(aj,ai,aa,&ndim2,&kflag));
	  
      k=0;
      for(i=0;i<*neq;i++){
	jref=aj[k];
	kstart=k;
	do{
	  k++;
	  if(k==ndim2) break;
	  if(aj[k]!=jref) break;
	}while(1);
	klen=k-kstart;
	FORTRAN(isortiid,(&ai[kstart],&aj[kstart],&aa[kstart],&klen,&kflag));
      }
	  
      /* npoint is the number of 3x3 matrices
	 ipoint contains the two-dimensional location of the matrices 
	 (row and column stored as column x neq3 + row)
	 ipoindex contains the one-dimensional location of the matrices in
	 fields aa3 */
	  
      npoint_=*neq;
      npoint=0;
      NNEW(ipoint,long long,npoint_);
      NNEW(ipoindex,ITG,npoint_);
	  
      neq3=*neq/3;
      index=0;
	  
      /* storing the matrix as a matrix of 3x3 submatrices */
	  
      for(i=0;i<ndim2;i++){
	il=ai[i];
	ic=aj[i];
	i3co=(ic-1)%3;
	i3c=((ic-1)-i3co)/3;
	i3lo=(il-1)%3;
	i3l=((il-1)-i3lo)/3;
	k=(long long)i3c*neq3+i3l;
	FORTRAN(nidentll,(ipoint,&k,&npoint,&id));
	if(npoint==0){
	  npoint++;
	  ipoint[npoint-1]=k;
	  ipoindex[npoint-1]=index;
	}
	else if(ipoint[id-1]!=k){
	  npoint++;
	  if(npoint>npoint_){
	    npoint_=(ITG)(1.1*npoint_);
	    RENEW(ipoint,long long,npoint_);
	    RENEW(ipoindex,ITG,npoint_);
	  }
	  index+=9;
	  ipoint[npoint-1]=k;
	  ipoindex[npoint-1]=index;
	}
	else{
	  index=ipoindex[id-1];
	}
	aa3[index+3*i3co+i3lo]=aa[i];
      }
	  
      /* defining the transformation matrix (diagonal matrix of
	 3x3 submatrices */
	  
      NNEW(trans,double,9*neq3);
      for (i=0;i<*nk;i++){
	idof=nactdof[mt*i+1];
	if(idof<=0) continue;
	itrans=inotr[2*i];
	if(itrans==0){
	  for(j=0;j<9;j++){
	    trans[3*(idof-1)+j]=0.;
	  }
	  trans[3*(idof-1)]=1.;
	  trans[3*(idof-1)+4]=1.;
	  trans[3*(idof-1)+8]=1.;
	}
	else{
	  FORTRAN(transformatrix,(&trab[7*itrans-7],&co[3*i],a));
	  for(j=0;j<9;j++){
	    trans[3*(idof-1)+j]=a[j];
	  }
	}
      }
	  
      /* postmultiplying the matrix with the transpose of the
	 transformation matrix */
	  
      for(i=0;i<npoint;i++){
	i3l=ipoint[i]%neq3;
	i3c=(ipoint[i]-i3l)/neq3;
	n=2;
	FORTRAN(mult,(&aa3[9*i],&trans[9*i3c],&n));
      }
	  
      /* premultiplying the matrix with the transformation matrix */
	  
      for(i=0;i<npoint;i++){
	i3l=ipoint[i]%neq3;
	n=1;
	FORTRAN(mult,(&aa3[9*i],&trans[9*i3l],&n));
      }
	  
      /* storing the new matrix in a linear format */
	  
      k=-1;
      for(i=0;i<npoint;i++){
	i3l=ipoint[i]%neq3;
	i3c=(ipoint[i]-i3l)/neq3;
	for(j=0;j<9;j++){
	  i3lo=j%3;
	  i3co=(j-i3lo)/3;
	  ic=3*i3c+i3co+1;
	  il=3*i3l+i3lo+1;
	  if(il>ic) continue;
	  k++;
	  ai[k]=il;
	  aj[k]=ic;
	  aa[k]=aa3[9*i+j];
	}
      }
      SFREE(aa3);SFREE(ipoint);SFREE(ipoindex);SFREE(trans);
    }
      
    FORTRAN(isortiid,(aj,ai,aa,&ndim,&kflag));
      
    k=0;
    for(i=0;i<*neq;i++){
      jref=aj[k];
      kstart=k;
      do{
	k++;
	if(k==ndim) break;
	if(aj[k]!=jref) break;
      }while(1);
      klen=k-kstart;
      FORTRAN(isortiid,(&ai[kstart],&aj[kstart],&aa[kstart],&klen,&kflag));
    }
      
    for(i=0;i<ndim;i++){
      fprintf(f3,"%" ITGFORMAT " %" ITGFORMAT " %20.13e\n",ai[i],aj[i],aa[i]);
    }
      
    fclose(f3);
      
    SFREE(ai);SFREE(aj);SFREE(aa);
      
  }  /* end mass matrix */

  *aup=au;
  *irowp=irow;

  /* node and global direction for each degree of freedom in the
     matrix (corresponding with a row or column number) */

  strcpy2(fdof,jobnamec,132);
  strcat(fdof,".dof");

  if((itrans==0)||(mei[0]==1)){
    printf(" Storing the node and global direction per entry (row or column)\n in the stiffness and mass matrices in the form node.direction in file %s \n\n",fdof);
  }else{
    printf(" Storing the node and local direction per entry (row or column)\n in the stiffness and mass matrices in the form node.direction in file %s \n\n",fdof);
  }

  if((f4=fopen(fdof,"wb"))==NULL){
    printf(" *ERROR in matrixstorage: cannot open %s for writing...\n",fdof);
    FORTRAN(stop,());
  }

  /* invert nactdof */

  NNEW(nactdofinv,ITG,mt**nk);
  NNEW(nodorig,ITG,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
  SFREE(nodorig);
  
  if((*mcs==0)||(cs[1]<0)){
    for(i=0;i<*neq;i++){
      inode=(ITG)((double)nactdofinv[(ITG)i]/mt)+1;
      idir=nactdofinv[(ITG)i]-mt*(inode-1);
      fprintf(f4,"%" ITGFORMAT ".%" ITGFORMAT "\n",inode,idir);
    }
  }else{
    for(i=0;i<*neq/2;i++){
      inode=(ITG)((double)nactdofinv[(ITG)i]/mt)+1;
      idir=nactdofinv[(ITG)i]-mt*(inode-1);
      fprintf(f4,"%" ITGFORMAT ".%" ITGFORMAT "\n",inode,idir);
    }
  }
  
  fclose(f4);

  printf("\n Job finished\n\n");
  
  //  FORTRAN(stopwithout201,());

  return;
}
