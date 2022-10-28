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

void mastructem(ITG *nk, ITG *kon, ITG *ipkon, char *lakon, ITG *ne,
	      ITG *nodeboun, ITG *ndirboun, ITG *nboun, ITG *ipompc,
	      ITG *nodempc, ITG *nmpc, ITG *nactdof, ITG *icol,
	      ITG *jq, ITG **mast1p, ITG **irowp, ITG *isolver, ITG *neq,
	      ITG *ikmpc, ITG *ilmpc,ITG *ipointer, ITG *nzs, 
	      ITG *ithermal,ITG *mi,ITG *ielmat, double *elcon, ITG *ncmat_, 
	      ITG *ntmat_,ITG *inomat,ITG *network){

  /* determines the structure of the thermo-electromagnetic matrices;
     (i.e. the location of the nonzeros */

  char lakonl[2]=" \0";

  ITG i,j,k,l,jj,ll,id,index,jdof1,jdof2,idof1,idof2,mpc1,mpc2,id1,id2,
    ist1,ist2,node1,node2,isubtract,nmast,ifree,istart,istartold,
    index1,index2,m,node,nzs_,ist,kflag,indexe,nope,isize,*mast1=NULL,
    *irow=NULL,mt=mi[1]+1,imat,idomain,jmin,jmax,*next=NULL,jstart;

  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irow=*irowp;

  kflag=1;
  nzs_=nzs[1];
  NNEW(next,ITG,nzs_);

  /* initialisation of nactmpc */

  for(i=0;i<mt**nk;++i){nactdof[i]=0;}

  /* determining the mechanical active degrees of freedom due to elements */
  if((*ithermal<2)||(*ithermal>=3)){
      for(i=0;i<*ne;++i){
	  
	  if(ipkon[i]<0) continue;
	  indexe=ipkon[i];
	  if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
	  else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
	  else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
	  else if (strcmp1(&lakon[8*i+3],"4")==0) nope=4;
	  else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
	  else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
	  else continue;
	  
          /* degrees of freedom:
             in domain 1: phi
             in domain 2: A and V
             in domain 3: A  */

	  imat=ielmat[i*mi[2]];
	  idomain=(ITG)elcon[(*ncmat_+1)**ntmat_*(imat-1)+2];
	  if(idomain==1){
	      jmin=5;jmax=6;
	  }else if(idomain==2){
	      jmin=1;jmax=5;
	  }else if(idomain==3){
	      jmin=1;jmax=4;
	  }else{
	      continue;
	  }

          /* displacement degrees of freedom */

	  for(j=0;j<nope;++j){
	      node=kon[indexe+j]-1;
	      for(k=jmin;k<jmax;++k){
		  nactdof[mt*node+k]=1;
	      }
	  }

      }
  }

  /* determining the thermal active degrees of freedom due to elements */

  if(*ithermal>1){
      for(i=0;i<*ne;++i){
	  
	  if(ipkon[i]<0) continue;

          /* only the A-V domain (domain 2) has temperature variables */

	  imat=ielmat[i*mi[2]];
	  idomain=(ITG)elcon[(*ncmat_+1)**ntmat_*(imat-1)+2];
	  if(idomain!=2) continue;

	  indexe=ipkon[i];
	  if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
	  else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
	  else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
	  else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
	  else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
	  else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
	  else if (strcmp1(&lakon[8*i],"E")==0){
	      lakonl[0]=lakon[8*i+7];
	      nope=atoi(lakonl)+1;}
	  else if ((strcmp1(&lakon[8*i],"D ")==0)||
                   ((strcmp1(&lakon[8*i],"D")==0)&&(*network==1))){

	      /* check for entry or exit element */
	      
	      if((kon[indexe]==0)||(kon[indexe+2]==0)) continue;
	      
	      /* generic network element */
	      
	      for(j=0;j<3;j=j+2){
		  node=kon[indexe+j]-1;
		  nactdof[mt*node]=1;
	      }
	      continue;}
	  else continue;
	  
	  for(j=0;j<nope;++j){
	      node=kon[indexe+j]-1;
	      nactdof[mt*node]=1;
	  }
      }
  }

  /* determining the active degrees of freedom due to mpc's */

  for(i=0;i<*nmpc;++i){
      index=ipompc[i]-1;
      do{
	  node=nodempc[3*index];
	  if(inomat[node-1]>0){
	      nactdof[mt*(node-1)+nodempc[3*index+1]]=1;
	  }
	  index=nodempc[3*index+2];
	  if(index==0) break;
	  index--;
      }while(1);
  }
  
  /* subtracting the SPC and MPC nodes */

  for(i=0;i<*nboun;++i){
      if(ndirboun[i]>mi[1]){continue;}
      nactdof[mt*(nodeboun[i]-1)+ndirboun[i]]=-2*(i+1);
  }

  for(i=0;i<*nmpc;++i){
    index=ipompc[i]-1;
    if(nodempc[3*index+1]>mi[1]) continue;
    nactdof[mt*(nodempc[3*index]-1)+nodempc[3*index+1]]=-2*i-1;
  }
 
  /* numbering the active degrees of freedom */
  
  neq[0]=0;
  for(i=0;i<*nk;++i){
    for(j=1;j<mt;++j){
	if(nactdof[mt*i+j]>0){
        if((*ithermal<2)||(*ithermal>=3)){
          ++neq[0];
          nactdof[mt*i+j]=neq[0];
        }
        else{
	    nactdof[mt*i+j]=0;
        }
      }
    }
  }
  neq[1]=neq[0];
  for(i=0;i<*nk;++i){
      if(nactdof[mt*i]>0){
      if(*ithermal>1){
        ++neq[1];
        nactdof[mt*i]=neq[1];
      }
      else{
	  nactdof[mt*i]=0;
      }
    }
  }
  neq[2]=neq[1];
  
  ifree=0;

    /* determining the position of each nonzero matrix element in
       the SUPERdiagonal matrix */

    /*   mast1(ipointer(i)) = first nonzero row in column i
	 next(ipointer(i))  points to further nonzero elements in 
                             column i */
      
    for(i=0;i<4**nk;++i){ipointer[i]=0;}

    /* electromagnetic entries */
    
    if((*ithermal<2)||(*ithermal>=3)){

    for(i=0;i<*ne;++i){
      
      if(ipkon[i]<0) continue;
      indexe=ipkon[i];
      if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
      else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
      else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
      else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
      else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
      else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
      else continue;
      
      for(jj=0;jj<mi[1]*nope;++jj){
	
	j=jj/mi[1];
	k=jj-mi[1]*j;
	
	node1=kon[indexe+j];
	jdof1=nactdof[mt*(node1-1)+k+1];
	
	for(ll=jj;ll<mi[1]*nope;++ll){
	  
	  l=ll/mi[1];
	  m=ll-mi[1]*l;
	  
	  node2=kon[indexe+l];
	  jdof2=nactdof[mt*(node2-1)+m+1];
	  
	  /* check whether one of the DOF belongs to a SPC or MPC */
	  
	  if((jdof1>0)&&(jdof2>0)){
	    insert(ipointer,&mast1,&next,&jdof1,&jdof2,&ifree,&nzs_);
	  }
	  else if((jdof1>0)||(jdof2>0)){
	    
	    /* idof1: genuine DOF
	       idof2: nominal DOF of the SPC/MPC */
	    
	    if(jdof1<=0){
	      idof1=jdof2;
	      idof2=jdof1;}
//	      idof2=8*node1+k-7;}
	    else{
	      idof1=jdof1;
	      idof2=jdof2;}
//	      idof2=8*node2+m-7;}
	    
	    if(*nmpc>0){
	      
//	      FORTRAN(nident,(ikmpc,&idof2,nmpc,&id));
//	      if((id>0)&&(ikmpc[id-1]==idof2)){
	      if(idof2!=2*(idof2/2)){
		
		/* regular DOF / MPC */
		
//		id=ilmpc[id-1];
		id=(-idof2+1)/2;
		ist=ipompc[id-1];
		index=nodempc[3*ist-1];
		if(index==0) continue;
		while(1){
		    idof2=nactdof[mt*(nodempc[3*index-3]-1)+nodempc[3*index-2]];
		  if(idof2>0){
		    insert(ipointer,&mast1,&next,&idof1,&idof2,&ifree,&nzs_);
		  }
		  index=nodempc[3*index-1];
		  if(index==0) break;
		}
		continue;
	      }
	    }
	  }
	  
	  else{
//	    idof1=8*node1+k-7;
//	    idof2=8*node2+m-7;
	    idof1=jdof1;
	    idof2=jdof2;
	    mpc1=0;
	    mpc2=0;
	    if(*nmpc>0){
/*	      FORTRAN(nident,(ikmpc,&idof1,nmpc,&id1));
	      if((id1>0)&&(ikmpc[id1-1]==idof1)) mpc1=1;
	      FORTRAN(nident,(ikmpc,&idof2,nmpc,&id2));
	      if((id2>0)&&(ikmpc[id2-1]==idof2)) mpc2=1;*/
	      if(idof1!=2*(idof1/2)) mpc1=1;
	      if(idof2!=2*(idof2/2)) mpc2=1;
	    }
	    if((mpc1==1)&&(mpc2==1)){
//	      id1=ilmpc[id1-1];
//	      id2=ilmpc[id2-1];
	      id1=(-idof1+1)/2;
	      id2=(-idof2+1)/2;
	      if(id1==id2){
		
		/* MPC id1 / MPC id1 */
		
		ist=ipompc[id1-1];
		index1=nodempc[3*ist-1];
		if(index1==0) continue;
		while(1){
		    idof1=nactdof[mt*(nodempc[3*index1-3]-1)+nodempc[3*index1-2]];
		  index2=index1;
		  while(1){
		      idof2=nactdof[mt*(nodempc[3*index2-3]-1)+nodempc[3*index2-2]];
		    if((idof1>0)&&(idof2>0)){
		      insert(ipointer,&mast1,&next,&idof1,&idof2,&ifree,&nzs_);}
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
		    idof1=nactdof[mt*(nodempc[3*index1-3]-1)+nodempc[3*index1-2]];
		  ist2=ipompc[id2-1];
		  index2=nodempc[3*ist2-1];
		  if(index2==0){
		    index1=nodempc[3*index1-1];
		    if(index1==0){break;}
		    else{continue;}
		  }
		  while(1){
		      idof2=nactdof[mt*(nodempc[3*index2-3]-1)+nodempc[3*index2-2]];
		    if((idof1>0)&&(idof2>0)){
		      insert(ipointer,&mast1,&next,&idof1,&idof2,&ifree,&nzs_);}
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

    }

    /* thermal entries*/

    if(*ithermal>1){

    for(i=0;i<*ne;++i){
      
      if(ipkon[i]<0) continue;
      
      /* only the A-V domain (domain 2) has temperature variables */
      
      imat=ielmat[i*mi[2]];
      idomain=(ITG)elcon[(*ncmat_+1)**ntmat_*(imat-1)+2];
      if(idomain!=2) continue;

      indexe=ipkon[i];
      if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
      else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
      else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
      else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
      else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
      else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
      else if (strcmp1(&lakon[8*i],"E")==0){
	  lakonl[0]=lakon[8*i+7];
	  nope=atoi(lakonl)+1;}
      else if ((strcmp1(&lakon[8*i],"D ")==0)||
               ((strcmp1(&lakon[8*i],"D")==0)&&(*network==1))){

	  /* check for entry or exit element */

	  if((kon[indexe]==0)||(kon[indexe+2]==0)) continue;
	  nope=3;}
      else continue;
      
      for(jj=0;jj<nope;++jj){
	
	j=jj;
	
	node1=kon[indexe+j];
	jdof1=nactdof[mt*(node1-1)];
	
	for(ll=jj;ll<nope;++ll){
	  
	  l=ll;
	  
	  node2=kon[indexe+l];
	  jdof2=nactdof[mt*(node2-1)];
	  
	  /* check whether one of the DOF belongs to a SPC or MPC */
	  
	  if((jdof1>0)&&(jdof2>0)){
	    insert(ipointer,&mast1,&next,&jdof1,&jdof2,&ifree,&nzs_);
	  }
	  else if((jdof1>0)||(jdof2>0)){
	    
	    /* idof1: genuine DOF
	       idof2: nominal DOF of the SPC/MPC */
	    
	    if(jdof1<=0){
	      idof1=jdof2;
//	      idof2=8*node1-8;}
	      idof2=jdof1;}
	    else{
	      idof1=jdof1;
//	      idof2=8*node2-8;}
	      idof2=jdof2;}
	    
	    if(*nmpc>0){
	      
//	      FORTRAN(nident,(ikmpc,&idof2,nmpc,&id));
//	      if((id>0)&&(ikmpc[id-1]==idof2)){
	      if(idof2!=2*(idof2/2)){
		
		/* regular DOF / MPC */
		
//		id=ilmpc[id-1];
		id=(-idof2+1)/2;
		ist=ipompc[id-1];
		index=nodempc[3*ist-1];
		if(index==0) continue;
		while(1){
		    idof2=nactdof[mt*(nodempc[3*index-3]-1)+nodempc[3*index-2]];
		  if(idof2>0){
		    insert(ipointer,&mast1,&next,&idof1,&idof2,&ifree,&nzs_);
		  }
		  index=nodempc[3*index-1];
		  if(index==0) break;
		}
		continue;
	      }
	    }
	  }
	  
	  else{
//	    idof1=8*node1-8;
//	    idof2=8*node2-8;
	    idof1=jdof1;
	    idof2=jdof2;
	    mpc1=0;
	    mpc2=0;
	    if(*nmpc>0){
/*	      FORTRAN(nident,(ikmpc,&idof1,nmpc,&id1));
	      if((id1>0)&&(ikmpc[id1-1]==idof1)) mpc1=1;
	      FORTRAN(nident,(ikmpc,&idof2,nmpc,&id2));
	      if((id2>0)&&(ikmpc[id2-1]==idof2)) mpc2=1;*/
		if(idof1!=2*(idof1/2)) mpc1=1;
		if(idof2!=2*(idof2/2)) mpc2=1;
	    }
	    if((mpc1==1)&&(mpc2==1)){
//	      id1=ilmpc[id1-1];
//	      id2=ilmpc[id2-1];
	      id1=(-idof1+1)/2;
	      id2=(-idof2+1)/2;
	      if(id1==id2){
		
		/* MPC id1 / MPC id1 */
		
		ist=ipompc[id1-1];
		index1=nodempc[3*ist-1];
		if(index1==0) continue;
		while(1){
		    idof1=nactdof[mt*(nodempc[3*index1-3]-1)+nodempc[3*index1-2]];
		  index2=index1;
		  while(1){
		      idof2=nactdof[mt*(nodempc[3*index2-3]-1)+nodempc[3*index2-2]];
		    if((idof1>0)&&(idof2>0)){
		      insert(ipointer,&mast1,&next,&idof1,&idof2,&ifree,&nzs_);}
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
		    idof1=nactdof[mt*(nodempc[3*index1-3]-1)+nodempc[3*index1-2]];
		  ist2=ipompc[id2-1];
		  index2=nodempc[3*ist2-1];
		  if(index2==0){
		    index1=nodempc[3*index1-1];
		    if(index1==0){break;}
		    else{continue;}
		  }
		  while(1){
		      idof2=nactdof[mt*(nodempc[3*index2-3]-1)+nodempc[3*index2-2]];
		    if((idof1>0)&&(idof2>0)){
		      insert(ipointer,&mast1,&next,&idof1,&idof2,&ifree,&nzs_);}
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

    }

    if(neq[1]==0){
      printf("\n *WARNING: no degrees of freedom in the model\n\n");
    }

    /*   determination of the following fields:       

       - irow: row numbers, column per column
       - icol(i)=# SUBdiagonal nonzero's in column i
       - jq(i)= location in field irow of the first SUBdiagonal
         nonzero in column i  */

    RENEW(irow,ITG,ifree);
    nmast=0;
    jq[0]=1;
    for(i=0;i<neq[1];i++){
	index=ipointer[i];
	do{
	    if(index==0) break;
	    irow[nmast++]=mast1[index-1];
	    index=next[index-1];
	}while(1);
	jq[i+1]=nmast+1;
//	icol[i]=jq[i+1]-jq[i];
    }

    /* sorting the row numbers within each column */

    for(i=0;i<neq[1];++i){
      if(jq[i+1]-jq[i]>0){
	isize=jq[i+1]-jq[i];
	FORTRAN(isortii,(&irow[jq[i]-1],&mast1[jq[i]-1],&isize,&kflag));
      }
    }

    /* removing duplicate entries */

    nmast=0;
    for(i=0;i<neq[1];i++){
	jstart=nmast+1;
	if(jq[i+1]-jq[i]>0){
	    irow[nmast++]=irow[jq[i]-1];
	    for(j=jq[i];j<jq[i+1]-1;j++){
		if(irow[j]==irow[nmast-1])continue;
		irow[nmast++]=irow[j];
	    }
	}
	jq[i]=jstart;
    }
    jq[neq[1]]=nmast+1;

    for(i=0;i<neq[1];i++){
	icol[i]=jq[i+1]-jq[i];
    }

    if(neq[0]==0){nzs[0]=0;}
    else{nzs[0]=jq[neq[0]]-1;}
    nzs[1]=jq[neq[1]]-1;

    nzs[2]=nzs[1];

    /* summary */

    printf(" number of equations\n");
    printf(" %" ITGFORMAT "\n",neq[1]);
    printf(" number of nonzero lower triangular matrix elements\n");
    printf(" %" ITGFORMAT "\n",nmast);
    printf("\n");
  
    SFREE(next);

    *mast1p=mast1;
    *irowp=irow;

    return;

}
