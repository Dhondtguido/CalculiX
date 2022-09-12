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
#include <pthread.h>
#include "CalculiX.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

static ITG *neq1,*irow1,*mast11,*jq1,kflag1,num_cpus;

void mastruct(ITG *nk, ITG *kon, ITG *ipkon, char *lakon, ITG *ne,
	      ITG *nodeboun, ITG *ndirboun, ITG *nboun, ITG *ipompc,
	      ITG *nodempc, ITG *nmpc, ITG *nactdof, ITG *icol,
	      ITG *jq, ITG **mast1p, ITG **irowp, ITG *isolver, ITG *neq,
	      ITG *ikmpc, ITG *ilmpc,ITG *ipointer, ITG *nzs, 
              ITG *nmethod,ITG *ithermal, ITG *ikboun, ITG *ilboun, 
              ITG *iperturb, ITG *mi, ITG *mortar, char *typeboun,
              char *labmpc, ITG *iit, ITG *icascade,ITG *network,
              ITG *iexpl){

  /* determines the structure of the thermo-mechanical matrices;
     (i.e. the location of the nonzeros */

  char lakonl[2]=" \0",lakonl2[3]="  \0";

  ITG i,j,k,l,jj,ll,id,index,jdof1,jdof2,idof1,idof2,mpc1,mpc2,id1,id2,
    ist1,ist2,node1,node2,isubtract,nmast,ifree,istart,istartold,
    index1,index2,m,node,nzs_,ist,kflag,indexe,nope,isize,*mast1=NULL,
    *irow=NULL,icolumn,nmastboun,mt=mi[1]+1,jmax,*next=NULL,nopeold=0,
    indexeold,identical,jstart,iatleastonenonzero,idof,ndof,*ithread=NULL;

  /* variables for multithreading procedure */

  ITG sys_cpus;
  char *env,*envloc,*envsys;
      
  num_cpus=0;
  sys_cpus=0;
  
  /* explicit user declaration prevails */
  
  envsys=getenv("NUMBER_OF_CPUS");
  if(envsys){
    sys_cpus=atoi(envsys);
    if(sys_cpus<0) sys_cpus=0;
  }
  
  /* automatic detection of available number of processors */
  
  if(sys_cpus==0){
    sys_cpus = getSystemCPUs();
    if(sys_cpus<1) sys_cpus=1;
  }
  
  /* local declaration prevails, if strictly positive */
  
  envloc = getenv("CCX_NPROC_CFD");
  if(envloc){
    num_cpus=atoi(envloc);
    if(num_cpus<0){
      num_cpus=0;
    }else if(num_cpus>sys_cpus){
      num_cpus=sys_cpus;
    }
  }
  
  /* else global declaration, if any, applies */
  
  env = getenv("OMP_NUM_THREADS");
  if(num_cpus==0){
    if (env)
      num_cpus = atoi(env);
    if (num_cpus < 1) {
      num_cpus=1;
    }else if(num_cpus>sys_cpus){
      num_cpus=sys_cpus;
    }
  }

  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irow=*irowp;

  kflag=1;

  /* determining nactdof (only at start of step or if MPC's
     changed */

  if((*iit<0)||(*icascade!=0)){

    /* initialisation of nactdof */
      
    for(i=0;i<mt**nk;++i){
      nactdof[i]=0;
    }
      
    /* determining the mechanical active degrees of freedom due to elements */

    if((*ithermal<2)||(*ithermal>=3)){
      for(i=0;i<*ne;++i){
	      
	if(ipkon[i]<0) continue;
	if(strcmp1(&lakon[8*i],"F")==0)continue;
	indexe=ipkon[i];
	if(strcmp1(&lakon[8*i],"U")==0){

	  /* user element
	     number of dofs: 7th entry of label
	     number of nodes: 8th entry of label */

	  ndof=lakon[8*i+6];
	  nope=lakon[8*i+7];}
	/* Bernhardi start */
	else if (strcmp1(&lakon[8*i+3],"8I")==0){nope=11;ndof=3;}
	else if(strcmp1(&lakon[8*i+3],"20")==0){nope=20;ndof=3;}
	/* Bernhardi end */
	else if (strcmp1(&lakon[8*i+3],"8")==0){nope=8;ndof=3;}
	else if (strcmp1(&lakon[8*i+3],"10")==0){nope=10;ndof=3;}
	else if ((strcmp1(&lakon[8*i+3],"4")==0)||
		 (strcmp1(&lakon[8*i+2],"4")==0)){nope=4;ndof=3;}
	else if (strcmp1(&lakon[8*i+3],"15")==0){nope=15;ndof=3;}
	else if (strcmp1(&lakon[8*i+3],"6")==0){nope=6;ndof=3;}
	else if (strcmp1(&lakon[8*i],"E")==0){
	  if((strcmp1(&lakon[8*i+6],"C")==0)&&(*mortar==1)){

	    /* face-to-face contact (all nodes already belong
	       to other elements */

	    continue;
	  }else if(strcmp1(&lakon[8*i+6],"F")!=0){

	    /* node-to-face contact */

	    lakonl[0]=lakon[8*i+7];
	    nope=atoi(lakonl)+1;
	    ndof=3;
	  }else{

	    /* advection elements */

	    continue;
	  }
	}else continue;
	      
	/* displacement degrees of freedom */
	      
	for(j=0;j<nope;++j){
	  node=kon[indexe+j]-1;
	  for(k=1;k<=ndof;++k){
	    nactdof[mt*node+k]=1;
	  }
	}
      }
    }
      
    /* determining the thermal active degrees of freedom due to elements */
      
    if(*ithermal>1){
      for(i=0;i<*ne;++i){
	      
	if(ipkon[i]<0) continue;
	if(strcmp1(&lakon[8*i],"F")==0)continue;
	indexe=ipkon[i];
	if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
	else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
	else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
	else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
	else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
	else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
	else if (strcmp1(&lakon[8*i],"E")==0){
	  if((strcmp1(&lakon[8*i+6],"C")==0)&&(*mortar==1)){
	    continue;
	  }else{
	    lakonl[0]=lakon[8*i+7];
	    nope=atoi(lakonl)+1;
	  }
	}else if ((strcmp1(&lakon[8*i],"D ")==0)||
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
      if (strcmp1(&labmpc[20*i],"FLUID")==0) continue;

      iatleastonenonzero=0;

      index=ipompc[i]-1;
      do{
	if(nodempc[3*index+1]<4){
	  idof=mt*(nodempc[3*index]-1)+nodempc[3*index+1];
	  if(nactdof[idof]==1){
	    iatleastonenonzero=1;
	  }else{
	    nactdof[idof]=1;
	  }
	}
	index=nodempc[3*index+2];
	if(index==0) break;
	index--;
      }while(1);

      if(iatleastonenonzero==1) continue;

      /* if all dofs in the MPC were inactive, keep then inactive
	 (may e.g. belong to network elements) */

      index=ipompc[i]-1;
      do{
	if(nodempc[3*index+1]<4){
	  nactdof[mt*(nodempc[3*index]-1)+nodempc[3*index+1]]=0;
	}
	index=nodempc[3*index+2];
	if(index==0) break;
	index--;
      }while(1);

    }
      
    /* subtracting the SPC and MPC nodes */
      
    for(i=0;i<*nboun;++i){
      if(ndirboun[i]>mi[1]) continue;
      if(strcmp1(&typeboun[i],"F")==0) continue;
      nactdof[mt*(nodeboun[i]-1)+ndirboun[i]]=-2*(i+1);
    }
      
    for(i=0;i<*nmpc;++i){
      if(strcmp1(&labmpc[20*i],"FLUID")==0) continue;
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
    if((*nmethod==2)||((*nmethod==4)&&(*iperturb<=1))||((*nmethod>=5)&&(*nmethod<=7))){
      neq[2]=neq[1]+*nboun;
    }
    else{neq[2]=neq[1];}
  }
  
  /* determining the subdiagonal nonzeros in the stiffness/mass matrix */

  ifree=0;
  nzs_=nzs[1];
  NNEW(next,ITG,nzs_);

  /* determining the position of each nonzero matrix element in
     the SUBdiagonal matrix (excluding diagonal): x-elements on
     the left of the vertical line */

  //               |x x x
  //        x      |x x x
  //        x x    |x x x
  //        x x x  |x x x

  /*   mast1(ipointer(i)) = first nonzero row in column i
       next(ipointer(i))  points to further nonzero elements in 
       column i */
      
  for(i=0;i<4**nk;++i){
    ipointer[i]=0;
  }

  /* mechanical entries */
    
  if((*ithermal<2)||(*ithermal>=3)){

    indexeold=0;

    for(i=0;i<*ne;++i){
      
      if(ipkon[i]<0) continue;
      if(strcmp1(&lakon[8*i],"F")==0)continue;
      indexe=ipkon[i];
      /* Bernhardi start */
      if (strcmp1(&lakon[8*i+3],"8I")==0){nope=11;ndof=3;}
      else if(strcmp1(&lakon[8*i+3],"20")==0){nope=20;ndof=3;}
      /* Bernhardi end */
      else if (strcmp1(&lakon[8*i+3],"8")==0){nope=8;ndof=3;}
      else if (strcmp1(&lakon[8*i+3],"10")==0){nope=10;ndof=3;}
      else if (strcmp1(&lakon[8*i+3],"4")==0){nope=4;ndof=3;}
      else if (strcmp1(&lakon[8*i+3],"15")==0){nope=15;ndof=3;}
      else if (strcmp1(&lakon[8*i+3],"6")==0){nope=6;ndof=3;}
      else if (strcmp1(&lakon[8*i],"E")==0){
	if((strcmp1(&lakon[8*i+6],"C")==0)&&(*mortar==1)){
	  nope=kon[indexe-1];
	  if(nope==nopeold){
	    identical=1;
	    for(j=0;j<nope;j++){
	      if(kon[indexe+j]!=kon[indexeold+j]){
		identical=0;
		break;
	      }
	    }
	    if(identical==1)continue;
	  }
	  nopeold=nope;
	  indexeold=indexe;
	  ndof=3;
	}else{
	  lakonl[0]=lakon[8*i+7];
	  nope=atoi(lakonl)+1;
	  ndof=3;
	}
      }else if(strcmp1(&lakon[8*i],"U")==0){
	  
	/* user element
	   number of dofs: 7th entry of label
	   number of nodes: 8th entry of label */
	  
	ndof=lakon[8*i+6];
	nope=lakon[8*i+7];
      }else continue;
      
      for(jj=0;jj<ndof*nope;++jj){
	
	j=jj/ndof;
	k=jj-ndof*j;
	
	node1=kon[indexe+j];
	jdof1=nactdof[mt*(node1-1)+k+1];
	
	for(ll=jj;ll<ndof*nope;++ll){
	  
	  l=ll/ndof;
	  m=ll-ndof*l;
	  
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

            /* regular DOF/SPC */

            /* boundary stiffness coefficients (for frequency
               and modal dynamic calculations) : x-elements
               on the right of the vertical line */

	    //               |x x x
	    //        x      |x x x
	    //        x x    |x x x
	    //        x x x  |x x x

	    if((*nmethod==2)||((*nmethod==4)&&(*iperturb<=1))||((*nmethod>=5)&&(*nmethod<=7))){
	      icolumn=neq[1]-idof2/2;
	      insertfreq(ipointer,&mast1,&next,&idof1,&icolumn,&ifree,&nzs_);
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

    indexeold=0;

    for(i=0;i<*ne;++i){
      
      if(ipkon[i]<0) continue;
      if(strcmp1(&lakon[8*i],"F")==0)continue;
      indexe=ipkon[i];
      if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
      else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
      else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
      else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
      else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
      else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
      else if (strcmp1(&lakon[8*i],"E")==0){
	if((strcmp1(&lakon[8*i+6],"C")==0)&&(*mortar==1)){
	  nope=kon[indexe-1];
	  if(nope==nopeold){
	    identical=1;
	    for(j=0;j<nope;j++){
	      if(kon[indexe+j]!=kon[indexeold+j]){
		identical=0;
		break;
	      }
	    }
	    if(identical==1)continue;
	  }
	  nopeold=nope;
	  indexeold=indexe;
	}else{
	  lakonl[0]=lakon[8*i+7];
	  nope=atoi(lakonl)+1;
	}
      }else if ((strcmp1(&lakon[8*i],"D ")==0)||
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

            /* regular DOF/SPC */

            /* boundary stiffness coefficients (for frequency and
               modal dynamic calculations */

	    if((*nmethod==2)||((*nmethod==4)&&(*iperturb<=1))||((*nmethod>=5)&&(*nmethod<=7))){
	      icolumn=neq[1]-idof2/2;
	      insertfreq(ipointer,&mast1,&next,&idof1,&icolumn,&ifree,&nzs_);
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

  /* subdiagonal elements of the regular stiffness/mass matrices
     (marked by X underneath) */ 

  //               |x x x
  //        X      |x x x
  //        X X    |x x x
  //        X X X  |x x x

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
  }

  /* sorting the row numbers within each column */
  
  if(neq[1]<num_cpus) num_cpus=neq[1];
  
  printf(" Using up to %d cpu(s) for setting up the structure of the matrix.\n", num_cpus);
  
  pthread_t tid[num_cpus];

  neq1=neq;irow1=irow;mast11=mast1;jq1=jq;kflag1=kflag;
	  
  NNEW(ithread,ITG,num_cpus);
  for(i=0; i<num_cpus; i++)  {
    ithread[i]=i;
    pthread_create(&tid[i], NULL, (void *)sortingmt, (void *)&ithread[i]);
  }
  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
  SFREE(ithread);

  /* for(i=0;i<neq[1];++i){
    if(jq[i+1]-jq[i]>0){
      isize=jq[i+1]-jq[i];
      FORTRAN(isortii,(&irow[jq[i]-1],&mast1[jq[i]-1],&isize,&kflag));
    }
    }*/

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

  if(neq[0]==0){
    nzs[0]=0;}
  else{
    nzs[0]=jq[neq[0]]-1;}
  nzs[1]=jq[neq[1]]-1;

  /* summary */

  if(*iexpl<=1){
    printf(" number of equations\n");
    printf(" %" ITGFORMAT "\n",neq[1]);
    printf(" number of nonzero lower triangular matrix elements\n");
    printf(" %" ITGFORMAT "\n",nmast);
    printf("\n");
  }

  /* determining irow, jq and icol for the boundary stiffness matrix (only
     for frequency and modal dynamic calculations) 
     (entries marked by X underneath) */

  //               |X X X
  //        x      |X X X
  //        x x    |X X X
  //        x x x  |X X X

  if((*nmethod==2)||((*nmethod==4)&&(*iperturb<=1))||((*nmethod>=5)&&(*nmethod<=7))){

    nmastboun=nmast;
    for(i=neq[1];i<neq[2];i++){
      index=ipointer[i];
      do{
	if(index==0) break;
	irow[nmastboun++]=mast1[index-1];
	index=next[index-1];
      }while(1);
      jq[i+1]=nmastboun+1;
    }

    /* sorting the row numbers within each column */
  
    if((neq[2]-neq[1])<num_cpus) num_cpus=neq[2]-neq[1];

    printf(" Using up to %d cpu(s) for setting up the structure of the matrix.\n", num_cpus);
  
    pthread_t tid[num_cpus];

    neq1=neq;irow1=irow;mast11=mast1;jq1=jq;kflag1=kflag;
	  
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)sortingfreqmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    SFREE(ithread);

    /*  for(i=neq[1];i<neq[2];++i){
      if(jq[i+1]-jq[i]>0){
	isize=jq[i+1]-jq[i];
	FORTRAN(isortii,(&irow[jq[i]-1],&mast1[jq[i]-1],&isize,&kflag));
      }
      }*/

    /* removing duplicate entries */

    nmastboun=nmast;
    for(i=neq[1];i<neq[2];i++){
      jstart=nmastboun+1;
      if(jq[i+1]-jq[i]>0){
	irow[nmastboun++]=irow[jq[i]-1];
	for(j=jq[i];j<jq[i+1]-1;j++){
	  if(irow[j]==irow[nmastboun-1])continue;
	  irow[nmastboun++]=irow[j];
	}
      }
      jq[i]=jstart;
    }
    jq[neq[2]]=nmastboun+1;
	
    for(i=neq[1];i<neq[2];i++){
      icol[i]=jq[i+1]-jq[i];
    }

    /* number of nonzero's in the boundary part*/

    nzs[2]=jq[neq[2]]-1;
  }
  else{nzs[2]=nzs[1];}

  SFREE(next);

  *mast1p=mast1;
  *irowp=irow;

  return;

}

/* sorting the row entries in each column */

void *sortingmt(ITG *i){

  ITG neqa,neqb,neqdelta,j,isize;
    
  neqdelta=(ITG)ceil(neq1[1]/(double)num_cpus);
  neqa=*i*neqdelta;
  neqb=neqa+neqdelta;
  if(neqb>neq1[1]) neqb=neq1[1];

  for(j=neqa;j<neqb;++j){
    if(jq1[j+1]-jq1[j]>0){
      isize=jq1[j+1]-jq1[j];
      FORTRAN(isortii,(&irow1[jq1[j]-1],&mast11[jq1[j]-1],&isize,&kflag1));
    }
  }
  
  return NULL;
}

/* sorting the row entries in each column */

void *sortingfreqmt(ITG *i){

  ITG neqa,neqb,neqdelta,j,isize;
    
  neqdelta=(ITG)ceil((neq1[2]-neq1[1])/(double)num_cpus);
  neqa=*i*neqdelta+neq1[1];
  neqb=neqa+neqdelta;
  if(neqb>neq1[2]) neqb=neq1[2];

  for(j=neqa;j<neqb;++j){
    if(jq1[j+1]-jq1[j]>0){
      isize=jq1[j+1]-jq1[j];
      FORTRAN(isortii,(&irow1[jq1[j]-1],&mast11[jq1[j]-1],&isize,&kflag1));
    }
  }
  
  return NULL;
}
