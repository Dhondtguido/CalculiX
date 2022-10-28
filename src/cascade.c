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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef SPOOLES
#include <misc.h>
#include <FrontMtx.h>
#include <SymbFac.h>
#endif

#include "CalculiX.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void cascade(ITG *ipompc, double **coefmpcp, ITG **nodempcp, ITG *nmpc,
	     ITG *mpcfree, ITG *nodeboun, ITG *ndirboun, ITG*nboun, ITG*ikmpc,
	     ITG *ilmpc, ITG *ikboun, ITG *ilboun, ITG *mpcend,
	     char *labmpc, ITG *nk, ITG *memmpc_, ITG *icascade, ITG *maxlenmpc,
	     ITG *callfrommain, ITG *iperturb, ITG *ithermal){

  /*   detects cascaded mpc's and decascades them; checks multiple
       occurrence of the same dependent DOF's in different mpc/spc's

       data structure of ipompc,coefmpc,nodempc:
       for each mpc, e.g. i, 
       -the nodes are stored in nodempc(1,ipompc(i)),
       nodempc(1,nodempc(3,ipompc(i))),
       nodempc(1,nodempc(3,nodempc(3,ipompc(i))))... till
       nodempc(3,nodempc(3,nodempc(3,.......))))))=0;
       -the corresponding directions in nodempc(2,ipompc(i)),
       nodempc(2,nodempc(3,ipompc(i))),.....
       -the corresponding coefficient in coefmpc(ipompc(i)),
       coefmpc(nodempc(3,ipompc(i))),.....
       the mpc is written as a(1)u(i1,j1)+a(2)u(i2,j2)+...
       +....a(k)u(ik,jk)=0, the first term is the dependent term,
       the others are independent, at least after execution of the
       present routine. The mpc's must be homogeneous, otherwise a
       error message is generated and the program stops. */

  ITG i,j,index,id,idof,nterm,idepend,*nodempc=NULL,
    ispooles,iexpand,ichange,indexold,ifluidmpc,
    mpc,indexnew,index1,index2,index1old,index2old,*jmpc=NULL,nl;

  double coef,*coefmpc=NULL;

  nodempc=*nodempcp;
  coefmpc=*coefmpcp;
    
  /*     for(i=0;i<*nmpc;i++){
	  j=i+1;
	  FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
	  }*/

  NNEW(jmpc,ITG,*nmpc);
  idepend=0;

  /*        check whether a node is used as a dependent node in a MPC
	    and in a SPC */

  for(i=0;i<*nmpc;i++){
    if(*nboun>0){
      FORTRAN(nident,(ikboun,&ikmpc[i],nboun,&id));}
    else{id=0;}
    if(id>0){
      if(ikboun[id-1]==ikmpc[i]){
	if(strcmp1(&labmpc[20*i],"FLUID")!=0){
	  printf(" *ERROR in cascade: the DOF corresponding to \n node %" ITGFORMAT " in direction %" ITGFORMAT " is detected on the \n dependent side of a MPC and a SPC\n\n",
		 (ikmpc[i])/8+1,ikmpc[i]-8*((ikmpc[i])/8));
	}else{
	  printf(" *ERROR in cascade: the DOF corresponding to \n face %" ITGFORMAT " in direction %" ITGFORMAT " is detected on the \n dependent side of a MPC and a SPC\n\n",
		 (-ikmpc[i])/8+1,-ikmpc[i]-8*((-ikmpc[i])/8));
	}
	FORTRAN(stop,());
      }
    }
  }

  /*     check whether there are user mpc's: in user MPC's the
	 dependent DOF can change, however, the number of terms
	 cannot change   */

  for(i=0;i<*nmpc;i++){

    /* linear mpc */

    /* because of the next line the size of field labmpc
       has to be defined as 20*nmpc+1: without "+1" an
       undefined field is accessed */

    if((strcmp1(&labmpc[20*i],"                    ")==0) ||
       (strcmp1(&labmpc[20*i],"CYCLIC")==0) ||
       (strcmp1(&labmpc[20*i],"SUBCYCLIC")==0)||
       (strcmp1(&labmpc[20*i],"PRETENSION")==0)||
       (strcmp1(&labmpc[20*i],"THERMALPRET")==0)||
       //	   (strcmp1(&labmpc[20*i],"CONTACT")==0)||
       (strcmp1(&labmpc[20*i],"FLUID")==0)||
       (*iperturb<2)) jmpc[i]=0;

    /* nonlinear mpc */

    else if((strcmp1(&labmpc[20*i],"RIGID")==0) ||
	    (strcmp1(&labmpc[20*i],"KNOT")==0) ||
	    (strcmp1(&labmpc[20*i],"PLANE")==0) ||
	    (strcmp1(&labmpc[20*i],"BEAM")==0) ||
	    (strcmp1(&labmpc[20*i],"STRAIGHT")==0)) jmpc[i]=1;

    /* user mpc */

    else{
      jmpc[i]=1;
      if(*icascade==0) *icascade=1;
    }
  }
    
  /*     decascading */

  ispooles=0;

  /* decascading using simple substitution */

  do{
    ichange=0;
    for(i=0;i<*nmpc;i++){

      if(strcmp1(&labmpc[20*i],"FLUID")!=0){
	ifluidmpc=0;
      }else{
	ifluidmpc=1;
      }

      if(jmpc[i]==1) nl=1;
      else nl=0;
      iexpand=0;
      index=nodempc[3*ipompc[i]-1];
      if(index==0) continue;
      do{
	if(ifluidmpc==0){
	  /* MPC on node */
	  idof=(nodempc[3*index-3]-1)*8+nodempc[3*index-2];
	}else{
	  if(nodempc[3*index-3]>0){
	    /* MPC on face */
	    idof=-((nodempc[3*index-3]-1)*8+nodempc[3*index-2]);
	  }else{
	    /* MPC on node 
	       SPC number: -nodempc[3*index-3]
	       node: nodeboun[-nodempc[3*index-3]-1] */

	    idof=(nodeboun[-nodempc[3*index-3]-1]-1)*8+nodempc[3*index-2];
	  }
	}
		
	FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
	if((id>0)&&(ikmpc[id-1]==idof)){
		    
	  /* a term on the independent side of the MPC is
	     detected as dependent node in another MPC */

	  //		    indexold=nodempc[3*index-1];
	  //		    coef=coefmpc[index-1];
	  mpc=ilmpc[id-1];

	  /* no expansion if there is a dependence of a
	     nonlinear MPC on another linear or nonlinear MPC
	     and the call is from main */ 

	  if((jmpc[mpc-1]==1)||(nl==1)){
	    *icascade=2;
	    if(idepend==0){
	      printf(" *INFO in cascade: linear MPCs and\n");
	      printf("       nonlinear MPCs depend on each other\n");
	      printf("       common node: %" ITGFORMAT " in direction %" ITGFORMAT "\n\n",nodempc[3*index-3],nodempc[3*index-2]);
	      idepend=1;}
	    if(*callfrommain==1){
	      index=nodempc[3*index-1];
	      if(index!=0) continue;
	      else break;}
	  }

	  /* collecting terms corresponding to the same DOF */
		    
	  index1=ipompc[i];
	  do{
	    index2old=index1;
	    index2=nodempc[3*index1-1];
	    if(index2==0) break;
	    do{
	      if((nodempc[3*index1-3]==nodempc[3*index2-3])&&
		 (nodempc[3*index1-2]==nodempc[3*index2-2])){
		coefmpc[index1-1]+=coefmpc[index2-1];
		nodempc[3*index2old-1]=nodempc[3*index2-1];
		nodempc[3*index2-1]=*mpcfree;
		*mpcfree=index2;
		index2=nodempc[3*index2old-1];
		if(index2==0) break;
	      }
	      else{
		index2old=index2;
		index2=nodempc[3*index2-1];
		if(index2==0) break;
	      }
	    }while(1);
	    index1=nodempc[3*index1-1];
	    if(index1==0) break;
	  }while(1);
		    
	  /* check for zero coefficients on the dependent side */
		    
	  index1=ipompc[i];
	  if(fabs(coefmpc[index1-1])<1.e-10){
	    printf(" *ERROR in cascade: zero coefficient on the\n");
	    printf("        dependent side of an equation\n");
	    printf("        dependent node: %" ITGFORMAT "",nodempc[3*index1-3]);
	    printf("        direction: %" ITGFORMAT "\n\n",nodempc[3*index1-2]);
	    FORTRAN(stop,());
	  }

	  /* a term in MPC i is being replaced; 
	     the coefficient of this term is coef;
	     the term following this term is at position indexold; */
			    
	  indexold=nodempc[3*index-1];
	  coef=coefmpc[index-1];
		    
	  ichange=1;iexpand=1;
	  if((strcmp1(&labmpc[20*i],"                    ")==0)&&
	     (strcmp1(&labmpc[20*(mpc-1)],"CYCLIC")==0))
	    strcpy1(&labmpc[20*i],"SUBCYCLIC",9);
	  indexnew=ipompc[mpc-1];
	  coef=-coef/coefmpc[indexnew-1];
	  indexnew=nodempc[3*indexnew-1];
	  if(indexnew!=0){
	    do{
	      coefmpc[index-1]=coef*coefmpc[indexnew-1];
	      nodempc[3*index-3]=nodempc[3*indexnew-3];
	      nodempc[3*index-2]=nodempc[3*indexnew-2];
	      indexnew=nodempc[3*indexnew-1];
	      if(indexnew!=0){
		nodempc[3*index-1]=*mpcfree;
		index=*mpcfree;
		*mpcfree=nodempc[3**mpcfree-1];
		if(*mpcfree==0){
		  *mpcfree=*memmpc_+1;
		  nodempc[3*index-1]=*mpcfree;
		  *memmpc_=(ITG)(1.1**memmpc_);
		  printf(" *INFO in cascade: reallocating nodempc; new size = %" ITGFORMAT "\n\n",*memmpc_);
		  RENEW(nodempc,ITG,3**memmpc_);
		  RENEW(coefmpc,double,*memmpc_);
		  for(j=*mpcfree;j<*memmpc_;j++){
		    nodempc[3*j-1]=j+1;
		  }
		  nodempc[3**memmpc_-1]=0;
		}
		continue;
	      }
	      else{
		nodempc[3*index-1]=indexold;
		break;
	      }
	    }while(1);
	  }else{
	    coefmpc[index-1]=0.;
	  }
	  break;
	}
	else{
	  index=nodempc[3*index-1];
	  if(index!=0) continue;
	  else break;
	}
      }while(1);
      if(iexpand==0) continue;
	    
      /* one term of the mpc was expanded 
	 collecting terms corresponding to the same DOF */
	    
      index1=ipompc[i];
      do{
	index2old=index1;
	index2=nodempc[3*index1-1];
	if(index2==0) break;
	do{
	  if((nodempc[3*index1-3]==nodempc[3*index2-3])&&
	     (nodempc[3*index1-2]==nodempc[3*index2-2])){
	    coefmpc[index1-1]+=coefmpc[index2-1];
	    nodempc[3*index2old-1]=nodempc[3*index2-1];
	    nodempc[3*index2-1]=*mpcfree;
	    *mpcfree=index2;
	    index2=nodempc[3*index2old-1];
	    if(index2==0) break;
	  }
	  else{
	    index2old=index2;
	    index2=nodempc[3*index2-1];
	    if(index2==0) break;
	  }
	}while(1);
	index1=nodempc[3*index1-1];
	if(index1==0) break;
      }while(1);
	    
      /* check for zero coefficients on the dependent and
	 independent side */
	    
      index1=ipompc[i];
      index1old=0;
      do {
	//		printf("index1=%d\n",index1);
	//		printf("coefmpc=%e\n",coefmpc[index1-1]);

	if(fabs(coefmpc[index1-1])<1.e-10){
	  if(index1old==0){
	    printf(" *ERROR in cascade: zero coefficient on the\n");
	    printf("        dependent side of an equation\n");
	    printf("        dependent node: %" ITGFORMAT "",nodempc[3*index1-3]);
	    printf("        direction: %" ITGFORMAT "\n\n",nodempc[3*index1-2]);
	    FORTRAN(stop,());
	  }
	  else{
	    nodempc[3*index1old-1]=nodempc[3*index1-1];
	    nodempc[3*index1-1]=*mpcfree;
	    *mpcfree=index1;
	    index1=nodempc[3*index1old-1];
	  }
	}
	else{
	  index1old=index1;
	  index1=nodempc[3*index1-1];
	}
	if(index1==0) break;
      }while(1);
    }
    if(ichange==0) break;
  }while(1);

  /*     determining the effective size of nodempc and coefmpc for
	 the reallocation*/

  *mpcend=0;
  *maxlenmpc=0;
  for(i=0;i<*nmpc;i++){
    index=ipompc[i];
    *mpcend=max(*mpcend,index);
    nterm=1;
    while(1){
      index=nodempc[3*index-1];
      if(index==0){
	*maxlenmpc=max(*maxlenmpc,nterm);
	break;
      }
      *mpcend=max(*mpcend,index);
      nterm++;
    }
  }

  SFREE(jmpc);

  *nodempcp=nodempc;
  *coefmpcp=coefmpc;
    
  /* for(i=0;i<*nmpc;i++){
    j=i+1;
    FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
    }*/
    
  return;
}
