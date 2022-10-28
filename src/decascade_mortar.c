/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

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
#include <time.h>
#include "CalculiX.h"
#include "mortar.h"
/**  function to decascade generated transformed MPCs
 * Author: Saskia Sitzmann
 *
 *  [in,out] nmpc		number of mpcs
 *  [in,out] ipompc       (i) pointer to nodempc and coeffmpc for MPC i
 *  [in,out] nodempcp      nodes and directions of MPCs
 *  [in,out] coefmpcp      coefficients of MPCs
 *  [in,out] ikmpc 	sorted dofs idof=8*(node-1)+dir for MPCs
 *  [in,out] ilmpc	SPC numbers for sorted dofs
 *  [in,out] memmpc_	size of nodempc/coefmpc
 *  [in,out] mpcfree	marking the next free space in nodempc
 **/
void decascade_mortar(ITG *nmpc,ITG *ipompc,ITG **nodempcp,double **coefmpcp,
		      ITG *ikmpc,ITG *ilmpc, ITG *memmpc_, ITG *mpcfree){
  
  ITG i,j,k,l,nl,index,indexold,index1,index2,index2old,indexnew,index1old,
    iexpand,ichange,idof,id,mpc,*nodempc=NULL;
  
  double coef,*coefmpc=NULL;
  
  nodempc=*nodempcp; coefmpc=*coefmpcp;
  /*     decascading */
  
  
  /* decascading using simple substitution */
  
  for(i=0;i<*nmpc;i++){
      //if(i+1==1060||i+1==1061)printf("mpc %"ITGFORMAT"\n",i+1);
      iexpand=0;
      index=nodempc[3*ipompc[i]-1];
      if(index==0) continue;
      do{
	//if(i+1==1060||i+1==1061)printf("\tindex %"ITGFORMAT" node %"ITGFORMAT" dir %"ITGFORMAT" \n",index,nodempc[3*index-3],nodempc[3*index-2]);
	idof=(nodempc[3*index-3]-1)*8+nodempc[3*index-2];
	FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
	indexold=nodempc[3*index-1];
	coef=coefmpc[index-1];
	if((id>0)&&(ikmpc[id-1]==idof)){
	  /* a term on the independent side of the MPC is
	     detected as dependent node in another MPC */ 
	  iexpand=1;
	  mpc=ilmpc[id-1];
	  
	  indexnew=ipompc[mpc-1];
	  coef=-coef/coefmpc[indexnew-1];
	  indexnew=nodempc[3*indexnew-1];
	  do{
	    coefmpc[index-1]=coef*coefmpc[indexnew-1];
	    nodempc[3*index-3]=nodempc[3*indexnew-3];
	    nodempc[3*index-2]=nodempc[3*indexnew-2];
	    
	    indexnew=nodempc[3*indexnew-1];
	    if(indexnew!=0){
	      if(*mpcfree==0){
		*mpcfree=*memmpc_+1;
		*memmpc_=(ITG)(1.1**memmpc_);
		//printf("*INFO in cascade: reallocating nodempc; new size = %" ITGFORMAT "\n\n",*memmpc_);
		RENEW(nodempc,ITG,3**memmpc_);
		RENEW(coefmpc,double,*memmpc_);
		for(j=*mpcfree;j<*memmpc_;j++){
		  nodempc[3*j-1]=j+1;
		}
		nodempc[3**memmpc_-1]=0;
	      }	      
	      nodempc[3*index-1]=*mpcfree;
	      //if(i+1==1060||i+1==1061)printf("\t\tindex %"ITGFORMAT" node %"ITGFORMAT" dir %"ITGFORMAT" nindex %"ITGFORMAT" \n",index,nodempc[3*index-3],nodempc[3*index-2],nodempc[3*index-1]); 
	      index=*mpcfree;
	      *mpcfree=nodempc[3**mpcfree-1];        
	      continue;
	    }
	    else{
	      nodempc[3*index-1]=indexold;
	      //if(i+1==1060||i+1==1061)printf("\t\tindex %"ITGFORMAT" node %"ITGFORMAT" dir %"ITGFORMAT" nindex %"ITGFORMAT" \n",index,nodempc[3*index-3],nodempc[3*index-2],nodempc[3*index-1]); 
	      break;
	    }
	  }while(1);
	}
	index=indexold;
	if(index!=0) continue;
	  else break;
      }while(1);
      //printf("\n");
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
      
  }

  
  *nodempcp=nodempc; *coefmpcp=coefmpc;
  return;
}
