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
/**
 * Calculates the transformation matrices \f$ T \f$ and \f$ T^{-1} \f$ for quad-quad or quad-lin method
 * see phd-thesis Sitzmann Chapter 4.1
 * Author: Saskia Sitzmann
 * 
 *  [in] ntie		number of contraints
 *  [in] ipkon		pointer into field kon...
 *  [in] kon 		.. for element i storing the connectivity list of elem. in succ. order 
 *  [in] nk 		number of nodes
 *  [in] lakon		(i) label for element i 
 *  [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 *  [in] itiefac		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i 
 *  [in] tieset           (i) name of tie i 
 *  [in] islavnode	field storing the nodes of the slave surface
 *  [in] islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i  
 *  [out] irowtp		field containing row numbers of aut
 *  [out] jqt	        pointer into field irowt
 *  [out] autp		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$ 
 *  [out] irowtinvp	field containing row numbers of autinv
 *  [out] jqtinv	pointer into field irowtinv
 *  [out] autinvp	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$ 
*/
void buildtquad(ITG *ntie,ITG *ipkon,ITG *kon,ITG *nk,char *lakon,
		ITG *nslavnode,ITG *itiefac,char *tieset,
		ITG *islavnode,ITG *islavsurf,
		ITG **irowtp,ITG *jqt,double **autp,
		ITG **irowtinvp,ITG *jqtinv,double **autinvp){  
  
  ITG i,j,l,nodesf,nodem,istart,icounter,ndim,ifree,ifree2,
    nzstloc,nzstlocinv,*krow=NULL,*kcol=NULL,
    *mast1=NULL,*mast2=NULL,*irowt=NULL,*irowtinv=NULL;
  
  double contribution,*contr=NULL,*aut=NULL,*autinv=NULL;
  
  irowt=*irowtp; aut=*autp;
  irowtinv=*irowtinvp; autinv=*autinvp;
  
  /** built T and T^-1 **/
  
  nzstloc=3*nslavnode[*ntie];
  NNEW(mast1,ITG,nzstloc);
  RENEW(aut,double,nzstloc);
  RENEW(irowt,ITG,nzstloc);
  
  nzstlocinv=3*nslavnode[*ntie];
  NNEW(mast2,ITG,nzstlocinv);
  RENEW(autinv,double,nzstloc);
  RENEW(irowtinv,ITG,nzstloc);
  
  ifree=1;ifree2=1;
  
  NNEW(contr,double,16);
  NNEW(krow,ITG,16);
  NNEW(kcol,ITG,16);
  
  for(i=0;i<*ntie;i++){	    
    if(tieset[i*(81*3)+80]=='C'){				
      for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){

	/* contribution for T */
	
	FORTRAN(create_t,(ipkon,kon,lakon,islavsurf,
			  contr,krow,kcol,&icounter,&l));
	
	for(j=0;j<icounter;j++){
	  contribution=contr[j];
	  nodesf=krow[j];				
	  nodem=kcol[j];				
	  insertas(&irowt,&mast1,&nodesf,&nodem,&ifree,&nzstloc,
		   &contribution,&aut);				
	}

	/* contribution for T^-1 */
	
	FORTRAN(create_tinv,(ipkon,kon,lakon,islavsurf,
			     contr,krow,kcol,&icounter,&l));
	
	for(j=0;j<icounter;j++){
	  contribution=contr[j];
	  nodesf=krow[j];				
	  nodem=kcol[j];				
	  insertas(&irowtinv,&mast2,&nodesf,&nodem,&ifree2,&nzstlocinv,
		   &contribution,&autinv);
	}
	
      }
    }
  }
  SFREE(contr);SFREE(krow);SFREE(kcol);
    
  nzstloc=ifree-1;
  ndim=*nk;
  matrixsort(aut,mast1,irowt,jqt,&nzstloc,&ndim);
  
  /* Getting rid of identical contributions in T
     (a node can belong to several slave faces) */
  
  icounter=0;
  for (i=0;i<*nk;i++){
    if(jqt[i]!=jqt[i+1]){
      irowt[icounter]=irowt[jqt[i]-1];
      aut[icounter]=aut[jqt[i]-1];
      icounter++;
      istart=icounter;
      for (j=jqt[i];j<jqt[i+1]-1;j++){
	if (irowt[j]==irowt[icounter-1]){
	  aut[icounter-1]=aut[j];   
	}else{
	  irowt[icounter]=irowt[j];
	  aut[icounter]=aut[j];
	  icounter++;
	  }
      }
    }else{ istart=icounter+1;}   
    jqt[i]=istart;
  }
  jqt[*nk]=icounter+1; 
  RENEW(irowt,ITG,icounter+1);
  RENEW(aut,double,icounter+1); 
  SFREE(mast1);	
  
  nzstlocinv=ifree2-1; 
  ndim=*nk;
  matrixsort(autinv,mast2,irowtinv,jqtinv,&nzstlocinv,&ndim);  
 
  /* Getting rid of identical contributions in T^-1
     (a node can belong to several slave faces) */
  
  icounter=0;
  for (i=0;i<*nk;i++){
    if(jqtinv[i]!=jqtinv[i+1]){
      irowtinv[icounter]=irowtinv[jqtinv[i]-1];
      autinv[icounter]=autinv[jqtinv[i]-1];
      icounter++;
      istart=icounter;
      for (j=jqtinv[i];j<jqtinv[i+1]-1;j++){
	if (irowtinv[j]==irowtinv[icounter-1]){
	  autinv[icounter-1]=autinv[j];   
	}else{
	  irowtinv[icounter]=irowtinv[j];
	  autinv[icounter]=autinv[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}   
    jqtinv[i]=istart;
  }
  jqtinv[*nk]=icounter+1;
  
  RENEW(irowtinv,ITG,icounter+1);
  RENEW(autinv,double,icounter+1); 
  SFREE(mast2);
  
  *irowtp=irowt;*autp=aut;
  *irowtinvp=irowtinv;*autinvp=autinv;
  
  return;
}	    

  
