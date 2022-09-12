/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

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
#include <time.h>
#include "CalculiX.h"
#include "mortar.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

/**  modify n and t  due to SPC's/MPC's
 * needed for 2D calculations
 * author: Saskia Sitzmann
 *  [in,out] n		slave normal
 *  [in,out] t		slave tangent
 *  [out] n2		transformed slave normal
 *  [out] that		transformed slave tangent
 *  [in] islavnodeentry	current slave node entry
 *  [in] nslavspc		(2*i) pointer to islavspc...
 *  [in] islavspc         ... which stores SPCs for slave node i
 *  [in] nsspc            number of SPC for slave nodes
 *  [in] nslavmpc		(2*i) pointer to islavmpc...
 *  [in] islavmpc		... which stores MPCs for slave node i
 *  [in] nsmpc		number of MPC for slave nodes 
 *  [in] nmastspc		(2*i) pointer to imastspc...
 *  [in] imastspc         ... which stores SPCs for master node i
 *  [in] nmspc            number of SPC for master nodes
 *  [in] nmastmpc		(2*i) pointer to imastmpc...
 *  [in] imastmpc		... which stores MPCs for master node i
 *  [in] nmmpc		number of MPC for master nodes
 *  [in] debug 		debug output flag 
 *  [in] node 	 	current slave node
 */

void trafontspcmpc( double *n,double *t,double *n2,double *that,
		    ITG *islavnodeentry,ITG *nboun,ITG *ndirboun,ITG *nodeboun,
		    double *xboun,ITG *nmpc,ITG *ipompc,ITG *nodempc,
		    double *coefmpc,ITG *ikboun,ITG *ilboun,ITG *ikmpc,
		    ITG *ilmpc,ITG *nslavspc,ITG *islavspc,ITG *nsspc,
		    ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,ITG *nmastspc,
		    ITG *imastspc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,
		    ITG *nmmpc,ITG *debug,ITG *node){
  
  ITG jj,kk,ist,dir,node1,index,dirind,dirdep,nt1;
  
  double dist,nn,coefdep;
  
  /** handle SPCs and MPCs -> modify tangentials and normal **/
  
  for(jj=nslavspc[2*(*islavnodeentry-1)];jj<nslavspc[2*(*islavnodeentry-1)+1];
      jj++){	    
    t[0]=0;t[1]=0;t[2]=0;	    
    ist=islavspc[2*jj];            
    dir=ndirboun[ist-1];	    
    node1=nodeboun[ist-1];	    
    if(*debug==1){printf("jj %" ITGFORMAT " ist %" ITGFORMAT " dir %"
			 ITGFORMAT " node %" ITGFORMAT "\n",jj,ist,dir,node1);}
    t[dir-1]=1.0;	    
    dist=t[0]*n[0]+t[1]*n[1]+t[2]*n[2];	    
    n[0]=n[0]-dist*t[0];	    
    n[1]=n[1]-dist*t[1];	    
    n[2]=n[2]-dist*t[2];	    
    nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);	    
    n[0]=n[0]/nn;n[1]=n[1]/nn;n[2]=n[2]/nn;	   
  }
  for(jj=nslavmpc[2*(*islavnodeentry-1)];jj<nslavmpc[2*(*islavnodeentry-1)+1];
      jj++){	    	   
    t[0]=0;t[1]=0;t[2]=0;             
    ist=islavmpc[2*jj];            
    dirdep=nodempc[3*(ist-1)+1];            
    coefdep=coefmpc[ist-1];            
    index=nodempc[3*(ist-1)+2];	    
    node1=nodempc[3*(ist-1)]; 	    
    t[dirdep-1]=coefdep;	    
    if(*debug==1){printf("jj %" ITGFORMAT " ist %" ITGFORMAT " dir %"
			 ITGFORMAT " node %" ITGFORMAT " coef %e\n",jj,ist,
			 dirdep,node1,coefdep);}	                
    while(index!=0){
      dirind=nodempc[3*(index-1)+1];
      node1=nodempc[3*(index-1)]; 
      if(node1==*node)t[dirind-1]=coefmpc[index-1];
      index=nodempc[3*(index-1)+2];
      if(*debug==1){printf("index %" ITGFORMAT " \n",index);}
    }
    nt1=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
    t[0]=t[0]/nt1;t[1]=t[1]/nt1;t[2]=t[2]/nt1;
    dist=t[0]*n[0]+t[1]*n[1]+t[2]*n[2];
    n[0]=n[0]-dist*t[0];
    n[1]=n[1]-dist*t[1];
    n[2]=n[2]-dist*t[2];
    nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    n[0]=n[0]/nn;n[1]=n[1]/nn;n[2]=n[2]/nn;	      
  }
  
  /** get second tangential and modify it due to SPC's/MPC's**/
  
  if(nslavspc[2*(*islavnodeentry-1)+1]-nslavspc[2*(*islavnodeentry-1)]>0 ||
     nslavmpc[2*(*islavnodeentry-1)+1]-nslavmpc[2*(*islavnodeentry-1)]>0){
    
    /** caution: t[3-5]= hat(t2)=n x hat(t1) **/
    
    t[3]=n[1]*t[2]-n[2]*t[1];            
    t[4]=n[2]*t[0]-n[0]*t[2];            
    t[5]=n[0]*t[1]-n[1]*t[0];	    
    for(kk=0;kk<6;kk++){that[kk]=t[kk];}	    
    for(jj=nslavmpc[2*(*islavnodeentry-1)];
	jj<nslavmpc[2*(*islavnodeentry-1)+1];jj++){             
      ist=islavmpc[2*jj];             
      dirdep=nodempc[3*(ist-1)+1];             
      coefdep=coefmpc[ist-1];             
      index=nodempc[3*(ist-1)+2];	     
      node1=nodempc[3*(ist-1)];                
      while(index!=0){               
	dirind=nodempc[3*(index-1)+1];	       
	node1=nodempc[3*(index-1)]; 	       
	if(node1==*node){		 
	  t[dirind-1+3]=t[dirind-1+3]-coefmpc[index-1]*t[dirdep-1+3]/coefdep;
	  n[dirind-1]=n[dirind-1]-coefmpc[index-1]*n[dirdep-1]/coefdep;
	}
	index=nodempc[3*(index-1)+2];	    
      }
      t[dirdep-1+3]=0.0;	    
      n[dirdep-1]=0.0;	    	   
    }	     
    nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);	
    nt1=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
    if(*debug==1){
      printf("n= %e %e %e \n",n[0],n[1],n[2]);
      printf("t1= %e %e %e \n",t[0],t[1],t[2]);
      printf("t2= %e %e %e \n",t[3],t[4],t[5]);
    }
    
    /** caution: t[3-5]= tilde(t2),hat(t2) with modifications due to SPC/MPCs
     *        that[0-5] =hat(t1) ,hat(t2)
     *	n[0-2]= n*hat(D)/D in case of directinal blocking
     *       n[0-2]= n otherwise 
     *       n2[0-2]= original n
     **/
    
  }
  if(*debug==1){
    printf("n= %e %e %e \n",n[0],n[1],n[2]);
    printf("t1= %e %e %e \n",t[0],t[1],t[2]);
    printf("t2= %e %e %e \n",t[3],t[4],t[5]);
  }
  return;
      
}
