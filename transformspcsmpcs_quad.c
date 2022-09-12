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
 *  transform SPCs/MPCs for quadratic slave elements needed for dual mortar contact
 * Author: Saskia Sitzmann
 * @todo: debug error with mid edge nodes
 *
 *  [out] nboun2           number of transformed SPCs
 *  [out] ndirboun2p	(i) direction of transformed SPC i 
 *  [out] nodeboun2p      (i) node of transformed SPC i
 *  [out] xboun2p         (i) value of transformed SPC i
 *  [out] nmpc2		number of transformed mpcs
 *  [out] ipompc2p        (i) pointer to nodempc and coeffmpc for transformed MPC i
 *  [out] nodempc2p       nodes and directions of transformed MPCs
 *  [out] coefmpc2p       coefficients of transformed MPCs
 *  [out] labmpc2p	transformed MPCs labels 
 *  [out] ikboun2p        sorted dofs idof=8*(node-1)+dir for transformed SPCs
 *  [out] ilboun2p        transformed SPC numbers for sorted dofs
 *  [out] ikmpc2p 	sorted dofs idof=8*(node-1)+dir for transformed MPCs
 *  [out] ilmpc2p		transformed SPC numbers for sorted dofs 
 *  [in] irowtlocinv	field containing row numbers of autlocinv
 *  [in] jqtlocinv	pointer into field irowtlocinv
 *  [in] autlocinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$  
 *  [out] nk2		number or generated points needed for transformed SPCs 
 *  [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, 
 =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad) 
 *  [out] nodeforc2p	transformed point force, node
 *  [out] ndirforc2p	transformed point force, dir
 *  [out] xforc2p		transformed point force, value
 *  [out] nforc2		number of transformed point forces  
 **/

void transformspcsmpcs_quad(ITG *nboun,ITG *ndirboun,ITG *nodeboun,
			    double *xboun,
			    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
			    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,ITG *nboun2,
			    ITG **ndirboun2p,ITG **nodeboun2p,double **xboun2p,
			    ITG *nmpc2,ITG **ipompc2p,ITG **nodempc2p,double **coefmpc2p,
			    char **labmpc2p,ITG **ikboun2p,ITG **ilboun2p,ITG **ikmpc2p,
			    ITG **ilmpc2p,ITG *irowtlocinv, ITG *jqtlocinv,double *autlocinv, 
			    ITG *nk,ITG *nk2,ITG *iflagdualquad,
			    ITG *ntie, char *tieset, ITG *itiefac,ITG *islavsurf,
			    char *lakon,ITG *ipkon,ITG *kon,ITG *mt,ITG *memmpc_,
			    ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
			    ITG **nodeforc2p,ITG **ndirforc2p,double **xforc2p,ITG *nforc2){
  
  char *labmpc2=NULL;

  ITG i,j,jj,ndimboun2,ndimmpc2a,ndimmpc2b,node,nodedep,nodeind,
    ndir,ndirdep,ndirind,ist,index,ifree,idof,id,nhelp,debug,
    mpcfree2,memmpc_2,nodevdep,nodewdep,nodevind,nodewind,newnode,
    *ndirboun2=NULL, *nodeboun2=NULL, *ipompc2=NULL, *nodempc2=NULL,
    ndimforc2,*ikboun2=NULL, *ilboun2=NULL, *ikmpc2=NULL, *ilmpc2=NULL,
    *nodeforc2=NULL,*ndirforc2=NULL;
  
  double alpha,a1,b1,fixed_disp,coeffdep,coeffind,fixed_forc,
    *xboun2=NULL, *coefmpc2=NULL, *xforc2=NULL; 

  debug=0;
  
  ndirboun2=*ndirboun2p; nodeboun2=*nodeboun2p; xboun2=*xboun2p;
  ipompc2=*ipompc2p; nodempc2=*nodempc2p; coefmpc2=*coefmpc2p;
  ikboun2=*ikboun2p; ilboun2=*ilboun2p; ikmpc2=*ikmpc2p; ilmpc2=*ilmpc2p;  
  labmpc2=*labmpc2p;nodeforc2=*nodeforc2p; ndirforc2=*ndirforc2p;
  xforc2=*xforc2p;
  
  // initialize
  
  *nk2=0;
  *nboun2=0;
  ndimboun2=*nboun;
  RENEW(ndirboun2,ITG, *nboun);
  RENEW(nodeboun2,ITG, *nboun);
  RENEW(xboun2,double,*nboun);
  *nmpc2=0;
  ndimmpc2a=*nmpc;
  ndimmpc2b=3**memmpc_;
  RENEW(ipompc2,ITG, ndimmpc2a);
  RENEW(nodempc2,ITG, 3*ndimmpc2b);
  RENEW(coefmpc2,double,ndimmpc2b);
  RENEW(ikboun2,ITG, *nboun);
  RENEW(ilboun2,ITG, *nboun);
  RENEW(ikmpc2,ITG, *nmpc);
  RENEW(ilmpc2,ITG, *nmpc);
  RENEW(labmpc2,char,20**nmpc+1);
  *nforc2=0; 
  ndimforc2=*nforc;
  RENEW(nodeforc2,ITG,2**nforc);
  RENEW(ndirforc2,ITG,*nforc);
  RENEW(xforc2,double,*nforc);  
  
  ifree=1;
  a1=0.0;
  b1=0.0;
  if(*iflagdualquad==1 || *iflagdualquad==3){
    alpha=1.0/2.0;
    b1=1.0;
    a1=alpha;
  }else{
    alpha=1.0/5.0;
    b1=(1-2*alpha);
    a1=alpha; 
  }
  
  if(debug==1)printf(" transformspcsmpcs_quad: nboun %" ITGFORMAT "\t",*nboun);
  
  // loop over all SPCs
  
  for(i=0;i<*nboun;i++){
    fixed_disp=xboun[i];
    ndir=ndirboun[i];
    node=nodeboun[i];
    idof=8*(node-1)+ndir;
    nhelp=0;
    nodevdep=0;
    nodewdep=0;
    
    if(nhelp==3 && ndirboun[i]<4){

      /* slave node & mid node for quadratic elements
	 spc is transformed to mpc
	 create new node */
	
      *nk2=*nk2+1;
      newnode=*nk+*nk2;
      
      // add spc's
      
      for(jj=0;jj<4;jj++){
	if(jj==ndir){
	  xboun2[*nboun2]=fixed_disp;
	}else{
	  xboun2[*nboun2]=0.0;   
	}
	ndirboun2[*nboun2]=jj;
	nodeboun2[*nboun2]=newnode;
	idof=8*(newnode-1)+jj;
	FORTRAN(nident,(ikboun2,&idof,nboun2,&id));
	
	for( j=*nboun2;j>id;j--){
	  ikboun2[j]=ikboun2[j-1];
	  ilboun2[j]=ilboun2[j-1];
	}
	ikboun2[id]=idof;
	ilboun2[id]=*nboun2+1;    
	*nboun2=*nboun2+1;
	if(*nboun2+1>ndimboun2){
	  ndimboun2=ndimboun2*1.5+1;
	  RENEW(xboun2,double,ndimboun2);
	  RENEW(ndirboun2,ITG, ndimboun2);
	  RENEW(nodeboun2,ITG, ndimboun2);
	  RENEW(ikboun2,ITG, ndimboun2);
	  RENEW(ilboun2,ITG, ndimboun2);
	}
      }
      
      // add mpc
      
      if(*nmpc2+1>ndimmpc2a){
	ndimmpc2a=ndimmpc2a*1.5+1;
	RENEW(ipompc2,ITG, ndimmpc2a);
	RENEW(ikmpc2,ITG, ndimmpc2a);
	RENEW(ilmpc2,ITG, ndimmpc2a);
	RENEW(labmpc2,char,20*ndimmpc2a+1);
      } 
      if(ifree+4>ndimmpc2b){
	ndimmpc2b=ndimmpc2b*1.5+4;
	RENEW(coefmpc2,double,ndimmpc2b);
	RENEW(nodempc2,ITG, 3*ndimmpc2b);
      }	
      for(jj=0;jj<20;jj++){
	labmpc2[20*(*nmpc2)+jj]=' ';}
      ipompc2[*nmpc2]=ifree;
      
      nodempc2[0+(ifree-1)*3]=node;
      nodempc2[1+(ifree-1)*3]=ndir;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=(b1);
      
      ifree++;
      nodempc2[0+(ifree-1)*3]=nodevdep;
      nodempc2[1+(ifree-1)*3]=ndir;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=(a1);
      
      ifree++;
      nodempc2[0+(ifree-1)*3]=nodewdep;
      nodempc2[1+(ifree-1)*3]=ndir;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=(a1);
      
      ifree++;
      nodempc2[0+(ifree-1)*3]=newnode;
      nodempc2[1+(ifree-1)*3]=ndir;
      nodempc2[2+(ifree-1)*3]=0;
      coefmpc2[ifree-1]=-1.0;
      
      ifree++;
      idof=8*(node-1)+ndir;
      FORTRAN(nident,(ikmpc2,&idof,nmpc2,&id));
      for( j=*nmpc2;j>id;j--){
	ikmpc2[j]=ikmpc2[j-1];
	ilmpc2[j]=ilmpc2[j-1];
      }
      ikmpc2[id]=idof;
      ilmpc2[id]=*nmpc2+1;    
      *nmpc2=*nmpc2+1; 
    }else{

      /* no slave node or corner slave node
	 nothing to do*/
	
      xboun2[*nboun2]=fixed_disp;
      ndirboun2[*nboun2]=ndir;
      nodeboun2[*nboun2]=node;
      FORTRAN(nident,(ikboun2,&idof,nboun2,&id));
      for( j=*nboun2;j>id;j--){
	ikboun2[j]=ikboun2[j-1];
	ilboun2[j]=ilboun2[j-1];
      }
      ikboun2[id]=idof;
      ilboun2[id]=*nboun2+1;    
      *nboun2=*nboun2+1;
      if(*nboun2+1>ndimboun2){
	ndimboun2=ndimboun2*1.5+1;
	RENEW(xboun2,double,ndimboun2);
	RENEW(ndirboun2,ITG, ndimboun2);
	RENEW(nodeboun2,ITG, ndimboun2);
	RENEW(ikboun2,ITG, ndimboun2);
	RENEW(ilboun2,ITG, ndimboun2);
      }
    }  
  }
  if(debug==1)printf("nboun2 %" ITGFORMAT "\n",*nboun2);
  RENEW(xboun2,double,*nboun2);
  RENEW(ndirboun2,ITG, *nboun2);
  RENEW(nodeboun2,ITG, *nboun2);
  
  // loop over all MPCs
  
  if(debug==1)printf(" transformspcsmpcs_quad: nmpc %" ITGFORMAT "\t",*nmpc);
  for(i=0;i<*nmpc;i++){ 
    ist=ipompc[i];
    nodedep=nodempc[0+(ist-1)*3];
    ndirdep=nodempc[1+(ist-1)*3];
    coeffdep=coefmpc[ist-1];
    nhelp=0;

    if(*nmpc2+1>ndimmpc2a){
      ndimmpc2a=ndimmpc2a*1.5+1;
      RENEW(ipompc2,ITG, ndimmpc2a);
      RENEW(ikmpc2,ITG, ndimmpc2a);
      RENEW(ilmpc2,ITG, ndimmpc2a);
      RENEW(labmpc2,char,20*ndimmpc2a+1);
    }    
    for(jj=0;jj<20;jj++){
      labmpc2[20*(*nmpc2)+jj]=labmpc[20*i+jj];}
    if(nhelp==3 && ndirdep<4){// slave node & mid node for quadratic elements

      // find adjacent nodes
	
      if(ifree+3>ndimmpc2b){
	ndimmpc2b=ndimmpc2b*1.5+3;
	RENEW(coefmpc2,double,ndimmpc2b);
	RENEW(nodempc2,ITG, 3*ndimmpc2b);
      }	
      ipompc2[*nmpc2]=ifree;
      nodempc2[0+(ifree-1)*3]=nodedep;
      nodempc2[1+(ifree-1)*3]=ndirdep;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=coeffdep*(b1);
      
      ifree++;
      nodempc2[0+(ifree-1)*3]=nodevdep;
      nodempc2[1+(ifree-1)*3]=ndirdep;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=coeffdep*(a1);
      
      ifree++;
      nodempc2[0+(ifree-1)*3]=nodewdep;
      nodempc2[1+(ifree-1)*3]=ndirdep;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=coeffdep*(a1);

      ifree++;
    }else{ // no slave node or corner slave node or mpc not for displacements 
      ipompc2[*nmpc2]=ifree;
      nodempc2[0+(ifree-1)*3]=nodedep;
      nodempc2[1+(ifree-1)*3]=ndirdep;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=coeffdep;

      ifree++;
      if(ifree>ndimmpc2b){
	ndimmpc2b=ndimmpc2b*1.5+1;
	RENEW(coefmpc2,double,ndimmpc2b);
	RENEW(nodempc2,ITG, 3*ndimmpc2b);
      }
    }
    index=nodempc[2+(ist-1)*3];
    if(index!=0){
      do{
	nodeind=nodempc[0+(index-1)*3];
	ndirind=nodempc[1+(index-1)*3];
	coeffind=coefmpc[index-1];
	nhelp=0;
	
	if(nhelp==3 && ndirdep<4){// slave node & mid node for quadratic elements
	  // find adjacent nodes
	    
	  if(ifree+3>ndimmpc2b){
	    ndimmpc2b=ndimmpc2b*1.5+3;
	    RENEW(coefmpc2,double,ndimmpc2b);
	    RENEW(nodempc2,ITG, 3*ndimmpc2b);
	  }
	  nodempc2[0+(ifree-1)*3]=nodeind;
	  nodempc2[1+(ifree-1)*3]=ndirind;
	  nodempc2[2+(ifree-1)*3]=ifree+1;   
	  coefmpc2[ifree-1]=coeffind*(b1);
	  ifree++;

	  nodempc2[0+(ifree-1)*3]=nodevind;
	  nodempc2[1+(ifree-1)*3]=ndirind;
	  nodempc2[2+(ifree-1)*3]=ifree+1;
	  coefmpc2[ifree-1]=coeffind*(a1);
	  ifree++;

	  nodempc2[0+(ifree-1)*3]=nodewind;
	  nodempc2[1+(ifree-1)*3]=ndirind;
	  coefmpc2[ifree-1]=coeffind*(a1);
	  ifree++;	    	    
	}else{// no slave node or corner slave node or mpc not for displacements   
	  nodempc2[0+(ifree-1)*3]=nodeind;
	  nodempc2[1+(ifree-1)*3]=ndirind;
	  coefmpc2[ifree-1]=coeffind;
	  ifree++;
	  if(ifree>ndimmpc2b){
	    ndimmpc2b=ndimmpc2b*1.5+1;
	    RENEW(coefmpc2,double,ndimmpc2b);
	    RENEW(nodempc2,ITG, 3*ndimmpc2b);
	  }	    
	  
	}	
	index=nodempc[2+(index-1)*3];
	if(index==0){
	  nodempc2[2+(ifree-2)*3]=0;
	  break;
	}else{
	  nodempc2[2+(ifree-2)*3]=ifree;
	}
      }while(1);	
    }
    idof=8*(nodedep-1)+ndirdep;
    FORTRAN(nident,(ikmpc2,&idof,nmpc2,&id));
    for( j=*nmpc2;j>id;j--){
      ikmpc2[j]=ikmpc2[j-1];
      ilmpc2[j]=ilmpc2[j-1];
    }
    ikmpc2[id]=idof;
    ilmpc2[id]=*nmpc2+1;    
    *nmpc2=*nmpc2+1;
    if(*nmpc2+1>ndimmpc2a){
      ndimmpc2a=ndimmpc2a*1.5+1;
      RENEW(ipompc2,ITG, ndimmpc2a);
      RENEW(ikmpc2,ITG, ndimmpc2a);
      RENEW(ilmpc2,ITG, ndimmpc2a);
      RENEW(labmpc2,char,20*ndimmpc2a+1);
    } 
  }
  RENEW(ipompc2,ITG, *nmpc2);
  RENEW(ikmpc2,ITG, *nmpc2);
  RENEW(ilmpc2,ITG, *nmpc2);
  RENEW(labmpc2,char,20**nmpc2+1);
  RENEW(coefmpc2,double,(ifree-1));
  RENEW(nodempc2,ITG, 3*(ifree-1));
  if(debug==1)printf("nmpc2 %" ITGFORMAT "\n",*nmpc2);
  
  if(debug==1)printf(" transformspcsmpcs_quad: nforc %" ITGFORMAT "\t",*nforc);
  
  // loop over all point forces
  
  for(i=0;i<*nforc;i++){
    fixed_forc=xforc[i];
    ndir=ndirforc[i];
    node=nodeforc[2*i];
    idof=8*(node-1)+ndir;
    nhelp=0;
    nodevdep=0;
    nodewdep=0;

    if(nhelp==3){
      if(*nforc2+4>ndimforc2){
	ndimforc2=ndimforc2*1.5+3;
	RENEW(xforc2,double,ndimforc2);
	RENEW(ndirforc2,ITG, ndimforc2);
	RENEW(nodeforc2,ITG, 2*ndimforc2);
      }
      
      //mid node
      
      xforc2[*nforc2]=fixed_forc*b1;
      ndirforc2[*nforc2]=ndir;
      nodeforc2[2**nforc2]=node;
      *nforc2=*nforc2+1;
      
      //corner node 1
      
      xforc2[*nforc2]=fixed_forc*a1;
      ndirforc2[*nforc2]=ndir;
      nodeforc2[2**nforc2]=nodevdep;
      *nforc2=*nforc2+1;
      
      //corner node 2
      
      xforc2[*nforc2]=fixed_forc*a1;
      ndirforc2[*nforc2]=ndir;
      nodeforc2[2**nforc2]=nodewdep;
      *nforc2=*nforc2+1;      
      
    }else{
	
      //nothing to do
	
      xforc2[*nforc2]=fixed_forc;
      ndirforc2[*nforc2]=ndir;
      nodeforc2[2**nforc2]=node;
      *nforc2=*nforc2+1;
      if(*nforc2+1>ndimforc2){
	ndimforc2=ndimforc2*1.5+1;
	RENEW(xforc2,double,ndimforc2);
	RENEW(ndirforc2,ITG, ndimforc2);
	RENEW(nodeforc2,ITG, 2*ndimforc2);
      }      
    }
  }
  if(debug==1)printf("nforc2 %" ITGFORMAT "\n\n",*nforc2); 
  
  //decascade mpcs
   
  mpcfree2=0;
  memmpc_2=ifree-1;
  decascade_mortar(nmpc2,ipompc2,&nodempc2,&coefmpc2,
		   ikmpc2,ilmpc2,&memmpc_2,&mpcfree2);

  *ndirboun2p=ndirboun2; *nodeboun2p=nodeboun2; *xboun2p=xboun2;
  *ipompc2p=ipompc2; *nodempc2p=nodempc2; *coefmpc2p=coefmpc2;
  *ikboun2p=ikboun2; *ilboun2p=ilboun2; *ikmpc2p=ikmpc2; *ilmpc2p=ilmpc2;
  *labmpc2p=labmpc2;
  *nodeforc2p=nodeforc2;*ndirforc2p=ndirforc2;
  *xforc2p=xforc2;
  
  return;
}
