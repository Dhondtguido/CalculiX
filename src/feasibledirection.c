/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                     */

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

#ifdef ARPACK

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "CalculiX.h"
#ifdef SPOOLES
#include "spooles.h"
#endif
#ifdef SGI
#include "sgi.h"
#endif
#ifdef TAUCS
#include "tau.h"
#endif
#ifdef MATRIXSTORAGE
#include "matrixstorage.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif

void feasibledirection(ITG *nobject,char **objectsetp,double **dgdxglobp,
		       double *g0,ITG *ndesi,ITG *nodedesi,ITG *nk,ITG *isolver,
		       ITG **ipkonp,ITG **konp,char **lakonp,ITG *ne,
		       ITG *nelemload,ITG *nload,ITG *nodeboun,ITG *nboun,
		       ITG *ndirboun,ITG *ithermal,double *co,double *vold,
		       ITG *mi,ITG **ielmatp,ITG *ielprop,double *prop,
		       ITG *kode,
		       ITG *nmethod,char *filab,ITG *nstate_,ITG *istep,
		       double *cs,char *set,ITG *nset,ITG *istartset,
		       ITG *iendset,
		       ITG *ialset,char *jobnamec,char *output,ITG *ntrans,
		       ITG *inotr,double *trab,char *orname,double *xdesi,
		       double *timepar,double *coini,ITG *ikboun,ITG *nactdof,
		       ITG *ne2d){
               
  /* finding a feasible direction based on the sensitivity information */
  
  char *objectset=NULL,cflag[1]=" ",description[13]="            ",*lakon=NULL,
    *gradprojname=NULL,*lakonfa=NULL;
  
  ITG i,*ipoacti=NULL,*iconstacti=NULL,nactive=0,nnlconst,iscaleflag,mode=-1,
    *inameacti=NULL,nconstraint=0,*inum=NULL,iinc=1,noddiam=-1,ngraph=1,
    idesvar=0,
    inorm=0,irand=0,ishape=0,icoordinate=1,*ipkon=NULL,iobject,istart,*kon=NULL,
    *ielmat=NULL,*nodedesiinv=NULL,ifeasd=0,methodfd,igeoconst=0,
    *nodedesipos=NULL,
    nkini,*ipoface=NULL,*nodface=NULL,*konfa=NULL,*ipkonfa=NULL,*iponoelfa=NULL,
    *inoelfa=NULL,nsurfs,ifreemax,inoelsize,iregion=0,*iponod2dto3d=NULL,inode,
    *nodedesiboun=NULL,addout,ipos;
         
  double *objnorm=NULL,*dgdxglob=NULL,*stn=NULL,ptime=0.,*gradproj=NULL,
    *extnorini=NULL,*extnor=NULL,distmin,*feasdir=NULL;
       
  objectset=*objectsetp;dgdxglob=*dgdxglobp;ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;ielmat=*ielmatp;

  /* if addout=1 additional output about the feasible direction 
     is written in the frd file */
  addout=0;
  
  methodfd=timepar[3];
  NNEW(feasdir,double,2**nk);
  NNEW(gradproj,double,3**nk);
  NNEW(gradprojname,char,405);
  NNEW(nodedesiinv,ITG,*nk);
  for(i=0;i<*ndesi;i++){
    nodedesiinv[nodedesi[i]-1]=1;
  }
  
  /* createinum is called in order to determine the nodes belonging
     to elements; this information is needed in frd_se */

  NNEW(inum,ITG,*nk);
  FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,nload,
                      nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
                      ielprop,prop));

  /***************************************************************************/
  /* assessment of geometrical constraints                                    /
  /***************************************************************************/
  
  for(i=0;i<*nobject;i++){
    if(strcmp1(&objectset[i*405+404],"G")==0){
      igeoconst=1;		          
    }
  }
  
  if((igeoconst==1)||(strcmp1(&objectset[86],"BOU")==0)){
    
    /*finding the normal directions */

    NNEW(nodedesipos,ITG,*nk);
    NNEW(ipoface,ITG,*nk);
    NNEW(nodface,ITG,5*6**ne);
    NNEW(konfa,ITG,8*6**ne);
    NNEW(ipkonfa,ITG,6**ne);
    NNEW(lakonfa,char,8*6**ne);
    
    FORTRAN(findextsurface,(nodface,ipoface,ne,ipkon,lakon,kon,
			    konfa,ipkonfa,nk,lakonfa,&nsurfs,
			    &ifreemax));
    
    RENEW(nodface,ITG,5*ifreemax);
    RENEW(konfa,ITG,8*nsurfs);
    RENEW(ipkonfa,ITG,nsurfs);
    RENEW(lakonfa,char,8*nsurfs);

    NNEW(iponoelfa,ITG,*nk);
    NNEW(inoelfa,ITG,3*8*6**ne);
   
    FORTRAN(extfacepernode,(iponoelfa,inoelfa,lakonfa,ipkonfa,konfa,
			    &nsurfs,&inoelsize));
   
    RENEW(inoelfa,ITG,3*inoelsize);
    NNEW(extnorini,double,3**nk);
    NNEW(extnor,double,3**nk);
    
    /* calculating the normals w.r.t. to the initial geometry */
    
    FORTRAN(normalsonsurface_se,(ipkon,kon,lakon,extnorini,coini,nk,ipoface,
			         nodface,nactdof,mi,nodedesiinv,&iregion,
			         iponoelfa,ndesi,nodedesi,iponod2dto3d,
			         ikboun,nboun,ne2d)); 		          

    /* calculating the normals w.r.t. the actual geoemtry */
    
    FORTRAN(normalsonsurface_se,(ipkon,kon,lakon,extnor,co,nk,ipoface,
			         nodface,nactdof,mi,nodedesiinv,&iregion,
			         iponoelfa,ndesi,nodedesi,iponod2dto3d,
			         ikboun,nboun,ne2d)); 		          

    SFREE(konfa);SFREE(ipkonfa);SFREE(lakonfa);SFREE(iponoelfa);SFREE(inoelfa);
    
    for(i=0;i<*nobject;i++){
      if(strcmp1(&objectset[i*405+3],"MEMBERSIZE")==0){
        
	printf("Computation of geometric constraint MAXMEMEBERSIZE/MINMEMBERSIZE\n\n");
	
        iobject=i+1;
        thicknessmain(co,nobject,nk,nodedesi,ndesi,objectset,set,nset,
                      istartset,iendset,ialset,&iobject,nodedesiinv,
                      dgdxglob,extnor,coini,g0); 

        ++*kode;
        frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,
		&mode,          
                &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,
		iendset,
                ialset,jobnamec,output,dgdxglob,&i,objectset,ntrans,inotr,
                trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,
		&ifeasd);

      }else if(strcmp1(&objectset[i*405],"PACKAGING")==0){
        
	printf("Computation of geometric constraint PACKAGING\n\n");
	
        iobject=i+1;
        packagingmain(co,nobject,nk,nodedesi,ndesi,objectset,set,nset,
                      istartset,iendset,ialset,&iobject,nodedesiinv,
                      dgdxglob,extnor,g0); 

        ++*kode;
        frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,
		&mode,          
                &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,
		iendset,
                ialset,jobnamec,output,dgdxglob,&i,objectset,ntrans,inotr,
                trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,
		&ifeasd);

      }else if((strcmp1(&objectset[i*405],"MAXGROWTH")==0)||
               (strcmp1(&objectset[i*405],"MAXSHRINKAGE")==0)){

        printf("Computation of geometric constraint MAXGROWTH/MAXSHRINKAGE\n\n");
      
        iobject=i+1;
        FORTRAN(maxdesvardisp,(nobject,nk,set,nset,istartset,iendset,ialset,
                               &iobject,nodedesiinv,dgdxglob,objectset,xdesi,
		               coini,co,nodedesipos,ndesi,nodedesi,g0,
			       extnorini)); 

        ++*kode;
        frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,
		&mode,          
                &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,
		iendset,
                ialset,jobnamec,output,dgdxglob,&i,objectset,ntrans,inotr,
                trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,
		&ifeasd);
		      
      }
    } 
    
    SFREE(nodedesipos);SFREE(extnorini);SFREE(extnor);
  }
  
  /* bring the sensitivities in the field dgdxglob in the right order
     --> sensitivity of the objective function always at first */
  
  FORTRAN(clonesensitivities,(nobject,nk,objectset,g0,dgdxglob));
  
  /* min or max optimization target 
     sensitivities always point in this driection that objective function 
     will be maximized
     --> target = 'MIN' --> sensitivities of the objective function will be 
     scaled with -1 */
  
  iscaleflag=3; 
  FORTRAN(scalesen,(dgdxglob,feasdir,nk,nodedesi,ndesi,objectset,&iscaleflag,
		    &istart));

  /* check number of constraints */

  for(i=0;i<*nobject;i++){
    if((strcmp1(&objectset[i*405+18],"LE")==0)||
       (strcmp1(&objectset[i*405+18],"GE")==0)){
      nconstraint++;
    }
  }
  
  /* ********************************************************************* */
  /*         GRADIENT DESCENT AKIN METHOD                                  */
  /* ********************************************************************* */
  
  if(methodfd==1){
    
    printf("Computation of feasible direction with Gradient Descent Akin method\n\n");
    
    /* Scaling all sensitivities to unit length */

    iscaleflag=1;   
    for(i=0;i<*nobject;i++){
      istart=i+1;
      FORTRAN(scalesen,(dgdxglob,feasdir,nk,nodedesi,ndesi,objectset
			,&iscaleflag,
                        &istart));
    }

    /* determining the assembled gradient vector */
    
    printf("Assemlby of combined constraint gradient vector\n\n");
    
    NNEW(nodedesiboun,ITG,*ndesi);
    FORTRAN(constassembly,(nobject,objectset,g0,ndesi,dgdxglob,nk,nodedesi,
                           gradproj,set,nset,nodedesiboun,istartset,iendset,
			   ialset,nodedesiinv));
    
    /* assembly of feasible direction */

    printf("Computation of feasible direction\n\n");
    
    FORTRAN(calcfeasibledirection_gd,(ndesi,nodedesi,dgdxglob,nactive,nobject,
				      nk,gradproj,objectset));

    if(addout==1){
      ifeasd=1;            
      ++*kode;
      ipos=0;
      frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,
	      &mode,          
              &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
              ialset,jobnamec,output,gradproj,&ipos,gradprojname,ntrans,inotr,
              trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,&ifeasd);
    }


    /*************************************************************************/
    /*        CALCULATION OF PROJECTED GRADIENT                              */
    /*************************************************************************/

  }else if(methodfd==2){                                    

    printf("Computation of feasible direction with Gradient Projection method\n\n");
    
    /* in the field "ipoacti" and the variable "nnlconst" the 
       characteristics of the N-matrix are stored:
       -nnlconst is the number of nonlinear constraints
       -in the first nnlconst entries of ipoacti the position of the 
       nonlinear constraints (e.g. mass) is inserted
       -in the following entries the position of the linear 
       constraints (wall thickness constraints and fixation of nodes) 
       in nodedesi(i) is stored (i=ipoacti(j))
       -the field "iconstacti" contains the information whether the  
       constraint is LE or GE 
       > -1 is LE
       >  1 is GE
       -in the field "inameacti" the name of the active constraint is 
       saved via a pointer to the field objectset */
  
    NNEW(objnorm,double,*nobject**ndesi);
    NNEW(ipoacti,ITG,*nobject**ndesi);
    NNEW(inameacti,ITG,*nobject**ndesi);  
    NNEW(iconstacti,ITG,*nobject**ndesi);

    /* estimate number of active constraints (nactive) on the basis of 
       the function values of the constraints */
        
    FORTRAN(checkconstraint,(nobject,objectset,g0,&nactive,&nnlconst,ipoacti,
                             ndesi,dgdxglob,nk,nodedesi,iconstacti,objnorm,
                             inameacti));   

    RENEW(objnorm,double,nactive);
    RENEW(ipoacti,ITG,nactive);
    RENEW(inameacti,ITG,nactive);
    RENEW(iconstacti,ITG,nactive);
        
    gradientprojection(nobject,objectset,dgdxglob,g0,ndesi,nodedesi,nk,
        	       isolver,set,&nset,istartset,iendset,ialset,gradproj,
		       gradprojname,&nactive,objnorm,ipoacti,iconstacti,
		       inameacti,&nnlconst);	                   

    if(addout==1){
      ifeasd=2;            
      ++*kode;
      ipos=0;
      frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,
	      &mode,          
              &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
              ialset,jobnamec,output,gradproj,&ipos,gradprojname,ntrans,inotr,
              trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,&ifeasd);
    }
 
  }


  /***************************************************************************/
  /*        FORWARD FILTERING OF SENSITIVITIES                                /
  /***************************************************************************/

  printf("Forward filtering of feasible direction\n\n");
  
  /*distmin is set to 1.0e-08 --> perturabtion which reduces the
    finite difference error the most */
  distmin=1.0e-08;
  
  filtermain_forward(co,gradproj,nk,nodedesi,ndesi,objectset,xdesi,
	             &distmin,feasdir);

  /* scaling the designnodes being in the transition between 
     the designspace and the non-designspace */
       
  transitionmain(co,feasdir,nobject,nk,nodedesi,ndesi,objectset,
        	 ipkon,kon,lakon,ipoface,nodface,nodedesiinv);
  
  printf("Normalization of forward filtered sensitivities feasdir\n\n");

  iscaleflag=4; 
  FORTRAN(scalesen,(dgdxglob,feasdir,nk,nodedesi,ndesi,objectset,&iscaleflag,
		    &istart));

  ifeasd=3;
  ++*kode;
  ipos=0;
  frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,&mode,
  	  &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
  	  ialset,jobnamec,output,feasdir,&ipos,gradprojname,ntrans,inotr,
  	  trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,&ifeasd);

  SFREE(inum);SFREE(objnorm);SFREE(ipoacti);SFREE(iconstacti);
  SFREE(nodedesiboun);
  SFREE(inameacti);SFREE(nodedesiinv);SFREE(gradproj);SFREE(gradprojname);
  SFREE(feasdir);SFREE(ipoface);SFREE(nodface);
  
  *objectsetp=objectset;*dgdxglobp=dgdxglob;*ipkonp=ipkon;*lakonp=lakon;
  *konp=kon;*ielmatp=ielmat;
 
  return;
  
} 

#endif
