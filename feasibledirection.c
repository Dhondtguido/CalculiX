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
                     ITG *mi,ITG **ielmatp,ITG *ielprop,double *prop,ITG *kode,
                     ITG *nmethod,char *filab,ITG *nstate_,ITG *istep,
                     double *cs,char *set,ITG *nset,ITG *istartset,ITG *iendset,
                     ITG *ialset,char *jobnamec,char *output,ITG *ntrans,
                     ITG *inotr,double *trab,char *orname,double *xdesi){
               
  /* finding a feasible direction based on the sensitivity information */
  
  char *objectset=NULL,cflag[1]=" ",description[13]="            ",*lakon=NULL;
  
  ITG nzss,*mast1=NULL,*irows=NULL,*icols=NULL,*jqs=NULL,*ipointer=NULL,
    symmetryflag=0,inputformat=0,i,iconst,iter,nactiveold,*ipoactiold=NULL,
    *iconstactiold=NULL,*ipoacti=NULL,*iconstacti=NULL,nactive=0,nnlconst,
    iscaleflag,nrhs=1,*inameacti=NULL,*inameactiold=NULL,nconstraint=0,
    *inum=NULL,iinc=1,mode=-1,noddiam=-1,ngraph=1,idesvar=0,inorm=0,irand=0,
    ishape=0,nfield,iforce,icoordinate=1,*ipkon=NULL,iobject,istart,
    *kon=NULL,*ielmat=NULL,*nodedesiinv=NULL;
             
  double *au=NULL,*ad=NULL,*adb=NULL,*aub=NULL,sigma=0,*rhs=NULL,
    *vector=NULL,*xlambd=NULL,*xtf=NULL,*objnorm=NULL,*dgdxglob=NULL,
    *stn=NULL,ptime=0.;  
  
  objectset=*objectsetp;dgdxglob=*dgdxglobp;ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;ielmat=*ielmatp;
  
  NNEW(nodedesiinv,ITG,*nk);
  for(i=0;i<*ndesi;i++){
     nodedesiinv[nodedesi[i]-1]=1;
  }
  
  /* assessment of geometrical constraints  */

  /* createinum is called in order to determine the nodes belonging
     to elements; this information is needed in frd_se */

  NNEW(inum,ITG,*nk);
  FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
                        nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
                        ielprop,prop));
  
  for(i=0;i<*nobject;i++){
    if(strcmp1(&objectset[i*405+3],"MEMBERSIZE")==0){

      iobject=i+1;
      thicknessmain(co,nobject,nk,nodedesi,ndesi,objectset,
                    ipkon,kon,lakon,set,nset,istartset,iendset,ialset,
                    &iobject,nodedesiinv,dgdxglob,xdesi); 

      ++*kode;
      nfield=2;
      iforce=0;
      frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,&mode,          
              &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
              ialset,jobnamec,output,dgdxglob,&i,objectset,ntrans,inotr,
              trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape);

    }else if((strcmp1(&objectset[i*405],"FIXGROWTH")==0)||
             (strcmp1(&objectset[i*405],"FIXSHRINKAGE")==0)){
      
      iobject=i+1;
      FORTRAN(fixnode,(nobject,nk,set,nset,istartset,iendset,ialset,
                       &iobject,nodedesiinv,dgdxglob,objectset)); 

      ++*kode;
      nfield=2;
      iforce=0;
      frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,&mode,          
              &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
              ialset,jobnamec,output,dgdxglob,&i,objectset,ntrans,inotr,
              trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape);
   
    }
    
  }
  
  
  FORTRAN(clonesensitivities,(nobject,nk,objectset,g0,dgdxglob));
  
  /* check if any constraints are active */

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

  /* check number of constraints */

  for(i=0;i<*nobject;i++){
    if((strcmp1(&objectset[i*405+18],"LE")==0)||
      (strcmp1(&objectset[i*405+18],"GE")==0)){
      nconstraint++;
    }
  }
  
  NNEW(objnorm,double,*nobject);
  NNEW(ipoacti,ITG,*nobject**ndesi);
  NNEW(inameacti,ITG,*nobject**ndesi);  
  NNEW(iconstacti,ITG,*nobject**ndesi);

  /* estimate nactive on the basis of the function values of the 
     constraints */
     
  if(nconstraint>0){
    
    FORTRAN(checkconstraint,(nobject,objectset,g0,&nactive,&nnlconst,ipoacti,
                             ndesi,dgdxglob,nk,nodedesi,iconstacti,objnorm,
                             inameacti));   
  }

  RENEW(ipoacti,ITG,nactive);
  RENEW(inameacti,ITG,nactive);
  RENEW(iconstacti,ITG,nactive);
     
  if(nactive>0){
     
    iscaleflag=1;   
    for(i=0;i<*nobject;i++){
      istart=i+1;
      FORTRAN(scalesen,(dgdxglob,nobject,nk,nodedesi,ndesi,objectset,
                      &iscaleflag,&istart));
    }
                         
    *nobject=*nobject+1; 
    RENEW(dgdxglob,double,2**nk**nobject);
    RENEW(objectset,char,405**nobject);
     
    nactiveold=nactive+1;
    iter=0;

    while((nactive<nactiveold)&&(nactive>0)){
  
      nactiveold=nactive;
      iter=iter+1;
     
      /* determining the structure of the N-Matrix */

      nzss=20000000;
      NNEW(mast1,ITG,nzss);
      NNEW(irows,ITG,1);
      NNEW(icols,ITG,nactive);
      NNEW(jqs,ITG,nactive+1);
      NNEW(ipointer,ITG,nactive);
      NNEW(rhs,double,nactive);
      NNEW(xlambd,double,nactive);
      NNEW(xtf,double,nactive);
      NNEW(vector,double,*ndesi);
     
      mastructnmatrix(icols,jqs,&mast1,&irows,ipointer,&nzss,&nactive,&nnlconst);

      RENEW(irows,ITG,nzss);
      SFREE(ipointer);
     
      /* determining the entries of the N-Matrix */
  
      NNEW(ad,double,nactive);
      NNEW(au,double,nzss);    
  
      FORTRAN(nmatrix,(ad,au,jqs,irows,ndesi,nodedesi,dgdxglob,&nactive,
                       nobject,&nnlconst,ipoacti,nk));

      /* Calculate inverse of the N-matrix */

      NNEW(adb,double,nactive);
      NNEW(aub,double,nzss); 

      for(i=0;i<nactive;i++){
        adb[i]=1.0;
      }
       
      /* LU decomposition of the left hand matrix */

      if(*isolver==0){
#ifdef SPOOLES
        spooles_factor(ad,au,adb,aub,&sigma,icols,irows,&nactive,&nzss,
                       &symmetryflag,&inputformat,&nzss);
#else
        printf(" *ERROR in projectgrad: the SPOOLES library is not linked\n\n");
        FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
        token=1;
        sgi_factor(ad,au,adb,aub,&sigma,icols,irows,&nactive,&nzss,token);
#else
        printf(" *ERROR in projectgrad: the SGI library is not linked\n\n");
        FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
        tau_factor(ad,&au,adb,aub,&sigma,icols,&irows,&nactive,&nzss);
#else
        printf(" *ERROR in projectgrad: the TAUCS library is not linked\n\n");
        FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
        pardiso_factor(ad,au,adb,aub,&sigma,icols,irows,&nactive,&nzss,
                       &symmetryflag,&inputformat,jqs,&nzss);
#else
        printf(" *ERROR in projectgrad: the PARDISO library is not linked\n\n");
        FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
        pastix_factor_main(ad,au,adb,aub,&sigma,icols,irows,&nactive,&nzss,
                      &symmetryflag,&inputformat,jqs,&nzss);
#else
        printf(" *ERROR in projectgrad: the PASTIX library is not linked\n\n");
        FORTRAN(stop,());
#endif
      }
        
      /* solve the system nactive-times */

      for(i=0;i<nactive;i++){
        xtf[i]=0.00;
      }
        
      FORTRAN(preprojectgrad,(vector,ndesi,nodedesi,dgdxglob,&nactive,nobject,
			      &nnlconst,ipoacti,nk,rhs,&iconst,objectset,xtf));
        
      /* Calculate the projected gradient and the lagrange multipliers */
        
      for(iconst=1;iconst<=nactive;iconst++){
            
        for(i=0;i<nactive;i++){
          rhs[i]=0.00;
        }
            
        rhs[iconst-1]=1.0;
            
        /* solve the system */
            
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(rhs,&nactive);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
          sgi_solve(rhs,token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau_solve(rhs,&nactive);
#endif
        }
        else if(*isolver==7){
#ifdef PARDISO
          pardiso_solve(rhs,&nactive,&symmetryflag,&inputformat,&nrhs);
#endif
        }
        else if(*isolver==8){
#ifdef PASTIX
          pastix_solve(rhs,&nactive,&symmetryflag,&nrhs);
#endif
        }
            
        for(i=0;i<*ndesi;i++){
          vector[i]=0.00;
        }
            
        /* carry out matrix multiplications */
           
        FORTRAN(projectgrad,(vector,ndesi,nodedesi,dgdxglob,&nactive,nobject,
                             &nnlconst,ipoacti,nk,rhs,&iconst,objectset,xlambd,
			     xtf,objnorm));
      }
        
      /* clean the system */
        
      if(*isolver==0){
#ifdef SPOOLES
        spooles_cleanup();
#endif
      }
      else if(*isolver==4){
#ifdef SGI
        sgi_cleanup(token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
        tau_cleanup();
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
        pardiso_cleanup(&nactive,&symmetryflag,&inputformat);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
#endif
      }
        
      /* write the results of the iteration in the dat-file */
        
      FORTRAN(writelm,(&iter,xlambd,&nactive,&nnlconst,objectset,nobject,
                       ipoacti,iconstacti,inameacti,nodedesi,dgdxglob,nk));
        
      /* check if the langrange multipliers of the active constraints
      have the "correct" sign ("correct" in the sense of active constraints) */
                 
      NNEW(ipoactiold,ITG,nactive);
      NNEW(iconstactiold,ITG,nactive);
      NNEW(inameactiold,ITG,nactive);
      
      FORTRAN(checkprojectgrad,(&nactiveold,&nactive,ipoacti,ipoactiold,
				objectset,xlambd,&nnlconst,iconstacti,
				iconstactiold,inameacti,inameactiold,g0,
				nobject,ndesi,nodedesi,dgdxglob,nk));

      SFREE(mast1);SFREE(irows);SFREE(icols);SFREE(jqs);
      SFREE(ad);SFREE(au);SFREE(adb);SFREE(aub);SFREE(rhs);SFREE(xlambd);
      SFREE(xtf);SFREE(vector);SFREE(ipoactiold);SFREE(iconstactiold);
      SFREE(inameactiold);
        
    }
     
  }
  
  /* normalization of projected gradient and documentation of results */
       
  FORTRAN(postprojectgrad,(ndesi,nodedesi,dgdxglob,&nactive,nobject,
                           &nnlconst,ipoacti,nk,objectset,inameacti));
                

  /* scaling the sensitivities: highest absolute value is scaled to 1 */
  /* storing the feasible direction in the frd file */

  if(nactive>0){
    iobject=*nobject;
  }else{
    iobject=1;
  }

  iscaleflag=2;
  FORTRAN(scalesen,(dgdxglob,nobject,nk,nodedesi,ndesi,objectset,
                    &iscaleflag,&iobject));
  ++*kode;
  nfield=2;
  iforce=0;
  if(nactive>0){
    iobject=*nobject-1;
  }else{
    iobject=0;
  }
    
  frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,&mode,          
          &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
          ialset,jobnamec,output,dgdxglob,&iobject,objectset,ntrans,inotr,
          trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape);
  
  SFREE(inum);
  SFREE(objnorm);
  SFREE(ipoacti);
  SFREE(iconstacti);
  SFREE(inameacti);
  SFREE(nodedesiinv);
  
  *objectsetp=objectset;*dgdxglobp=dgdxglob;*ipkonp=ipkon;*lakonp=lakon;
  *konp=kon;*ielmatp=ielmat;
  
  return;
  
} 

#endif
