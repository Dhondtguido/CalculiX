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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "CalculiX.h"

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
		       double *timepar){
               
  /* finding a feasible direction based on the sensitivity information */
  
  char *objectset=NULL,cflag[1]=" ",description[13]="            ",*lakon=NULL,
    *gradprojname=NULL,*lakonfa=NULL;
  
  ITG i,*ipoacti=NULL,*iconstacti=NULL,nactive=0,nnlconst,iscaleflag,
    *inameacti=NULL,nconstraint=0,*inum=NULL,iinc=1,mode=-1,noddiam=-1,ngraph=1,
    idesvar=0,inorm=0,irand=0,ishape=0,icoordinate=1,*ipkon=NULL,
    iobject,istart,*kon=NULL,*ielmat=NULL,*nodedesiinv=NULL,nev,ncv,
    *select=NULL,ido,ldz,iparam[11],info,mei[4],lworkl,ipntr[14],jrow,rvec=1,
    nzss,*mast1=NULL,*irows=NULL,*icols=NULL,*jqs=NULL,*ipointer=NULL,
    symmetryflag=0,inputformat=0,nrhs=1,*istartdesi=NULL,*ialdesi=NULL,
    iregion=0,*iponoel=NULL,*inoel=NULL,*ipoface=NULL,*nodface=NULL,*konfa=NULL,
    *ipkonfa=NULL,nsurfs,ifreemax,ifeasd,fdmethod;
       
  double *objnorm=NULL,*dgdxglob=NULL,*stn=NULL,ptime=0.,distmin,
    *mshupdate=NULL,
    *gradproj=NULL,*constgrad=NULL,pi,mxiter,tol,trace=0,*resid=NULL,*zz=NULL,
    *workd=NULL,*workl=NULL,*temp_array=NULL,*d=NULL,*stdesc=NULL,*au=NULL,
    *ad=NULL,*adb=NULL,*aub=NULL,sigma=0,*rhs=NULL,*weightformgrad=NULL;   
  
  objectset=*objectsetp;dgdxglob=*dgdxglobp;ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;ielmat=*ielmatp;

  fdmethod=timepar[3];
  fdmethod=2;
  NNEW(gradproj,double,3**nk);
  NNEW(gradprojname,char,405);
  NNEW(nodedesiinv,ITG,*nk);
  for(i=0;i<*ndesi;i++){
    nodedesiinv[nodedesi[i]-1]=1;
  }
  
  /* createinum is called in order to determine the nodes belonging
     to elements; this information is needed in frd_se */

  NNEW(inum,ITG,*nk);
  FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
		      nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,
		      ielmat,
		      ielprop,prop));

  /****************************************************************************/
  /* assessment of geometrical constraints                                     /
  /****************************************************************************/
  
  for(i=0;i<*nobject;i++){
    if(strcmp1(&objectset[i*405+3],"MEMBERSIZE")==0){

      iobject=i+1;
      thicknessmain(co,nobject,nk,nodedesi,ndesi,objectset,
                    ipkon,kon,lakon,set,nset,istartset,iendset,ialset,
                    &iobject,nodedesiinv,dgdxglob,xdesi); 

      ++*kode;
      frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,
	      &mode,          
              &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
              ialset,jobnamec,output,dgdxglob,&i,objectset,ntrans,inotr,
              trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,&ifeasd);

    }else if((strcmp1(&objectset[i*405],"FIXGROWTH")==0)||
             (strcmp1(&objectset[i*405],"FIXSHRINKAGE")==0)){
      
      iobject=i+1;
      FORTRAN(fixnode,(nobject,nk,set,nset,istartset,iendset,ialset,
                       &iobject,nodedesiinv,dgdxglob,objectset)); 

      ++*kode;
      frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,
	      &mode,          
              &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
              ialset,jobnamec,output,dgdxglob,&i,objectset,ntrans,inotr,
              trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,&ifeasd);
    } 
  }
  
  FORTRAN(clonesensitivities,(nobject,nk,objectset,g0,dgdxglob));
  
  /* min or max optimization target */
  
  /*iscaleflag=3;   
    FORTRAN(scalesen,(dgdxglob,nk,nodedesi,ndesi,objectset,&iscaleflag,&istart));*/

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
  
  if((fdmethod==1)&&(nconstraint>0)){

    fdmethod=1;


    /**************************************************************************/
    /*        CALCULATION OF PROJECTED GRADIENT                               */
    /**************************************************************************/

  }else if((fdmethod==2)&&(nconstraint>0)){                                    

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
    
    ifeasd=2;            
    ++*kode;
    *nobject=0;

    frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,&iinc,
	    &mode,          
            &noddiam,description,mi,&ngraph,ne,cs,set,nset,istartset,iendset,
            ialset,jobnamec,output,gradproj,nobject,gradprojname,ntrans,inotr,
            trab,&idesvar,orname,&icoordinate,&inorm,&irand,&ishape,&ifeasd);

  }

  SFREE(mshupdate);SFREE(inum);SFREE(objnorm); SFREE(ipoacti);SFREE(iconstacti);
  SFREE(inameacti);SFREE(nodedesiinv);SFREE(gradproj);SFREE(gradprojname);
  SFREE(stdesc);
  
  *objectsetp=objectset;*dgdxglobp=dgdxglob;*ipkonp=ipkon;*lakonp=lakon;
  *konp=kon;*ielmatp=ielmat;
 
  return;
  
} 
