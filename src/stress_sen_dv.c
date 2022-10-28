/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "CalculiX.h"

void stress_sen_dv(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
		   double *dstn,double *elcon,ITG *nelcon,double *rhcon,
		   ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
		   ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,
		   ITG *ntmat_,double *t0,double *t1,ITG *ithermal,
		   double *prestr,ITG *iprestr,char *filab,ITG *iperturb,
		   double *dv, ITG *nmethod,double *dtime,double *time,
		   double *ttime,double *plicon,ITG *nplicon,double *plkcon,
		   ITG *nplkcon,double *xstateini,double *xstate,ITG *npmat_,
		   char *matname,ITG *mi,ITG *ielas,ITG *ncmat_,ITG *nstate_,
		   double *stiini,double *vini,double *emeini,double *enerini,
		   ITG *istep,ITG *iinc,double *springarea,double *reltime,
		   ITG *ne0,double *thicke,double *pslavsurf,double *pmastsurf,
		   ITG *mortar,double *clearini,ITG *ielprop,double *prop,
		   ITG *kscale,ITG *iobject,double *g0,ITG *nea,ITG *neb,
		   ITG *nasym,double *distmin,double *dstx,ITG *ialnk,
		   double *dgdu,ITG *ialeneigh,ITG *neaneigh,ITG *nebneigh,
		   ITG *ialnneigh,ITG *naneigh,ITG *nbneigh,double *stn,
		   double *expks,char *objectset,ITG *idof,ITG *node,ITG *idir,
		   double *vold,double *dispmin,double *physcon){         
            
  ITG symmetryflag=0,mt=mi[1]+1,i,iactpos,calcul_fn,list,
    calcul_qa,calcul_cauchy,ikin=0,nal,iout=2,icmd=3,nener=0,
    *inum=NULL,nprintl=0,unperturbflag,nfield,ndim,iorienglob,
    force,mscalmethod=0,*islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,
    mortartrafoflag=0,intscheme=0;

  double *fn=NULL,*eei=NULL,qa[4]={0.,0.,-1.,0.},*xstiff=NULL,*ener=NULL,    
    *eme=NULL,kspart,dksper,*smscale=NULL,enerscal=0.,*t0g=NULL,*t1g=NULL,
    *autloc=NULL;
  
  char cflag[1];
    
  if(*nasym!=0){symmetryflag=2;}
      
  NNEW(eme,double,6*mi[0]**ne);
  NNEW(inum,ITG,*nk);
  
  /* calculating the perturbed stresses */
  
  NNEW(fn,double,mt**nk);
  NNEW(eei,double,6*mi[0]**ne);

  /* setting the output variables */
  
  calcul_fn=0;
  calcul_qa=0;
  if(iperturb[1]==1){calcul_cauchy=1;}else{calcul_cauchy=0;}
  
  list=1;
  FORTRAN(resultsmech,(co,kon,ipkon,lakon,ne,dv,
		       dstx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		       ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
		       iprestr,eme,iperturb,fn,&iout,qa,dv,
		       nmethod,
		       dv,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
		       xstateini,xstiff,xstate,npmat_,matname,mi,ielas,&icmd,
		       ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
		       springarea,reltime,&calcul_fn,&calcul_qa,&calcul_cauchy,
		       &nener,&ikin,&nal,ne0,thicke,emeini,
		       pslavsurf,pmastsurf,mortar,clearini,nea,neb,ielprop,
		       prop,kscale,&list,ialnk,smscale,&mscalmethod,&enerscal,
		       t0g,t1g,islavelinv,autloc,irowtloc,jqtloc,
		       &mortartrafoflag,&intscheme,physcon));

  /* extrapolating the stresses */

  nfield=6;
  ndim=6;
  iorienglob=0;
  force=0;
  //  strcpy1(&cflag[0],&filab[2962],1);
  strcpy1(&cflag[0],&filab[178],1);
  
  FORTRAN(extrapolate_se,(dstx,dstn,ipkon,inum,kon,lakon,
			  &nfield,nk,ne,mi,&ndim,orab,ielorien,co,&iorienglob,
			  cflag,dv,&force,ielmat,thicke,ielprop,prop,ialeneigh,
			  neaneigh,nebneigh,ialnneigh,naneigh,nbneigh));

  /* Calculate KS-function and sensitivity */

  FORTRAN(objective_stress_se,(nk,iobject,mi,dstn,objectset,
			       ialnneigh,naneigh,nbneigh,stn,&dksper));

  dgdu[*idof]=dksper/(*dispmin**expks);
      
  SFREE(fn);SFREE(eei);SFREE(eme);SFREE(inum);

}
