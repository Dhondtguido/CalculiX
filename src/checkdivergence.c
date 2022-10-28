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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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

void checkdivergence(double *co, ITG *nk, ITG *kon, ITG *ipkon, char *lakon,
	  ITG *ne, double *stn, ITG *nmethod, 
	  ITG *kode, char *filab, double *een, double *t1act,
          double *time, double *epn,ITG *ielmat,char *matname,
          double *enern, double *xstaten, ITG *nstate_, ITG *istep,
          ITG *iinc, ITG *iperturb, double *ener, ITG *mi, char *output,
          ITG *ithermal, double *qfn, ITG *mode, ITG *noddiam, double *trab,
          ITG *inotr, ITG *ntrans, double *orab, ITG *ielorien, ITG *norien,
          char *description,double *sti,
	  ITG *icutb, ITG *iit, double *dtime, double *qa, double *vold,
          double *qam, double *ram1, double *ram2, double *ram,
          double *cam, double *uam, ITG *ntg, double *ttime,
          ITG *icntrl, double *theta, double *dtheta, double *veold,
          double *vini, ITG *idrct, double *tper,ITG *istab, double *tmax, 
          ITG *nactdof, double *b, double *tmin, double *ctrl, double *amta,
          ITG *namta, ITG *itpamp, ITG *inext, double *dthetaref, ITG *itp,
          ITG *jprint, ITG *jout, ITG *uncoupled, double *t1, ITG *iitterm,
          ITG *nelemload, ITG *nload, ITG *nodeboun, ITG *nboun, ITG *itg,
          ITG *ndirboun, double *deltmx, ITG *iflagact,char *set,ITG *nset,
	  ITG *istartset,ITG *iendset,ITG *ialset, double *emn, double *thicke,
	  char *jobnamec,ITG *mortar,ITG *nmat,ITG *ielprop,double *prop,
	  ITG *ialeatoric,ITG *kscale,
	  double *energy, double *allwk, double *energyref, double *emax,
	  double *enres, double *enetoll, double *energyini, double *allwkini,
	  double *temax, double *sizemaxinc, ITG* ne0, ITG* neini,
	  double *dampwk, double *dampwkini, double *energystartstep) {

    ITG ia,ngraph=1,k,*ipneigh=NULL,*neigh=NULL,*inum=NULL,mt=mi[1]+1,kscalemax;

    double *vr=NULL,*vi=NULL,*stnr=NULL,
	*stni=NULL,*vmax=NULL,*stnmax=NULL,*cs=NULL,
        *fn=NULL,*eenmax=NULL,*fnr=NULL,*fni=NULL,*qfx=NULL,*cdn=NULL,
        *cdnr=NULL,*cdni=NULL;

    ia=ctrl[7];kscalemax=ctrl[55];

    /* check whether divergence was signaled in radflowload.c
       => repeat the increment with a smaller size */

	*dtheta=*dtheta*qa[2];
	printf("increment size decrease requested by the network\n\n");
	*dthetaref=*dtheta;
	printf(" divergence; the increment size is decreased to %e\n",*dtheta**tper);
	printf(" the increment is reattempted\n\n");
	
	FORTRAN(writestadiv,(istep,iinc,icutb,iit,ttime,time,
				 dtime));
	
	*istab=0;
	if(*itp==1){
	    *itp=0;
	    (*inext)--;
	}
	
	/* check whether new increment size is smaller than minimum */
	
	if(*dtheta<*tmin){
	    printf("\n *ERROR: increment size smaller than minimum\n");
	    printf(" best solution and residuals are in the frd file\n\n");
	    NNEW(fn,double,mt**nk);
	    NNEW(inum,ITG,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	    FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,
				   nk,sti,stn,ipkon,inum,kon,lakon,ne,mi,orab,
				   ielorien,co,itg,ntg,vold,ielmat,thicke,ielprop,prop));
	    ++*kode;
	    
	    (*ttime)+=(*time);
	    frd(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,nmethod,
		kode,filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
		xstaten,nstate_,istep,iinc,ithermal,qfn,mode,noddiam,
		trab,inotr,ntrans,orab,ielorien,norien,description,
		ipneigh,neigh,mi,sti,vr,vi,stnr,stni,vmax,stnmax,
		&ngraph,veold,ener,ne,cs,set,nset,istartset,iendset,
		ialset,eenmax,fnr,fni,emn,thicke,jobnamec,output,qfx,cdn,
		mortar,cdnr,cdni,nmat,ielprop,prop,sti);
	    
	    FORTRAN(stop,());
	}
	*icntrl=1;
	(*icutb)++;
	if(*mortar==1){
	    *kscale=kscalemax;
	    printf("\n reducing the constant stiffnesses by a factor of %d \n\n",*kscale);
	}
	
	/* check whether too many cutbacks */
	
	if(*icutb>ia){
	    printf("\n *ERROR: too many cutbacks\n");
	    printf(" best solution and residuals are in the frd file\n\n");
	    NNEW(fn,double,mt**nk);
	    NNEW(inum,ITG,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	    FORTRAN(storeresidual,(nactdof,b,fn,filab,ithermal,
				   nk,sti,stn,ipkon,inum,kon,lakon,ne,mi,orab,
				   ielorien,co,itg,ntg,vold,ielmat,thicke,ielprop,prop));
	    ++*kode;
	    
	    (*ttime)+=(*time);
	    frd(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,nmethod,
		kode,filab,een,t1act,fn,ttime,epn,ielmat,matname,enern,
		xstaten,nstate_,istep,iinc,ithermal,qfn,mode,noddiam,
		trab,inotr,ntrans,orab,ielorien,norien,description,
		ipneigh,neigh,mi,sti,vr,vi,stnr,stni,vmax,stnmax,
		&ngraph,veold,ener,ne,cs,set,nset,istartset,iendset,
		ialset,eenmax,fnr,fni,emn,thicke,jobnamec,output,qfx,cdn,
		mortar,cdnr,cdni,nmat,ielprop,prop,sti);
	    
	    FORTRAN(stop,());
	}
	if(*uncoupled){
	    if(*ithermal==1){
		(ctrl[0])/=4;
	    }
	    *ithermal=3;
	}
	
	/* default value for qa[2] */
	
	qa[2]=-1.;
	
	*iflagact=0;

return;
}
