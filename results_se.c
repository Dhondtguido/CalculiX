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

static char *lakon1,*matname1;

static ITG *kon1,*ipkon1,*ne1,*nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,
    *norien1,*ntmat1_,*ithermal1,*iprestr1,*iperturb1,*iout1,*nmethod1,
    *nplicon1,*nplkcon1,*npmat1_,*mi1,*ielas1,*icmd1,*ncmat1_,*nstate1_,
    *istep1,*iinc1,calcul_fn1,calcul_cauchy1,nener1,ikin1,
    num_cpus,mt1,*nk1,*ne01,*mortar1,*ielprop1,idesvar1,*nodedesi1,*nkon1,
    *icoordinate1,*ialdesi1,*neapar=NULL,*nebpar=NULL;

static double *co1,*v1,*stx1,*elcon1,*rhcon1,*alcon1,*alzero1,*orab1,*t01,*t11,
    *prestr1,*eme1,*fn1=NULL,*vold1,*veold1,*dtime1,*time1,
    *ttime1,*plicon1,*plkcon1,*xstateini1,*xstiff1,*xstate1,*stiini1,
    *vini1,*ener1,*eei1,*enerini1,*springarea1,*reltime1,
    *thicke1,*emeini1,*prop1,*dxstiff1,
    *pslavsurf1,*pmastsurf1,*clearini1,*dfn1,*fn01,
    *sti1,*xdesi1;

void results_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
       double *v,double *stn,ITG *inum,double *stx,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,
       double *t1,ITG *ithermal,double *prestr,ITG *iprestr,char *filab,
       double *eme,double *emn,
       double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,ITG *iout,
       double *qa,double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,double *veold,
       double *accold,double *bet,double *gam,double *dtime,double *time,
       double *ttime,double *plicon,ITG *nplicon,double *plkcon,
       ITG *nplkcon,double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
       double *epn,char *matname,ITG *mi,ITG *ielas,ITG *icmd,ITG *ncmat_,
       ITG *nstate_,
       double *stiini,double *vini,ITG *ikboun,ITG *ilboun,double *ener,
       double *enern,double *emeini,double *xstaten,double *eei,double *enerini,
       double *cocon,ITG *ncocon,char *set,ITG *nset,ITG *istartset,
       ITG *iendset,
       ITG *ialset,ITG *nprint,char *prlab,char *prset,double *qfx,double *qfn,
       double *trab,
       ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,ITG *nload,
       ITG *ikmpc,ITG *ilmpc,
       ITG *istep,ITG *iinc,double *springarea,double *reltime, ITG *ne0,
       double *xforc, ITG *nforc, double *thicke,
       double *shcon,ITG *nshcon,char *sideload,double *xload,
       double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
       double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
       ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
       ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
       double *energy,double *df,double *distmin,ITG *ndesi,ITG *nodedesi,
       double *sti,ITG *nkon,ITG *jqs,ITG *irows,
       ITG *nactdofinv,ITG *icoordinate,double *dxstiff,ITG *istartdesi,
       ITG *ialdesi,double *xdesi,ITG *ieigenfrequency,double *fint,
       ITG *ishapeenergy,char *typeboun){

    ITG intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,nener=0,ikin,
        intpointvart,mt=mi[1]+1,i,j,idesvar,iorien,idir,im,
        nea,neb;
    
    double *dfn=NULL,*fn0=NULL,a[9],pgauss[3],rotvec[3],orabsav[7];

    /* calculating the sensitivity of the internal forces */
      
    /* variables for multithreading procedure */
    
    ITG sys_cpus,*ithread=NULL;
    char *env,*envloc,*envsys;
    
    num_cpus = 0;
    sys_cpus=0;

    /* explicit user declaration prevails */

    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
        sys_cpus=atoi(envsys);
        if(sys_cpus<0) sys_cpus=0;
    }

    /* automatic detection of available number of processors */

    if(sys_cpus==0){
        sys_cpus = getSystemCPUs();
        if(sys_cpus<1) sys_cpus=1;
    }

    /* local declaration prevails, if strictly positive */

    envloc = getenv("CCX_NPROC_SENS");
    if(envloc){
        num_cpus=atoi(envloc);
        if(num_cpus<0){
            num_cpus=0;
        }else if(num_cpus>sys_cpus){
            num_cpus=sys_cpus;
        }
        
    }

    /* else global declaration, if any, applies */

    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
        if (env)
            num_cpus = atoi(env);
        if (num_cpus < 1) {
            num_cpus=1;
        }else if(num_cpus>sys_cpus){
            num_cpus=sys_cpus;
        }
    }

// next line is to be inserted in a similar way for all other paralell parts

    if(*ne<num_cpus) num_cpus=*ne;
    
    pthread_t tid[num_cpus];
    
    /* 1. nodewise storage of the primary variables
       2. determination which derived variables have to be calculated */

    resultsini(nk,v,ithermal,filab,iperturb,f,fn,
	       nactdof,iout,qa,vold,b,nodeboun,ndirboun,
	       xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,
	       veold,accold,bet,gam,dtime,mi,vini,nprint,prlab,
	       &intpointvarm,&calcul_fn,&calcul_f,&calcul_qa,&calcul_cauchy,
	       &ikin,&intpointvart,typeboun,&num_cpus,mortar,&nener);

    NNEW(fn0,double,mt**nkon);
    NNEW(dfn,double,mt**nk);
	    
    if(((*nmethod!=4)&&(*nmethod!=5))||(iperturb[0]>1)){
	printf(" Using up to %" ITGFORMAT " cpu(s) for the sensitivity of the internal forces.\n\n", num_cpus);
    }

    /* nodal forces without perturbation */

    idesvar=0;
    
    /* calculating the stresses and material tangent at the 
       integration points; calculating the internal forces */
    
    if(((ithermal[0]<=1)||(ithermal[0]>=3))&&(intpointvarm==1)){
    
        /* determining the element bounds in each thread */

	NNEW(neapar,ITG,num_cpus);
	NNEW(nebpar,ITG,num_cpus);
	elementcpuload(neapar,nebpar,ne,ipkon,&num_cpus);
	
	NNEW(fn01,double,num_cpus*mt**nkon);
	
	co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;v1=v;
	stx1=stx;elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;
	nrhcon1=nrhcon;alcon1=alcon;nalcon1=nalcon;alzero1=alzero;
	ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;
	ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;prestr1=prestr;
	iprestr1=iprestr;eme1=eme;iperturb1=iperturb;iout1=iout;
	vold1=vold;nmethod1=nmethod;veold1=veold;dtime1=dtime;
	time1=time;ttime1=ttime;plicon1=plicon;nplicon1=nplicon;
	plkcon1=plkcon;nplkcon1=nplkcon;xstateini1=xstateini;
	xstiff1=xstiff;xstate1=xstate;npmat1_=npmat_;matname1=matname;
	mi1=mi;ielas1=ielas;icmd1=icmd;ncmat1_=ncmat_;nstate1_=nstate_;
	stiini1=stiini;vini1=vini;ener1=ener;eei1=eei;enerini1=enerini;
	istep1=istep;iinc1=iinc;springarea1=springarea;reltime1=reltime;
	calcul_fn1=calcul_fn;calcul_cauchy1=calcul_cauchy;
	nener1=nener;ikin1=ikin;mt1=mt;nk1=nk;ne01=ne0;thicke1=thicke;
	emeini1=emeini;pslavsurf1=pslavsurf;clearini1=clearini;
	pmastsurf1=pmastsurf;mortar1=mortar;ielprop1=ielprop;prop1=prop;
	idesvar1=idesvar;nodedesi1=nodedesi;
	sti1=sti;nkon1=nkon;icoordinate1=icoordinate;
	dxstiff1=dxstiff;ialdesi1=ialdesi;xdesi1=xdesi;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,num_cpus);
	for(i=0; i<num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)resultsmechmt_se, (void *)&ithread[i]);
	}
	for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	
	/* Assembling fn0 and dfn */
	
	for(i=0;i<mt**nkon;i++){
	    fn0[i]=fn01[i];
	}
	for(i=0;i<mt**nkon;i++){
	    for(j=1;j<num_cpus;j++){
		fn0[i]+=fn01[i+j*mt**nkon];
	    }
	}
	
	SFREE(fn01);SFREE(ithread);SFREE(neapar);SFREE(nebpar); 
    }
    
    /* in case of nonlinear geometry calculate vector fint */
    
    if((iperturb[1]==1)&&(*ishapeenergy==1)){
       FORTRAN(createfint,(ne,ipkon,lakon,kon,nactdof,mi,fn0,fint));
    }
    
    /* loop over the design variables (perturbation) */
    
    for(idesvar=1;idesvar<=*ndesi;idesvar++){
	
	DMEMSET(dfn,0,mt**nk,0.);
	   
        /* calculate a delta in the orientation
           in case the material orientation is the design variable */
	
	if(*icoordinate!=1){
	    iorien=(idesvar-1)/3;
	    
	    /* save nominal orientation */
	    
	    memcpy(&orabsav[0],&orab[7*iorien],sizeof(double)*7);
	    
	    /* calculate the transformation matrix */
	    
	    FORTRAN(transformatrix,(&orab[7*iorien],pgauss,a));
	    
	    /* calculate the rotation vector from the transformation matrix */
	    
	    FORTRAN(rotationvector,(a,rotvec));
	    idir=(idesvar-1)-iorien*3;
	    
	    /* add a small variation to the rotation vector component */
	    
	    rotvec[idir]+=*distmin;
	    
	    /* determine the new transformation matrix */
	    
	    FORTRAN(rotationvectorinv,(a,rotvec));
	    
	    /* determine two new points in the x-y plane */
	    
	    for(i=0;i<6;i++){orab[7*iorien+i]=a[i];}
	}
	
	/* calculating the stresses and material tangent at the 
	   integration points; calculating the internal forces */
	
	if(((ithermal[0]<=1)||(ithermal[0]>=3))&&(intpointvarm==1)){
	    nea=istartdesi[idesvar-1];
	    neb=istartdesi[idesvar]-1;
	    
	    FORTRAN(resultsmech_se,(co,kon,ipkon,lakon,ne,v,
		stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		ielmat,ielorien,norien,orab,ntmat1_,t0,t1,ithermal,prestr,
		iprestr,eme,iperturb,fn,iout,vold,nmethod,
	        veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
                xstateini,xstiff,xstate,npmat1_,matname,mi,ielas,icmd,
                ncmat1_,nstate1_,stiini,vini,ener,eei,enerini,istep,iinc,
                springarea,reltime,&calcul_fn,&calcul_cauchy,&nener,
                &ikin,ne0,thicke,emeini,
                pslavsurf,pmastsurf,mortar,clearini,&nea,&neb,ielprop,prop,
                dfn,&idesvar,nodedesi,
	        fn0,sti,icoordinate,dxstiff,ialdesi,xdesi));
	}
	
	/* calculating the matrix system internal force vector
           for nonlinear geometrical calculations */
	
//	if((iperturb[1]==1)&&(*ieigenfrequency!=1)){
	if(*ieigenfrequency!=1){
	   FORTRAN(resultsforc_se,(nk,dfn,nactdofinv,ipompc,nodempc,
				coefmpc,nmpc,mi,fmpc,&calcul_fn,&calcul_f,
				&idesvar,df,jqs,irows,distmin));
	   }
	   
        /* restoring the nominal orientation (in case the design variables
           are the orientations */
	
	if(*icoordinate!=1){
	    if(idesvar>0){
		memcpy(&orab[7*iorien],&orabsav[0],sizeof(double)*7);
	    }
	}
	
    }   /* end loop over design variables */
    
    
    SFREE(fn0);SFREE(dfn);
    
    return;
    
}

/* subroutine for multithreading of resultsmech_se */

void *resultsmechmt_se(ITG *i){

    ITG indexfn0,indexdfn,nea,neb;

    if(idesvar1==0){
	indexfn0=*i*mt1**nkon1;
	indexdfn=0;
    }else{
	indexfn0=0;
	indexdfn=*i*mt1**nk1;
    }

/*    nedelta=(ITG)floor(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if((*i==num_cpus-1)&&(neb<*ne1)) neb=*ne1;*/

    nea=neapar[*i]+1;
    neb=nebpar[*i]+1;

    FORTRAN(resultsmech_se,(co1,kon1,ipkon1,lakon1,ne1,v1,
          stx1,elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
          ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,prestr1,
          iprestr1,eme1,iperturb1,fn1,iout1,vold1,
          nmethod1,
          veold1,dtime1,time1,ttime1,plicon1,nplicon1,plkcon1,nplkcon1,
          xstateini1,xstiff1,xstate1,npmat1_,matname1,mi1,ielas1,icmd1,
          ncmat1_,nstate1_,stiini1,vini1,ener1,eei1,enerini1,istep1,iinc1,
          springarea1,reltime1,&calcul_fn1,&calcul_cauchy1,&nener1,
          &ikin1,ne01,thicke1,emeini1,
          pslavsurf1,pmastsurf1,mortar1,clearini1,&nea,&neb,ielprop1,prop1,
          &dfn1[indexdfn],&idesvar1,nodedesi1,
	  &fn01[indexfn0],sti1,icoordinate1,dxstiff1,ialdesi1,xdesi1));

    return NULL;
}
