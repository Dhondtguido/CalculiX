
/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

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
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>

#include "CalculiX.h"
#include "readfrd.h"

//ITG nkold,mpcrfna,mpcrfnb;

static ITG *nkapar=NULL,*nkbpar=NULL,*integerglob1,nkold1,*nk1,*iprfn1,
  *konrfn1,ideltamax;

static double *co1,*doubleglob1,*ratiorfn1; 

void readnewmesh(char *jobnamec,ITG *nboun,ITG *nodeboun,ITG *iamboun,
		 double *xboun,ITG *nload,char *sideload,ITG *iamload,
		 ITG *nforc,ITG *nodeforc,
		 ITG *iamforc,double *xforc,ITG *ithermal,ITG *nk,
		 double **t1p,ITG **iamt1p,ITG *ne,char **lakonp,ITG **ipkonp,
		 ITG **konp,ITG *istartset,ITG *iendset,ITG *ialset,
		 char *set,ITG *nset,char *filab,double **cop,ITG **ipompcp,
		 ITG **nodempcp,double **coefmpcp,ITG *nmpc,ITG *nmpc_,
		 char **labmpcp,ITG *mpcfree,ITG *memmpc_,ITG **ikmpcp,
		 ITG **ilmpcp,ITG *nk_,ITG *ne_,ITG *nkon_,ITG *istep,
		 ITG *nprop_,ITG **ielpropp,ITG *ne1d,ITG *ne2d,ITG **iponorp,
		 double **thicknp,double **thickep,ITG *mi,double **offsetp,
		 ITG **iponoelp,ITG **rigp,ITG **ne2bounp,ITG **ielorienp,
		 ITG **inotrp,double **t0p,double **t0gp,double **t1gp,
		 double **prestrp,double **voldp,double **veoldp,ITG **ielmatp,
		 ITG *irobustdesign,ITG **irandomtypep,double **randomvalp,
		 ITG *nalset,ITG *nalset_,ITG *nkon,double *xnor,
		 ITG *iaxial,ITG *network,ITG *nlabel,ITG *iuel,ITG *iperturb,
		 ITG *iprestr,ITG *ntie,char *tieset,ITG **iparentelp,
		 ITG *ikboun,ITG *ifreebody,ITG **ipobodyp,ITG *nbody,
		 ITG **iprfnp,ITG **konrfnp,double **ratiorfnp,ITG *nodempcref,
		 double *coefmpcref,ITG *memmpcref_,ITG *mpcfreeref,
		 ITG *maxlenmpcref,ITG *maxlenmpc,ITG *norien,double *tietol,
		 ITG *ntrans,ITG *nam)
{

  char masterrfnfile[132]="",fnrfn[132]="",*inpc=NULL,*labmpc=NULL,*lakon=NULL,
    masterurffile[132]="";

  ITG *integerglob=NULL,iglob=1,irefine=1,*inodestet=NULL,nnodestet=0,i,
    istart,j,iquadratic,nenew,nline,nset_=0,*ipoinp=NULL,nkold,nkmt,
    *inp=NULL,*ipoinpc=NULL,idummy[2]={0,0},nuel_=0,inp_size,nentries=19,
    nkold_,neold_,nkonold_,*ielprop=NULL,*iponor=NULL,*iponoel=NULL,
    *rig=NULL,*ne2boun=NULL,*ielorien=NULL,*inotr=NULL,im,mt=mi[1]+1,
    *ielmat=NULL,*irandomtype=NULL,*iparentel=NULL,*ipompc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,*nodempc=NULL,*kon=NULL,*ipkon=NULL,*iamt1=NULL,*ipobody=NULL,
    limit,*ipobody2=NULL,ifreebody2,index,index2,*iprfn=NULL,*konrfn=NULL,
    *jq=NULL,*irow=NULL,*icol=NULL,*loc=NULL,*irowt=NULL,*jqt=NULL,
    *itemp=NULL,*ixcol=NULL,*ipoface=NULL,*nodface=NULL,iamplitude,
    num_cpus,idelta,isum,num_cpus_loc,jstart,kstart,k,m,jstartn,
    kstartn;

  double *doubleglob=NULL,sigma=0.,*thickn=NULL,*thicke=NULL,*offset=NULL,
    *t0=NULL,*t0g=NULL,*t1g=NULL,*prestr=NULL,*vold=NULL,*veold=NULL,
    *randomval=NULL,*t1=NULL,*coefmpc=NULL,*co=NULL,*ratiorfn=NULL,
    *au=NULL,*t0new=NULL,*t1new=NULL,*vnew=NULL,*venew=NULL;
      
  /* variables for multithreading procedure */
  
  ITG sys_cpus,*ithread=NULL;
  char *env,*envsys;

  ielprop=*ielpropp;iponor=*iponorp;thickn=*thicknp;thicke=*thickep;
  offset=*offsetp;iponoel=*iponoelp;rig=*rigp;ne2boun=*ne2bounp;
  ielorien=*ielorienp;inotr=*inotrp;t0=*t0p;t0g=*t0gp;t1g=*t1gp;
  prestr=*prestrp;vold=*voldp;veold=*veoldp;ielmat=*ielmatp;
  irandomtype=*irandomtypep;randomval=*randomvalp;t1=*t1p;
  iparentel=*iparentelp;ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;
  ilmpc=*ilmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;co=*cop;kon=*konp;
  ipkon=*ipkonp;lakon=*lakonp;iamt1=*iamt1p;ipobody=*ipobodyp;
  iprfn=*iprfnp;konrfn=*konrfnp;ratiorfn=*ratiorfnp;
  
  /* adding the refined mesh (only in the first step) */
  
  if(*istep==1){

    num_cpus=0;
    sys_cpus=0;
  
    /* explicit user declaration prevails */
  
    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
      sys_cpus=atoi(envsys);
      if(sys_cpus<0) sys_cpus=0;
    }
  
    /* automatic detection of available number of processors */
  
    if(sys_cpus==0){
      sys_cpus=getSystemCPUs();
      if(sys_cpus<1) sys_cpus=1;
    }
  
    /* else global declaration, if any, applies */
  
    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
      if(env)
	num_cpus=atoi(env);
      if(num_cpus<1) {
	num_cpus=1;
      }else if(num_cpus>sys_cpus){
	num_cpus=sys_cpus;
      }
    }

    nkold=*nk;
  
    /* reading the input file with information on the refined mesh */

    NNEW(ipoinp,ITG,2*nentries);

    strcpy2(fnrfn,jobnamec,132);
    strcat(fnrfn,".rfn");

    readinput(fnrfn,&inpc,&nline,&nset_,ipoinp,&inp,&ipoinpc,idummy,&nuel_,
	      &inp_size); 

    /* determine new values for nk_, ne_ and nkon_ */

    nkold_=*nk_;
    neold_=*ne_;
    nkonold_=*nkon_;
    FORTRAN(allocation_rfn,(nk_,ne_,nkon_,ipoinp,ipoinpc,inpc,inp));

    /* reallocating the fields depending on nk_, ne_ and nkon_ */

    RENEW(co,double,3**nk_);
    RENEW(kon,ITG,*nkon_);
    RENEW(ipkon,ITG,*ne_);
    for(i=neold_;i<*ne_;i++) ipkon[i]=-1;
    RENEW(lakon,char,8**ne_);
    if(*nprop_>0){
      RENEW(ielprop,ITG,*ne_);
      for(i=neold_;i<*ne_;i++) ielprop[i]=-1;
    }
    if((*ne1d!=0)||(*ne2d!=0)){
      RENEW(iponor,ITG,2**nkon_);
      for(i=2*nkonold_;i<2**nkon_;i++) iponor[i]=-1;
      RENEW(thickn,double,2**nk_);
      RENEW(thicke,double,mi[2]**nkon_);
      RENEW(offset,double,2**ne_);
      RENEW(iponoel,ITG,*nk_);
      RENEW(rig,ITG,*nk_);
      RENEW(ne2boun,ITG,2**nk_);
    }

    if(*norien>0){
      RENEW(ielorien,ITG,mi[2]**ne_);
      ITGMEMSET(ielorien,mi[2]*neold_,mi[2]**ne_,0);
    }

    if(*ntrans>0){
      RENEW(inotr,ITG,2**nk_);
      ITGMEMSET(inotr,2*nkold_,2**nk_,0);
    }

    if(ithermal[0]>0){
      RENEW(t0,double,*nk_);
      RENEW(t1,double,*nk_);
      if((*ne1d!=0)||(*ne2d!=0)){
	RENEW(t0g,double,2**nk_);
	RENEW(t1g,double,2**nk_);
      }
      DMEMSET(t0,nkold_,*nk_,1.2357111319);
      DMEMSET(t1,nkold_,*nk_,1.2357111319);
      if(*nam>0){
	RENEW(iamt1,ITG,*nk_);
	ITGMEMSET(iamt1,nkold_,*nk_,0);
      }
    }

    if(*iprestr>0){
      RENEW(prestr,double,6*mi[0]**ne_);
      DMEMSET(prestr,6*mi[0]*neold_,6*mi[0]**ne_,0.);
    }
      
    RENEW(vold,double,mt**nk_);
    DMEMSET(vold,mt*nkold_,mt**nk_,0.);
    RENEW(veold,double,mt**nk_);
    DMEMSET(veold,mt*nkold_,mt**nk_,0.);
    
    RENEW(ielmat,ITG,mi[2]**ne_);
    ITGMEMSET(ielmat,mi[2]*neold_,mi[2]**ne_,0);
    
    if(irobustdesign[0]>0){
      RENEW(irandomtype,ITG,*nk_);
      RENEW(randomval,double,2**nk_);
    }

    NNEW(iparentel,ITG,*ne_);

    FORTRAN(calinput_rfn,(co,filab,set,istartset,iendset,ialset,nset,&nset_,
			  nalset,nalset_,mi,kon,ipkon,lakon,nkon,ne,ne_,
			  iponor,xnor,istep,ipoinp,inp,iaxial,ipoinpc,
			  network,nlabel,iuel,&nuel_,ielmat,inpc,iperturb,
			  iprestr,nk,nk_,ntie,tieset,iparentel,tietol));

    SFREE(inp);SFREE(inpc);SFREE(ipoinp);SFREE(ipoinpc);
    RENEW(iparentel,ITG,*ne);

    /* transferring the material and orientation information from the
       parent elements */
    
    for(i=0;i<*ne_;i++){
      if(iparentel[i]>0){
	ielmat[i]=ielmat[iparentel[i]-1];
      }
    }

    if(*norien>0){
      for(i=0;i<*ne_;i++){
	if(iparentel[i]>0){
	  ielorien[i]=ielorien[iparentel[i]-1];
	}
      }
    }
    
    /* create MPC's to determine the temperatures t0, t1, vold,
       and veold in the new mesh based on the values in the old mesh */

    strcpy2(masterurffile,jobnamec,132);
    strcat(masterurffile,".urf.frd");

    /* reading the old mesh */

    getglobalresults(masterurffile,&integerglob,&doubleglob,nboun,iamboun,
		     xboun,nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
		     ithermal,nk,t1,iamt1,&sigma,&irefine);

    /* parallellizing the calculation of the interpolation coefficients of
       the new nodes within the old mesh */
    
    pthread_t tid[num_cpus];
    
    if(num_cpus>*nk){
      num_cpus_loc=*nk;
    }else{
      num_cpus_loc=num_cpus;
    }

    /* allocating the fields for the cpu load
       cpu i covers the nodes j: nkapar[i]<=j<nkbpar[i] (C-convention) or
                                 nkapar[i]<j<=nkbpar[i] (FORTRAN-convention)  */
    
    NNEW(nkapar,ITG,num_cpus_loc);
    NNEW(nkbpar,ITG,num_cpus_loc);

    /* dividing the new nodes among the cpus */
    
    idelta=(ITG)floor(*nk)/(double)(num_cpus_loc);
    ideltamax=idelta;
    isum=0;
    for(i=0;i<num_cpus_loc;i++){
      nkapar[i]=isum;
      if(i!=num_cpus_loc-1){
	isum+=idelta;
      }else{
	isum=*nk;
	if(ideltamax<*nk-nkapar[i]) ideltamax=*nk-nkapar[i];
      }
      nkbpar[i]=isum;
    }

    /* performing the coefficient calculation in parallel */
    
    NNEW(iprfn,ITG,num_cpus_loc*(ideltamax+1));
    NNEW(konrfn,ITG,num_cpus_loc*20*ideltamax);
    NNEW(ratiorfn,double,num_cpus_loc*20*ideltamax);

    co1=co;doubleglob1=doubleglob;integerglob1=integerglob;
    iprfn1=iprfn;konrfn1=konrfn;ratiorfn1=ratiorfn;

    NNEW(ithread,ITG,num_cpus_loc);

    for(i=0; i<num_cpus_loc; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL,
		     (void *)genratiomt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus_loc; i++)  pthread_join(tid[i], NULL);
    
    /* reordering the information in serial form */

    for(i=0;i<num_cpus_loc;i++){
      jstart=i*(ideltamax+1);
      jstartn=nkapar[i];
      kstart=i*20*ideltamax;
      if(i==0){
	kstartn=0;
      }else{
	kstartn=iprfn[nkbpar[i-1]];
      }
      for(k=0;k<nkbpar[i]-nkapar[i];k++){
	iprfn[jstartn+k]=kstartn+iprfn[jstart+k];
	for(m=0;m<iprfn[jstart+k+1]-iprfn[jstart+k];m++){
	  konrfn[iprfn[jstartn+k]+m]=konrfn[kstart+iprfn[jstart+k]+m];
	  ratiorfn[iprfn[jstartn+k]+m]=ratiorfn[kstart+iprfn[jstart+k]+m];
	}
      }
      iprfn[jstartn+nkbpar[i]-nkapar[i]]=kstartn+iprfn[jstart+nkbpar[i]-nkapar[i]];
    }

    RENEW(iprfn,ITG,*nk+1);
    RENEW(konrfn,ITG,iprfn[*nk]);
    RENEW(ratiorfn,double,iprfn[*nk]);

    /* interpolating the initial temperatures t0 */

    if(ithermal[0]>0){
      NNEW(t0new,double,*nk);
      for(i=0;i<*nk;i++){
	for(j=0;j<iprfn[i+1]-iprfn[i];j++){
	  t0new[i]+=ratiorfn[iprfn[i]+j]*t0[konrfn[iprfn[i]+j]-1];
	}
      }
      cpypardou(t0,t0new,nk,&num_cpus_loc);
      SFREE(t0new);
      if(*nam>0){
	for(i=0;i<*nk;i++){
	  iamplitude=iamt1[konrfn[iprfn[i]]-1];
	  for(j=1;j<iprfn[i+1]-iprfn[i];j++){
	    if(iamt1[konrfn[iprfn[i]+j]-1]!=iamplitude){
	      iamplitude=-1;
	      break;
	    }
	  }
	  if(iamplitude==-1){
	    printf(" *ERROR in readnewmesh: temperature amplitude\n");
	    printf("        is interpolated in refined node %" ITGFORMAT " \n",i);
	    printf("        but amplitude in the surrounding nodes is\n");
	    printf("        not the same.\n");
	    FORTRAN(stop,());
	  }
	  iamt1[i]=iamplitude;
	}
      }
    }

    /* interpolating vold and veold */

    nkmt=*nk*mt;
    NNEW(vnew,double,nkmt);
    NNEW(venew,double,nkmt);
    for(i=0;i<*nk;i++){
      for(j=0;j<iprfn[i+1]-iprfn[i];j++){
	for(k=0;k<mt;k++){
	  vnew[i*mt+k]+=ratiorfn[iprfn[i]+j]*vold[(konrfn[iprfn[i]+j]-1)*mt+k];
	  venew[i*mt+k]+=ratiorfn[iprfn[i]+j]*veold[(konrfn[iprfn[i]+j]-1)*mt+k];
	}
      }
    }
    cpypardou(vold,vnew,&nkmt,&num_cpus_loc);
    cpypardou(veold,venew,&nkmt,&num_cpus_loc);
    SFREE(vnew);
    SFREE(venew);
    
  }

  /* the remaining lines are executed in each step */

  /* transfering the body load information from the
     parent elements */

   if(*nbody>0){
    limit=(ITG)(1.1**ne);
    if(limit<100) limit=100;
    NNEW(ipobody2,ITG,2*limit);
    ifreebody2=*ne+1;
  
    for(i=0;i<*ne;i++){
      if(ipkon[i]<0) continue;
      
      if(iparentel[i]==0){
	index=i+1;
      }else{
	index=iparentel[i];
      }
      
      index2=i+1;
      ipobody2[2*index2-2]=ipobody[2*index-2];
      index=ipobody[2*index-1];

      do{
	if(index!=0){
	  ipobody2[2*index2-1]=ifreebody2;
	  index2=ifreebody2;
	  
	  ifreebody2++;
	  if(ifreebody2>limit){
	    limit=(ITG)(1.1*limit);
	    RENEW(ipobody2,ITG,2*limit);
	  }

	  ipobody2[2*index2-2]=ipobody[2*index-2];
	  index=ipobody[2*index-1];
	}else{
	  ipobody2[2*index2-1]=0;
	  break;
	}
      }while(1);
    }
  
    RENEW(ipobody,ITG,2*(ifreebody2-1));
    memcpy(ipobody,ipobody2,sizeof(ITG)*2*(ifreebody2-1));
    *ifreebody=ifreebody2;
    SFREE(ipobody2);
    }

  /* remove the refine request, if any */
  
  if(strcmp1(&filab[4089],"RM")==0){
    strcpy1(&filab[4089],"  ",2);
  }

  if(ithermal[0]>0){

    /* interpolating the initial temperatures t1 */

    NNEW(t1new,double,*nk);
    for(i=0;i<*nk;i++){
      for(j=0;j<iprfn[i+1]-iprfn[i];j++){
	t1new[i]+=ratiorfn[iprfn[i]+j]*t1[konrfn[iprfn[i]+j]-1];
      }
    }
    cpypardou(t1,t1new,nk,&num_cpus_loc);
    SFREE(t1new);

  }

  *ielpropp=ielprop;*iponorp=iponor;*thicknp=thickn;*thickep=thicke;
  *offsetp=offset;*iponoelp=iponoel;*rigp=rig;*ne2bounp=ne2boun;
  *ielorienp=ielorien;*inotrp=inotr;*t0p=t0;*t0gp=t0g;*t1gp=t1g;
  *prestrp=prestr;*voldp=vold;*veoldp=veold;*ielmatp=ielmat;
  *irandomtypep=irandomtype;*randomvalp=randomval;*t1p=t1;
  *iparentelp=iparentel;*ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;
  *ilmpcp=ilmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;*cop=co;*konp=kon;
  *ipkonp=ipkon;*lakonp=lakon;*iamt1p=iamt1;*ipobodyp=ipobody;
  *iprfnp=iprfn;*konrfnp=konrfn;*ratiorfnp=ratiorfn;
  
  return NULL;

}

/* subroutine for multithreading of calcenergy */

void *genratiomt(ITG *i){

  ITG nka,nkb,index1,index2;

  nka=nkapar[*i]+1;
  nkb=nkbpar[*i];

  index1=*i*(ideltamax+1);
  index2=*i*20*ideltamax;

  FORTRAN(genratio,(co1,doubleglob1,integerglob1,&nka,&nkb,
		    &iprfn1[index1],&konrfn1[index2],&ratiorfn1[index2]));

  return NULL;
}
