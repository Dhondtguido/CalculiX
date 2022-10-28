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
#include "mortar.h"

void dealloc_cal(ITG *ncs_,ITG **icsp,ITG *mcs,double **csp,
		 char **tiesetp,double **tietolp,double **cop,
		 ITG **konp,ITG **ipkonp,char **lakonp,ITG **nodebounp,
		 ITG **ndirbounp,char **typebounp,double **xbounp,
		 ITG **ikbounp,ITG **ilbounp,ITG **nodebounoldp,
		 ITG **ndirbounoldp,double **xbounoldp,ITG **ipompcp,
		 char **labmpcp,ITG **ikmpcp,ITG **ilmpcp,
		 double **fmpcp,ITG **nodempcp,double **coefmpcp,
		 ITG **nodempcrefp,double **coefmpcrefp,ITG **ikmpcrefp,
		 ITG **nodeforcp,ITG **ndirforcp,double **xforcp,
		 ITG **ikforcp,ITG **ilforcp,double **xforcoldp,
		 ITG **nelemloadp,char **sideloadp,double **xloadp,
		 double **xloadoldp,char **cbodyp,ITG **ibodyp,
		 double **xbodyp,double **xbodyoldp,ITG *nam,
		 ITG **iambounp,ITG **iamforcp,ITG **iamloadp,
		 char **amnamep,double **amtap,ITG **namtap,char **setp,
		 ITG **istartsetp,ITG **iendsetp,ITG **ialsetp,
		 double **elconp,ITG **nelconp,double **rhconp,
		 ITG **nrhconp,double **shconp,ITG **nshconp,
		 double **coconp,ITG **ncoconp,double **alconp,
		 ITG **nalconp,double **alzerop,ITG *nprop,
		 ITG **ielpropp,double **propp,ITG *npmat_,
		 double **pliconp,ITG **npliconp,double **plkconp,
		 ITG **nplkconp,ITG *ndamp,double **daconp,ITG *norien,
		 char **ornamep,double **orabp,ITG **ielorienp,
		 ITG *ntrans,double **trabp,ITG **inotrp,ITG *iprestr,
		 double **prestrp,ITG *ithermal,double **t0p,
		 double **t1p,double **t1oldp,ITG **iamt1p,ITG *ne1d,
		 ITG *ne2d,double **t0gp,double **t1gp,
		 ITG *irobustdesign,ITG **irandomtypep,
		 double **randomvalp,char **prlabp,char **prsetp,
		 char **filabp,double **xmodalp,ITG **ielmatp,
		 char **matnamep,double **stip,double **emep,
		 double **enerp,double **xstatep,double **voldp,
		 double **veoldp,double **velp,double **velop,
		 double **veloop,ITG **iponorp,double **xnorp,
		 ITG **knorp,double **thickep,double **offsetp,
		 ITG **iponoelp,ITG **inoelp,ITG **rigp,
		 ITG **ne2bounp,ITG **islavsurfp,ITG *mortar,
		 double **pslavsurfp,double **clearinip,
		 ITG *nobject_,char **objectsetp,ITG *nmethod,ITG *iperturb,
		 ITG *irefineloop,ITG **iparentelp,ITG **iprfnp,
		 ITG **konrfnp,double **ratiorfnp,char **headingp,
		 ITG **nodedesip,double **dgdxglobp,double **g0p,
		 ITG *nuel_,double **xdesip,ITG *nfc,double **coeffcp,
		 ITG **ikdcp,double **edcp,double **coinip){

  char *tieset=NULL,*lakon=NULL,*typeboun=NULL,*labmpc=NULL,*sideload=NULL,
    *cbody=NULL,*amname=NULL,*set=NULL,*orname=NULL,*prlab=NULL,*prset=NULL,
    *filab=NULL,*matname=NULL,*objectset=NULL,*heading=NULL;
  
  ITG *ics=NULL,*kon=NULL,*ipkon=NULL,*nodeboun=NULL,*ndirboun=NULL,
    *ikboun=NULL,*ilboun=NULL,*nodebounold=NULL,*ndirbounold=NULL,
    *ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,*nodempc=NULL,*nodempcref=NULL,
    *ikmpcref=NULL,*nodeforc=NULL,*ndirforc=NULL,*ikforc=NULL,*ilforc=NULL,
    *nelemload=NULL,*ibody=NULL,*iamboun=NULL,*iamforc=NULL,*iamload=NULL,
    *namta=NULL,*istartset=NULL,*iendset=NULL,*ialset=NULL,*nelcon=NULL,
    *nrhcon=NULL,*nshcon=NULL,*ncocon=NULL,*nalcon=NULL,*ielprop=NULL,
    *nplicon=NULL,*nplkcon=NULL,*ielorien=NULL,*inotr=NULL,*iamt1=NULL,
    *irandomtype=NULL,*ielmat=NULL,*iponor=NULL,*knor=NULL,*iponoel=NULL,
    *inoel=NULL,*ne2boun=NULL,*islavsurf=NULL,*rig=NULL,*iparentel=NULL,
    *iprfn=NULL,*konrfn=NULL,*nodedesi=NULL,*ikdc=NULL;

  double *cs=NULL,*tietol=NULL,*co=NULL,*xboun=NULL,*xbounold=NULL,
    *fmpc=NULL,*coefmpc=NULL,*coefmpcref=NULL,*xforc=NULL,*xforcold=NULL,
    *xload=NULL,*xloadold=NULL,*xbody=NULL,*xbodyold=NULL,*amta=NULL,
    *elcon=NULL,*rhcon=NULL,*shcon=NULL,*cocon=NULL,*alcon=NULL,*alzero=NULL,
    *prop=NULL,*plicon=NULL,*plkcon=NULL,*dacon=NULL,*orab=NULL,*trab=NULL,
    *prestr=NULL,*t0=NULL,*t1=NULL,*t1old=NULL,*t0g=NULL,*t1g=NULL,
    *randomval=NULL,*xmodal=NULL,*sti=NULL,*eme=NULL,*ener=NULL,*xstate=NULL,
    *vold=NULL,*veold=NULL,*vel=NULL,*velo=NULL,*veloo=NULL,*xnor=NULL,
    *thicke=NULL,*offset=NULL,*pslavsurf=NULL,*clearini=NULL,*ratiorfn=NULL,
    *dgdxglob=NULL,*g0=NULL,*xdesi=NULL,*coeffc=NULL,*edc=NULL,*coini=NULL;

  ics=*icsp;cs=*csp;tieset=*tiesetp;tietol=*tietolp;co=*cop;kon=*konp;
  ipkon=*ipkonp;lakon=*lakonp;nodeboun=*nodebounp;ndirboun=*ndirbounp;
  typeboun=*typebounp;xboun=*xbounp;ikboun=*ikbounp;ilboun=*ilbounp;
  nodebounold=*nodebounoldp;ndirbounold=*ndirbounoldp;xbounold=*xbounoldp;
  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;fmpc=*fmpcp;
  nodempc=*nodempcp;coefmpc=*coefmpcp;nodeforc=*nodeforcp;ndirforc=*ndirforcp;
  xforc=*xforcp;ikforc=*ikforcp;ilforc=*ilforcp;xforcold=*xforcoldp;
  nelemload=*nelemloadp;sideload=*sideloadp;xload=*xloadp;xloadold=*xloadoldp;
  cbody=*cbodyp;ibody=*ibodyp;xbody=*xbodyp;xbodyold=*xbodyoldp;
  iamboun=*iambounp;iamforc=*iamforcp;iamload=*iamloadp;amname=*amnamep;
  amta=*amtap;namta=*namtap;set=*setp;istartset=*istartsetp;iendset=*iendsetp;
  ialset=*ialsetp;elcon=*elconp;nelcon=*nelconp;rhcon=*rhconp;nrhcon=*nrhconp;
  shcon=*shconp;nshcon=*nshconp;cocon=*coconp;ncocon=*ncoconp;alcon=*alconp;
  nalcon=*nalconp;alzero=*alzerop;ielprop=*ielpropp;prop=*propp;
  plicon=*pliconp;nplicon=*npliconp;plkcon=*plkconp;nplkcon=*nplkconp;
  dacon=*daconp;orname=*ornamep;orab=*orabp;ielorien=*ielorienp;trab=*trabp;
  inotr=*inotrp;prestr=*prestrp;t0=*t0p;t1=*t1p;t1old=*t1oldp;iamt1=*iamt1p;
  t0g=*t0gp;t1g=*t1gp;irandomtype=*irandomtypep;randomval=*randomvalp;
  prlab=*prlabp;prset=*prsetp;filab=*filabp;xmodal=*xmodalp;ielmat=*ielmatp;
  matname=*matnamep;sti=*stip;eme=*emep;ener=*enerp;xstate=*xstatep;
  vold=*voldp;veold=*veoldp;vel=*velp;velo=*velop;veloo=*veloop;iponor=*iponorp;
  xnor=*xnorp;knor=*knorp;thicke=*thickep;offset=*offsetp;iponoel=*iponoelp;
  inoel=*inoelp;rig=*rigp;ne2boun=*ne2bounp;islavsurf=*islavsurfp;
  pslavsurf=*pslavsurfp;clearini=*clearinip;objectset=*objectsetp;
  iparentel=*iparentelp;coefmpcref=*coefmpcrefp;nodempcref=*nodempcrefp;
  ikmpcref=*ikmpcrefp;iprfn=*iprfnp;konrfn=*konrfnp;ratiorfn=*ratiorfnp;
  heading=*headingp;nodedesi=*nodedesip;g0=*g0p;dgdxglob=*dgdxglobp;
  xdesi=*xdesip;coeffc=*coeffcp;ikdc=*ikdcp;edc=*edcp;coini=*coinip;
								 
  /* deallocating all fields except the *inp fields */
			 
  if(*ncs_>0) SFREE(ics);
  if(*mcs>0) SFREE(cs);
  SFREE(tieset);SFREE(tietol);

  SFREE(co);SFREE(kon);SFREE(ipkon);SFREE(lakon);

  SFREE(nodeboun);SFREE(ndirboun);SFREE(typeboun);SFREE(xboun);SFREE(ikboun);
  SFREE(ilboun);SFREE(nodebounold);SFREE(ndirbounold);SFREE(xbounold);

  SFREE(ipompc);SFREE(labmpc);SFREE(ikmpc);SFREE(ilmpc);SFREE(fmpc);
  SFREE(nodempc);SFREE(coefmpc);

  if(*nfc>0){
    SFREE(coeffc);SFREE(ikdc);SFREE(edc);
  }

  SFREE(nodempcref);SFREE(coefmpcref);SFREE(ikmpcref);

  SFREE(nodeforc);SFREE(ndirforc);SFREE(xforc);SFREE(ikforc);SFREE(ilforc);
  SFREE(xforcold);

  SFREE(nelemload);SFREE(sideload);SFREE(xload);SFREE(xloadold);

  SFREE(cbody);SFREE(ibody);SFREE(xbody);SFREE(xbodyold);

  if(*nam>0){SFREE(iamboun);SFREE(iamforc);SFREE(iamload);SFREE(amname);
    SFREE(amta);SFREE(namta);}

  SFREE(set);SFREE(istartset);SFREE(iendset);SFREE(ialset);

  SFREE(elcon);SFREE(nelcon);SFREE(rhcon);SFREE(nrhcon);SFREE(shcon);
  SFREE(nshcon);
  SFREE(cocon);SFREE(ncocon);SFREE(alcon);SFREE(nalcon);SFREE(alzero);
  if(*nprop>0){SFREE(ielprop);SFREE(prop);}
  if(*npmat_>0){SFREE(plicon);SFREE(nplicon);SFREE(plkcon);SFREE(nplkcon);}
  if(*ndamp>0){SFREE(dacon);}

  if(*norien>0){SFREE(orname);SFREE(orab);SFREE(ielorien);}
  if(*ntrans>0){SFREE(trab);SFREE(inotr);}
  if(*iprestr>0){SFREE(prestr);}

  if(ithermal[0]!=0){
    SFREE(t0);SFREE(t1);SFREE(t1old);
    if(*nam>0) SFREE(iamt1);
    if((*ne1d!=0)||(*ne2d!=0)||(*nuel_!=0)){SFREE(t0g);SFREE(t1g);}
  }

  if(irobustdesign[0]>0){SFREE(irandomtype);SFREE(randomval);}
  
  SFREE(prlab);SFREE(prset);SFREE(filab);SFREE(xmodal);

  SFREE(ielmat);SFREE(matname);

  SFREE(sti);SFREE(eme);SFREE(ener);SFREE(xstate);

  SFREE(vold);
  if(*irefineloop==1){
    if((*nmethod==4)||(*nmethod==5)||(*nmethod==8)||(*nmethod==9)||
       ((abs(*nmethod)==1)&&(iperturb[0]>=2))){SFREE(veold);}
  }else{
    SFREE(veold);
  }
  SFREE(vel);
  SFREE(velo);
  SFREE(veloo);

  if((*ne1d!=0)||(*ne2d!=0)){
    SFREE(iponor);SFREE(xnor);SFREE(knor);SFREE(thicke);SFREE(offset);
    SFREE(iponoel);SFREE(inoel);SFREE(rig);SFREE(ne2boun);
  }

  SFREE(islavsurf);
  if(*mortar==1){SFREE(pslavsurf);SFREE(clearini);}

  if(*nobject_>0){SFREE(objectset);SFREE(nodedesi);SFREE(dgdxglob);
    SFREE(g0);SFREE(xdesi);SFREE(coini);}

  if(*irefineloop>0){SFREE(iparentel);SFREE(iprfn);SFREE(konrfn);
    SFREE(ratiorfn);}

  SFREE(heading);

  *icsp=ics;*csp=cs;*tiesetp=tieset;*tietolp=tietol;*cop=co;*konp=kon;
  *ipkonp=ipkon;*lakonp=lakon;*nodebounp=nodeboun;*ndirbounp=ndirboun;
  *typebounp=typeboun;*xbounp=xboun;*ikbounp=ikboun;*ilbounp=ilboun;
  *nodebounoldp=nodebounold;*ndirbounoldp=ndirbounold;*xbounoldp=xbounold;
  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;*fmpcp=fmpc;
  *nodempcp=nodempc;*coefmpcp=coefmpc;*nodeforcp=nodeforc;*ndirforcp=ndirforc;
  *xforcp=xforc;*ikforcp=ikforc;*ilforcp=ilforc;*xforcoldp=xforcold;
  *nelemloadp=nelemload;*sideloadp=sideload;*xloadp=xload;*xloadoldp=xloadold;
  *cbodyp=cbody;*ibodyp=ibody;*xbodyp=xbody;*xbodyoldp=xbodyold;
  *iambounp=iamboun;*iamforcp=iamforc;*iamloadp=iamload;*amnamep=amname;
  *amtap=amta;*namtap=namta;*setp=set;*istartsetp=istartset;*iendsetp=iendset;
  *ialsetp=ialset;*elconp=elcon;*nelconp=nelcon;*rhconp=rhcon;*nrhconp=nrhcon;
  *shconp=shcon;*nshconp=nshcon;*coconp=cocon;*ncoconp=ncocon;*alconp=alcon;
  *nalconp=nalcon;*alzerop=alzero;*ielpropp=ielprop;*propp=prop;
  *pliconp=plicon;*npliconp=nplicon;*plkconp=plkcon;*nplkconp=nplkcon;
  *daconp=dacon;*ornamep=orname;*orabp=orab;*ielorienp=ielorien;*trabp=trab;
  *inotrp=inotr;*prestrp=prestr;*t0p=t0;*t1p=t1;*t1oldp=t1old;*iamt1p=iamt1;
  *t0gp=t0g;*t1gp=t1g;*irandomtypep=irandomtype;*randomvalp=randomval;
  *prlabp=prlab;*prsetp=prset;*filabp=filab;*xmodalp=xmodal;*ielmatp=ielmat;
  *matnamep=matname;*stip=sti;*emep=eme;*enerp=ener;*xstatep=xstate;
  *voldp=vold;*veoldp=veold;*velp=vel;*velop=velo;*veloop=veloo;*iponorp=iponor;
  *xnorp=xnor;*knorp=knor;*thickep=thicke;*offsetp=offset;*iponoelp=iponoel;
  *inoelp=inoel;*rigp=rig;*ne2bounp=ne2boun;*islavsurfp=islavsurf;
  *pslavsurfp=pslavsurf;*clearinip=clearini;*objectsetp=objectset;
  *iparentelp=iparentel;*coefmpcrefp=coefmpcref;*nodempcrefp=nodempcref;
  *ikmpcrefp=ikmpcref;*iprfnp=iprfn;*konrfnp=konrfn;*ratiorfnp=ratiorfn;
  *headingp=heading;*nodedesip=nodedesi;*dgdxglobp=dgdxglob;*g0p=g0;
  *xdesip=xdesi;*coeffcp=coeffc;*ikdcp=ikdc;*edcp=edc;*coinip=coini;
  
  return;
}
