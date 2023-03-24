/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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

/**
 * function called before solver transforming the SPCs/MPCs, the matrix and the right hand side for quadratic elements 
 (see Phd-thesis Sitzmann,Chapter 4)

 * * Author: Saskia Sitzmann
 *
 *  [out] iflagact	flag indicating, whether the coupling matrices have to 
 be recalculated
 *  [out] nzsc2		number of nonzero,nondiagonal entries of intermediate 
 system matrix
 *  [out] auc2p		intermediate system matrix
 *  [out] adc2p         intermediate system matrix, diagonal terms
 *  [out] irowc2p       rows for intermediate system matrix
 *  [out] icolc2p	columns for intermediate system matrix
 *  [out] jqc2p		pointer to irowc 
 *  [out] aubdp		coupling matrix \f$ B_d[nactdof(i,p),nactdof(j,q)]\f$ 
 for all active degrees od freedoms 
 *  [out] irowbdp	field containing row numbers of aubd
 *  [out] jqbdp		pointer into field irowbd
 *  [out] aubdtilp	matrix \f$ \tilde{D}^{-1}\tilde{B}_d[nactdof(i,p),
 nactdof(j,q)]\f$ for all active degrees od freedoms
 *  [out] irowbdtilp	field containing row numbers of aubd
 *  [out] jqbdtilp	pointer into field irowbdtil 
 *  [out] aubdtil2p	coupling matrix \f$ \tilde{D}$ and 
 $\tilde{B}^2_d[nactdof(i,p),nactdof(j,q)]\f$ for all 
 active degrees od freedoms 
 *  [out] irowbdtil2p	field containing row numbers of aubdtil2
 *  [out] jqbdtil2p	pointer into field irowbdtil2
 *  [out] auddp		coupling matrix \f$ D_d[nactdof(i,p),nactdof(j,q)]\f$ 
 for all active degrees od freedoms
 *  [out] irowddp	field containing row numbers of audd
 *  [out] jqddp		pointer into field irowdd
 *  [out] auddtilp	coupling matrix \f$ \tilde{D}_d[nactdof(i,p),nactdof(j,
 q)]\f$ for all active degrees od freedoms
 *  [out] irowddtilp	field containing row numbers of audd
 *  [out] jqddtilp	pointer into field irowdd
 *  [out] auddtil2p	matrix \f$ Id_d[nactdof(i,p),nactdof(j,q)]\f$ for all 
 active degrees od freedoms
 *  [out] irowddtil2p	field containing row numbers of audd
 *  [out] jqddtil2p	pointer into field irowdd
 *  [out] auddinvp	coupling matrix \f$ \tilde{D}^{-1}_d[nactdof(i,p),
 nactdof(j,q)]\f$ for all active degrees od freedoms
 *  [out] irowddinvp	field containing row numbers of auddinv
 *  [out] jqddinvp	pointer into field irowddinv
 *  [out] jqtempp	field storing the untransformed stiffness matrix 
 representation
 *  [out] irowtempp	field storing the untransformed stiffness matrix 
 representation
 *  [out] icoltempp	field storing the untransformed stiffness matrix 
 representation
 *  [out] nzstemp	field storing the untransformed stiffness matrix size
 *  [out] slavnor	slave normal
 *  [out] slavtan	slave tangent 
 *  [out] nboun2        number of transformed SPCs
 *  [out] ndirboun2p	(i) direction of transformed SPC i 
 *  [out] nodeboun2p    (i) node of transformed SPC i
 *  [out] xboun2p       (i) value of transformed SPC i
 *  [out] nmpc2		number of transformed mpcs
 *  [out] ipompc2p      (i) pointer to nodempc and coeffmpc for transformed 
 MPC i
 *  [out] nodempc2p     i and directions of transformed MPCs
 *  [out] coefmpc2p     coefficients of transformed MPCs
 *  [out] labmpc2p 	transformed mpc labels
 *  [out] ikboun2p      sorted dofs idof=8*(node-1)+dir for transformed SPCs
 *  [out] ilboun2p      transformed SPC numbers for sorted dofs
 *  [out] ikmpc2p 	sorted dofs idof=8*(node-1)+dir for transformed MPCs
 *  [out] ilmpc2p	transformed SPC numbers for sorted dofs
 *  [out] nslavspcp	(2*i) pointer to islavspc...
 *  [out] islavspcp     ... which stores SPCs for slave node i
 *  [out] nsspc         number of SPC for slave i
 *  [out] nslavmpcp	(2*i) pointer to islavmpc...
 *  [out] islavmpcp	... which stores MPCs for slave node i
 *  [out] nsmpc		number of MPC for slave i
 *  [out] nslavspc2p	(2*i) pointer to islavspc2...
 *  [out] islavspc2p    ... which stores transformed SPCs for slave node i
 *  [out] nsspc2        number of transformed SPC for slave i
 *  [out] nslavmpc2p	(2*i) pointer to islavmpc2...
 *  [out] islavmpc2p	... which stores transformed MPCs for slave node i
 *  [out] nsmpc2	number of transformed MPC for slave i 
 *  [in] imastnode	field storing the i of the master surfaces
 *  [in] nmastnode	(i)pointer into field imastnode for contact tie i 
 *  [out] nmastspcp	(2*i) pointer to imastspc...
 *  [out] imastspcp     ... which stores SPCs for master node i
 *  [out] nmspc         number of SPC for master i
 *  [out] nmastmpcp	(2*i) pointer to imastmpc...
 *  [out] imastmpcp	... which stores MPCs for master node i
 *  [out] nmmpc		number of MPC for master i 
 *  [in] islavelinv       (i)==0 if there is no slave node in the element, 
 >0 otherwise
 *  [in] islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration
 points generated before looking at face i 
 *  [in] autloc		transformation matrix \f$ T[p,q]\f$ for slave i \f$ p,q
 \f$  
 *  [in] irowtloc		field containing row numbers of autloc
 *  [in] jqtloc	        pointer into field irowtloc
 *  [in] autlocinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave i 
 \f$ p,q \f$ 
 *  [in] irowtlocinv	field containing row numbers of autlocinv
 *  [in] jqtlocinv	pointer into field irowtlocinv
 *  [in] nk2		number or generated points needed for transformed SPCs 
 *  [in]  iflagdualquad   flag indicating what mortar contact is used 
 (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG 
 quad-quad) 
 *  [out] f_cmp		not used any more
 *  [out] f_csp		contact force for active degrees of freedom
 *  [in]	auxtil2		auxilary field
 *  [in] pslavsurf	field storing  position xil, etal and weight for 
 integration point on slave side
 *  [in] pmastsurf 	field storing position and etal for integration points 
 on master side
 *  [in] islavact	(i) indicates, if slave node i is active 
 (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, 
 =0 inactive node, =1 sticky node, =2 slipping/active 
 node) 
 *  [in] cdn		?
 *  [in] cvtilini		C*v at start of the increment
 *  [in] cvtil		C*v
 *  [in] idamping		flag indicating whether damping is used
 *  [in] ilin		flag indicating wheter fist iteration is calculated 
 linear geometrically 
 *  [in] iperturb_sav	saved iperturb values 
 *  [out] nodeforc2p	transformed point force, node
 *  [out] ndirforc2p	transformed point force, dir
 *  [out] xforc2p		transformed point force, value
 *  [out] nforc2		number of transformed point forces   
 **/

void premortar(ITG *iflagact,ITG *ismallsliding,ITG *nzs,ITG *nzsc2,
	       double **auc2p,double **adc2p,ITG **irowc2p,ITG **icolc2p,
	       ITG **jqc2p,
	       double **aubdp,ITG **irowbdp,ITG **jqbdp,
	       double **aubdtilp,ITG **irowbdtilp,ITG **jqbdtilp,
	       double **aubdtil2p,ITG **irowbdtil2p,ITG **jqbdtil2p,
	       double **auddp,ITG **irowddp,ITG **jqddp,
	       double **auddtilp,ITG **irowddtilp,ITG **jqddtilp,
	       double **auddtil2p,ITG **irowddtil2p,ITG **jqddtil2p,
	       double **auddinvp,ITG **irowddinvp,ITG **jqddinvp,
	       ITG **jqtempp,ITG **irowtempp,ITG **icoltempp,ITG *nzstemp,
	       ITG *iit,double *slavnor,double *slavtan,
	       ITG *icol,ITG *irow,ITG *jq,
	       ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       ITG *nboun2,ITG **ndirboun2p,ITG **nodeboun2p,
	       double **xboun2p,
	       ITG *nmpc2,ITG **ipompc2p,ITG **nodempc2p,double **coefmpc2p,
	       char **labmpc2p,
	       ITG **ikboun2p,ITG **ilboun2p,ITG **ikmpc2p,ITG **ilmpc2p,
	       ITG **nslavspcp,ITG **islavspcp,ITG **nslavmpcp,
	       ITG **islavmpcp,
	       ITG **nslavspc2p,ITG **islavspc2p,ITG **nslavmpc2p,
	       ITG **islavmpc2p,
	       ITG **nmastspcp,ITG **imastspcp,ITG **nmastmpcp,
	       ITG **imastmpcp,
	       ITG **nmastmpc2p,ITG **imastmpc2p,ITG *nmmpc2,
	       ITG *nsspc,ITG *nsspc2,ITG *nsmpc,ITG *nsmpc2,
	       ITG *imastnode,ITG *nmastnode,ITG *nmspc,ITG *nmmpc,
	       double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,double *stn,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,double *prestr,
	       ITG *iprestr,char *filab,double *eme,double *emn,
	       double *een,ITG *iperturb,double *f,ITG *nactdof,
	       ITG *iout,double *qa,
	       double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
	       double *xbounact,double *xboun,ITG *nboun,ITG *ipompc,
	       ITG *nodempc,
	       double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
	       ITG *neq,double *veold,double *accold,
	       double *dtime,double *time,
	       double *ttime,double *plicon,
	       ITG *nplicon,double *plkcon,ITG *nplkcon,
	       double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
	       char *matname,ITG *mi,ITG *ielas,
	       ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
	       double *vini,double *ener,
	       double *enern,double *emeini,double *xstaten,double *eei,
	       double *enerini,double *cocon,ITG *ncocon,char *set,
	       ITG *nset,ITG *istartset,
	       ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
	       char *prset,double *qfx,double *qfn,double *trab,
	       ITG *inotr,ITG *ntrans,ITG *nelemload,
	       ITG *nload,ITG *istep,ITG *iinc,
	       double *springarea,double *reltime,ITG *ne0,double *xforc,
	       ITG *nforc,double *thicke,
	       double *shcon,ITG *nshcon,char *sideload,double *xload,
	       double *xloadold,ITG *icfd,ITG *inomat,
	       ITG *islavelinv,ITG *islavsurf,
	       ITG *iponoels,ITG *inoels,                      
	       ITG *mortar,ITG *nslavnode,ITG *islavnode,ITG *nslavs,
	       ITG *ntie,
	       double *autloc,ITG *irowtloc,ITG *jqtloc,
	       double *autlocinv,ITG *irowtlocinv,ITG *jqtlocinv,
	       ITG *nk2,ITG *iflagdualquad,
	       char *tieset,ITG *itiefac  ,ITG *rhsi,
	       double *au,double *ad,double **f_cmp,double **f_csp,
	       double *t1act,double *cam,double *bet,double *gam,
	       double *epn,
	       double *xloadact,ITG *nodeforc,ITG *ndirforc,double *xforcact,
	       double *xbodyact,ITG *ipobody,ITG *nbody,double *cgr,
	       ITG *nzl,double *sti,ITG *iexpl,ITG *mass,ITG *buckling,
	       ITG *stiffness,
	       ITG *intscheme,double *physcon,ITG *coriolis,ITG *ibody,
	       ITG *integerglob,double *doubleglob,ITG *nasym,
	       double *alpham,double *betam,double *auxtil2,
	       double *pslavsurf,double *pmastsurf,
	       double *clearini,ITG *ielprop,double *prop,
	       ITG *islavact,double *cdn,ITG *memmpc_,
	       double *cvinitil,double *cvtil,ITG *idamping,
	       ITG *ilin,ITG *iperturb_sav,double *adb,double *aub,
	       ITG **nodeforc2p,ITG **ndirforc2p,double **xforc2p,
	       ITG *nforc2,
	       ITG *itietri,double *cg,double *straight,ITG *koncont,
	       double *energyini,
	       double *energy,ITG *kscale,ITG *iponoel,ITG *inoel,ITG *nener,
	       char *orname,ITG *network,
	       char *typeboun,ITG *num_cpus,double *t0g,double *t1g,
	       double *smscale,ITG *mscalmethod,char *jobnamef){
  
  ITG im,i,ii,j,jj,k,l,mt=mi[1]+1,node,jfaces,nelems,ifaces,nope,nopes,idummy,
    *ndirboun2=NULL,*nodeboun2=NULL,*ipompc2=NULL,*nodempc2=NULL,nodes[8],
    konl[20],jj2,ifac,debug,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,
    *nslavspc=NULL,*islavspc=NULL,*nslavmpc=NULL,*islavmpc=NULL,
    *nslavspc2=NULL,*islavspc2=NULL,*nslavmpc2=NULL,*islavmpc2=NULL,
    *nmastspc=NULL,*imastspc=NULL,*nmastmpc=NULL,*imastmpc2=NULL,
    *nmastmpc2=NULL,*imastmpc=NULL,
    *irowc2=NULL,*icolc2=NULL,*jqc2=NULL,*irowbd=NULL,*jqbd=NULL,
    *irowbdtil=NULL,*jqbdtil=NULL,
    *irowbdtil2=NULL,*jqbdtil2=NULL,*irowdd=NULL,*jqdd=NULL,*irowddtil=NULL,
    *jqddtil=NULL,
    *irowddinv=NULL,*jqddinv=NULL,*irowddtil2=NULL,*jqddtil2=NULL,
    *irowtemp=NULL,*icoltemp=NULL,*jqtemp=NULL,*inum=NULL,*icoltil=NULL,
    *irowtil=NULL,*jqtil=NULL,
    *nodeforc2=NULL,*ndirforc2=NULL,mortartrafoflag=1;
  
  double alpha,*xboun2=NULL,*coefmpc2=NULL,*auc2=NULL,*adc2=NULL,*aubd=NULL,
    *aubdtil=NULL,*aubdtil2=NULL,*ftil=NULL,*fexttil=NULL,
    *audd=NULL,*auddtil=NULL,*auddinv=NULL,*auddtil2=NULL,*v=NULL,
    *stx=NULL,*fn=NULL,*autil=NULL,*finitil=NULL,*fextinitil=NULL,
    *adtil=NULL,*fmpc2=NULL,*voldtil=NULL,*vinitil=NULL,*accoldtil=NULL,
    *btil=NULL,*aubtil=NULL,*adbtil=NULL,
    *adctil=NULL,*auctil=NULL,*veoldtil=NULL,*volddummy=NULL,
    *vectornull=NULL,*f_cs=NULL,*f_cm=NULL,
    *xforc2=NULL,*fnext=NULL,*veolddummy=NULL,*accolddummy=NULL;
  
  char *labmpc2=NULL;

  debug=0;
    
  alpha=1-2*sqrt(*bet);
  
  ndirboun2=*ndirboun2p;nodeboun2=*nodeboun2p;xboun2=*xboun2p;
  ipompc2=*ipompc2p;nodempc2=*nodempc2p;coefmpc2=*coefmpc2p;
  labmpc2=*labmpc2p;ikboun2=*ikboun2p;ilboun2=*ilboun2p;ikmpc2=*ikmpc2p;
  ilmpc2=*ilmpc2p;
  nslavspc=*nslavspcp;islavspc=*islavspcp;nslavmpc=*nslavmpcp;
  islavmpc=*islavmpcp;
  nslavspc2=*nslavspc2p;islavspc2=*islavspc2p;nslavmpc2=*nslavmpc2p;
  islavmpc2=*islavmpc2p;
  nmastspc=*nmastspcp;imastspc=*imastspcp;nmastmpc=*nmastmpcp;
  imastmpc=*imastmpcp;nmastmpc2=*nmastmpc2p;imastmpc2=*imastmpc2p;
  auc2=*auc2p;adc2=*adc2p;irowc2=*irowc2p;icolc2=*icolc2p;jqc2=*jqc2p;
  aubd=*aubdp;irowbd=*irowbdp;jqbd=*jqbdp;
  aubdtil=*aubdtilp;irowbdtil=*irowbdtilp;jqbdtil=*jqbdtilp;
  aubdtil2=*aubdtil2p;irowbdtil2=*irowbdtil2p;jqbdtil2=*jqbdtil2p;
  audd=*auddp;irowdd=*irowddp;jqdd=*jqddp;
  auddtil=*auddtilp;irowddtil=*irowddtilp;jqddtil=*jqddtilp;
  auddtil2=*auddtil2p;irowddtil2=*irowddtil2p;jqddtil2=*jqddtil2p;
  auddinv=*auddinvp;irowddinv=*irowddinvp;jqddinv=*jqddinvp;
  irowtemp=*irowtempp;icoltemp=*icoltempp;jqtemp=*jqtempp;
  f_cs=*f_csp;f_cm=*f_cmp;
  ndirforc2=*ndirforc2p;nodeforc2=*nodeforc2p;xforc2=*xforc2p;
  
  fflush(stdout);
  
  NNEW(f_cs,double,neq[1]);
  NNEW(f_cm,double,neq[1]);
  
  // check for coupled thermo-mechanical calculation
  
  if(ithermal[0]>1){
    printf("\tprecontactmortar: coupled thermo-mechanical calculations NOT");
    printf("supported yet!\n \tPlease use surface-to-surface penalty contact");
    printf("instead.\n\n STOP!\n");
    fflush(stdout);
    FORTRAN(stop,());
  }
  
  // fix for linear calculation in first iteration of first increment
  
  if((*ilin==1)&&(*iit==1)&&(*iinc==1)){  
    *ielas=1;  
    iperturb[0]=-1;  
    iperturb[1]=0;	  
  }  
  
  /* small sliding is autometically set active due to combined fix-point
     Newton approach 
     do NOT change this unless the additional derivates neglected here 
     have been implemented */
  
  *iflagact=0;
  *ismallsliding=1;
  *nzsc2=nzs[1];
  NNEW(auc2,double,*nzsc2);
  NNEW(adc2,double,neq[1]);
  NNEW(irowc2,ITG,*nzsc2);
  NNEW(icolc2,ITG,neq[1]);
  NNEW(jqc2,ITG,neq[1]+1); 
  if(*iit==1 || *ismallsliding==0){
    NNEW(aubd,double,6*nslavnode[*ntie]);
    NNEW(irowbd,ITG,6*nslavnode[*ntie]);
    NNEW(jqbd,ITG,neq[1]+1);
    NNEW(aubdtil,double,6*nslavnode[*ntie]);
    NNEW(irowbdtil,ITG,6*nslavnode[*ntie]);
    NNEW(jqbdtil,ITG,neq[1]+1);
    NNEW(aubdtil2,double,6*nslavnode[*ntie]);
    NNEW(irowbdtil2,ITG,6*nslavnode[*ntie]);
    NNEW(jqbdtil2,ITG,neq[1]+1);
    NNEW(audd,double,3*nslavnode[*ntie]);
    NNEW(irowdd,ITG,3*nslavnode[*ntie]);
    NNEW(jqdd,ITG,neq[1]+1);
    NNEW(auddtil,double,3*nslavnode[*ntie]);
    NNEW(irowddtil,ITG,3*nslavnode[*ntie]);
    NNEW(jqddtil,ITG,neq[1]+1);
    NNEW(auddtil2,double,3*nslavnode[*ntie]);
    NNEW(irowddtil2,ITG,3*nslavnode[*ntie]);
    NNEW(jqddtil2,ITG,neq[1]+1);
    NNEW(auddinv,double,3*nslavnode[*ntie]);
    NNEW(irowddinv,ITG,3*nslavnode[*ntie]);
    NNEW(jqddinv,ITG,neq[1]+1);	    
  }
  
  /* storing the original stiffness matrix */
  
  NNEW(jqtemp,ITG,neq[1]+1);
  NNEW(irowtemp,ITG,nzs[1]);
  NNEW(icoltemp,ITG,neq[1]);
  for(i=0;i<3;i++){
    nzstemp[i]=nzs[i];}
  for (i=0;i<neq[1];i++){jqtemp[i]=jq[i];icoltemp[i]=icol[i];}
  jqtemp[neq[1]]=jq[neq[1]];
  for (i=0;i<nzs[1];i++){irowtemp[i]=irow[i];}
  
  if(*iit==1 || *ismallsliding==0){
    DMEMSET(slavnor,0,3*nslavnode[*ntie],0.);
    DMEMSET(slavtan,0,6*nslavnode[*ntie],0.);
  }
  
  /* setting iflagact=1 before calling contactmortar invokes
     combined fix-point Newton approach */
  
  if(*iit>1 && *ismallsliding==1){*iflagact=1;}
  
  /* transform SPCs/MPCs in case of quadratic finite elements
     Caution: there is still a problem with MPCs on slave mid nodes,
     avoid this!!! */
  
  //wird immer aufgerufen,da sich xbounact geaendert haben kann
  
  transformspcsmpcs_quad(nboun,ndirboun,nodeboun,xbounact,
			 nmpc,ipompc,nodempc,coefmpc,labmpc,
			 ikboun,ilboun,ikmpc,ilmpc,
			 nboun2,&ndirboun2,&nodeboun2,&xboun2,
			 nmpc2,&ipompc2,&nodempc2,&coefmpc2,&labmpc2,
			 &ikboun2,&ilboun2,&ikmpc2,&ilmpc2,
			 irowtlocinv,jqtlocinv,autlocinv,
			 nk,nk2,iflagdualquad,
			 ntie,tieset,itiefac,islavsurf,
			 lakon,ipkon,kon,&mt,memmpc_,
			 nodeforc,ndirforc,xforcact,nforc,
			 &nodeforc2,&ndirforc2,&xforc2,nforc2);

  RENEW(islavspc,ITG,2**nboun);
  RENEW(islavmpc,ITG,2**nmpc);
  RENEW(islavspc2,ITG,2**nboun2);
  RENEW(islavmpc2,ITG,2**nmpc2);
  RENEW(imastspc,ITG,2**nboun);
  RENEW(imastmpc,ITG,2**nmpc);
  RENEW(imastmpc2,ITG,2**nmpc);
  
  /* cataloque SPCs/MPCs */
  
  FORTRAN(catsmpcslavno,(ntie,islavnode,imastnode,nslavnode,
			 nmastnode,nboun,ndirboun,nodeboun,nmpc,
			 ipompc,nodempc,ikboun,ilboun,ikmpc,ilmpc,
			 nboun2,nmpc2,ipompc2,nodempc2,ikboun2,
			 ilboun2,ikmpc2,ilmpc2,nslavspc,islavspc,
			 nsspc,nslavmpc,islavmpc,nsmpc,nslavspc2,
			 islavspc2,nsspc2,nslavmpc2,islavmpc2,
			 nsmpc2,nmastspc,imastspc,nmspc,nmastmpc,
			 imastmpc,nmmpc,nmastmpc2,imastmpc2,
			 nmmpc2,jobnamef));
  
  RENEW(islavspc,ITG,2**nsspc+1);
  RENEW(islavmpc,ITG,2**nsmpc+1);
  RENEW(islavspc2,ITG,2**nsspc2+1);
  RENEW(islavmpc2,ITG,2**nsmpc2+1);
  RENEW(imastspc,ITG,2**nmspc+1);
  RENEW(imastmpc,ITG,2**nmmpc+1);
  RENEW(imastmpc2,ITG,2**nmmpc2+1);
  
  /* Check for additional MPCs on slave mid nodes */
  
  if(*iit==1 && *iinc==1){
      
    FORTRAN(nortanslav,(tieset,ntie,ipkon,kon,lakon,set,co,vold,nset,
			islavsurf,itiefac,islavnode,nslavnode,slavnor,slavtan,
			mi));
    
    // call checkspsmpc
    
    FORTRAN(checkspcmpc,(ntie,tieset,islavnode,imastnode,nslavnode,nmastnode,
			 slavnor,islavact,nboun,ndirboun,xboun,
			 nodempc,coefmpc,ikboun,ilboun,nmpc2,ipompc2,nodempc2,
			 nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
			 nmspc,nmastmpc,imastmpc,nmmpc));
        
    for (i=0;i<*ntie;i++){  
      if(tieset[i*(81*3)+80]=='C'){      
	for(j=nslavnode[i];j<nslavnode[i+1];j++){     
	  node=islavnode[j];
	  int checkformidnode=0;
	  for(jj=nslavmpc[2*(j)];jj<nslavmpc[2*(j)+1];jj++){
	    if(islavmpc[2*jj+1]==-2){
		  
	      // MPC cannot be identified
		  
	      checkformidnode=1;
	      if(debug==1)printf(" check for mid node %"ITGFORMAT"\n",node);
	    }	      
	  }
	    
	  if(checkformidnode==1){
	
	    //check for mid node
	
	    do{
	    	      
	      for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){
		ifaces = islavsurf[2*(l-1)+0];
		nelems = (ITG)(ifaces/10);
		jfaces = ifaces - nelems*10;
		FORTRAN(getnumberofnodes,(&nelems,&jfaces,lakon,&nope,&nopes,
					  &idummy)); 
		for(jj=0;jj<nope;jj++){
		  konl[jj]=kon[ipkon[nelems-1]+jj];
		}
		for(jj=0;jj<nopes;jj++){
		  jj2=jj+1;
		  ifac=FORTRAN(getlocno,(&jj2,&jfaces,&nope));
		  nodes[jj]=konl[ifac-1]; 
		}
		ii=-1;
		for(jj=0;jj<nopes;jj++){
		  if(nodes[jj]==node){ii=jj;}
		}
		if(ii>-1){
		  break;
		}
	      }

	      break;
	    }while(1);
      
	    if((ii>2 && nopes==6)||(ii>3 && nopes==8)){
	  
	      // mid node found with extra not supported mpc ->error
	  
	      printf(" precontactmortar: Problem with slave mid node  \n\n");
	      printf(" *ERROR: Slave mid node %"ITGFORMAT" has",node); 
	      printf(" additional MPC which is not a directional blocking ");
	      printf("MPC or a 1-to-1 cyclic symmetry MPC.  \n");
	      printf("\t\t This is not supported yet!!!!!!!!!!!!!\n");
	      fflush(stdout);
	      FORTRAN(stop,());
	    }
	  }        
	}
      }
    }
      
    if(*iit==1 || *ismallsliding==0){
      DMEMSET(slavnor,0,3*nslavnode[*ntie],0.);
      DMEMSET(slavtan,0,6*nslavnode[*ntie],0.);
    } 
    if(debug==1)printf(" precontactmortar: MPC-check OK\n\n");
  }
    
  /* fix for quadratic FE */
    
  NNEW(v,double,mt**nk);
  NNEW(stx,double,6*mi[0]**ne);
  NNEW(fn,double,mt*(*nk+*nk2));
  NNEW(fmpc2,double,*nmpc2);
  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
  *iout=-1;
  NNEW(ftil,double,neq[1]);
  NNEW(fexttil,double,neq[1]);
  NNEW(inum,ITG,*nk);
  
  /* calculating the internal forces and tangent stiffness using
     modified shape functions for quadratic elements */
  
  if(debug==1)printf(" precontactmortar: call results\n");
  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
        elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
        ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
        prestr,iprestr,filab,eme,emn,een,iperturb,
        ftil,fn,nactdof,iout,qa,vold,b,nodeboun,
        ndirboun,xbounact,nboun,ipompc,
        nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
        bet,gam,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
        xstateini,xstiff,xstate,npmat_,epn,matname,mi,ielas,icmd,
        ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
        xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
        ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc2,
        nelemload,nload,ikmpc,ilmpc,istep,iinc,springarea,
        reltime,ne0,thicke,shcon,nshcon,
        sideload,xloadact,xloadold,icfd,inomat,pslavsurf,pmastsurf,
        mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
        islavsurf,ielprop,prop,energyini,energy,kscale,iponoel,
        inoel,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	itiefac,tieset,smscale,mscalmethod,nbody,t0g,t1g,
        islavelinv,autloc,irowtloc,jqtloc,nboun2,
        ndirboun2,nodeboun2,xboun2,nmpc2,ipompc2,nodempc2,coefmpc2,
        labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	intscheme,physcon);
  
  if(debug==1)printf(" precontactmortar: results_dstil finished\n");
  fflush(stdout);
  
  SFREE(v);SFREE(stx);SFREE(fn);SFREE(inum);SFREE(fmpc2);
  *iout=0;	    
  
  *rhsi=1;
  NNEW(adtil,double,neq[1]);
  NNEW(autil,double,nzs[1]);
  NNEW(irowtil,ITG,nzs[1]);
  NNEW(jqtil,ITG,neq[1]+1);
  NNEW(icoltil,ITG,neq[1]);
  for(i=0;i<neq[1]+1;i++){
    jqtil[i]=jq[i];}
  for(i=0;i<neq[1];i++){
    icoltil[i]=icol[i];}
  for(i=0;i<nzs[1];i++){
    irowtil[i]=irow[i];}
  
  if(debug==1)printf(" precontactmortar: call mafillsmmain\n");
  fflush(stdout);
  
  /* calculating the external forces fext and stiffness matrix au/ad using
     modified shape functions for quadratic elements */
  
  mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun2,ndirboun2,xboun2,nboun2,
	       ipompc2,nodempc2,coefmpc2,nmpc2,nodeforc2,ndirforc2,xforc2,
	       nforc2,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
	       nbody,cgr,adtil,autil,fexttil,nactdof,icol,jqtil,irowtil,
	       neq,nzl,
	       nmethod,ikmpc2,ilmpc2,ikboun2,ilboun2,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	       ielmat,ielorien,norien,orab,ntmat_,
	       t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
	       nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	       xstiff,npmat_,dtime,matname,mi,
	       ncmat_,mass,stiffness,buckling,rhsi,intscheme,
	       physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
	       coriolis,ibody,xloadold,reltime,veold,springarea,nstate_,
	       xstateini,xstate,thicke,integerglob,doubleglob,
	       tieset,istartset,iendset,ialset,ntie,nasym,pslavsurf,
	       pmastsurf,mortar,clearini,ielprop,prop,ne0,fnext,kscale,
	       iponoel,inoel,network,ntrans,inotr,trab,smscale,
	       mscalmethod,set,nset,islavelinv,autloc,irowtloc,jqtloc,
	       &mortartrafoflag);
  
  if(debug==1)printf(" precontactmortar: mafillsmmain_dstil finished\n");
  fflush(stdout);
  
  for(i=0;i<neq[1];i++){
    ad[i]=adtil[i];
  }
  for(i=0;i<nzs[1];i++){
    au[i]=autil[i];
  }
  
  /* transform vold,vini,accold for dynamic calculations */
  
  NNEW(voldtil,double,mt**nk);
  NNEW(veoldtil,double,mt**nk);
  NNEW(vinitil,double,mt**nk);
  NNEW(accoldtil,double,mt**nk);
  NNEW(btil,double,neq[1]);
  
  if(debug==1)printf(" precontactmortar: call calcresidual\n");

  /* calculating the residual b */
  
  calcresidual(nmethod,neq,btil,fexttil,ftil,iexpl,nactdof,auxtil2,voldtil,
	       vinitil,dtime,accoldtil,nk,adbtil,aubtil,jqtil,irowtil,nzl,
	       &alpha,fextinitil,finitil,islavnode,nslavnode,mortar,ntie,f_cm,
	       f_cs,mi,nzs,nasym,idamping,veoldtil,adctil,auctil,cvinitil,cvtil,
	       alpham,num_cpus);
  
  if(debug==1)printf(" precontactmortar: calcresidual finished\n");
  fflush(stdout);  
	
  for(k=0;k<neq[1];k++){
    b[k]=btil[k];}
  
  SFREE(btil);SFREE(ftil);SFREE(fexttil);

  SFREE(voldtil);SFREE(vinitil);SFREE(accoldtil);SFREE(veoldtil);

  SFREE(adtil);SFREE(autil);SFREE(irowtil);SFREE(jqtil);SFREE(icoltil);
  
  /* update vold due to spcs to get gap right for rigid body movements */
  
  if(*iinc==1 && *iit==1 &&*nmethod!=4){       
    NNEW(v,double,mt**nk);	       
    NNEW(volddummy,double,mt**nk);
    NNEW(veolddummy,double,mt**nk);
    NNEW(accolddummy,double,mt**nk);	       
    for(k=0;k<mt**nk;k++){
      volddummy[k]=0.0;
      veolddummy[k]=0.0;
      accolddummy[k]=0.0;}           
    memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);	       
    NNEW(vectornull,double,neq[1]);	       
    *iout=-1;
    
    FORTRAN(resultsini_mortar,(nk,v,ithermal,iperturb,
			       nactdof,iout,volddummy,vectornull,nodeboun,ndirboun,
			       xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,
			       bet,gam,dtime,mi));
    
    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);	     	
    SFREE(v);SFREE(vectornull);SFREE(volddummy);SFREE(veolddummy);
    SFREE(accolddummy);	    
  }
  
  *ielas=0;
  *iout=0;
  iperturb[0]=iperturb_sav[0];
  iperturb[1]=iperturb_sav[1];
  mortartrafoflag=0;
  
  *ndirboun2p=ndirboun2;*nodeboun2p=nodeboun2;*xboun2p=xboun2;
  *ipompc2p=ipompc2;*nodempc2p=nodempc2;*coefmpc2p=coefmpc2;
  *labmpc2p=labmpc2;*ikboun2p=ikboun2;*ilboun2p=ilboun2;*ikmpc2p=ikmpc2;
  *ilmpc2p=ilmpc2;
  *nslavspcp=nslavspc;*islavspcp=islavspc;*nslavmpcp=nslavmpc;
  *islavmpcp=islavmpc;
  *nslavspc2p=nslavspc2;*islavspc2p=islavspc2;*nslavmpc2p=nslavmpc2;
  *islavmpc2p=islavmpc2;
  *nmastspcp=nmastspc;*imastspcp=imastspc;*nmastmpcp=nmastmpc;
  *nmastmpc2p=nmastmpc2;*imastmpcp=imastmpc;*imastmpc2p=imastmpc2;   
  *auc2p=auc2;*adc2p=adc2;*irowc2p=irowc2;*icolc2p=icolc2;*jqc2p=jqc2;
  *aubdp=aubd;*irowbdp=irowbd;*jqbdp=jqbd;
  *aubdtilp=aubdtil;*irowbdtilp=irowbdtil;*jqbdtilp=jqbdtil;
  *aubdtil2p=aubdtil2;*irowbdtil2p=irowbdtil2;*jqbdtil2p=jqbdtil2;
  *auddp=audd;*irowddp=irowdd;*jqddp=jqdd;
  *auddtilp=auddtil;*irowddtilp=irowddtil;*jqddtilp=jqddtil;
  *auddtil2p=auddtil2;*irowddtil2p=irowddtil2;*jqddtil2p=jqddtil2;
  *auddinvp=auddinv;*irowddinvp=irowddinv;*jqddinvp=jqddinv;
  *irowtempp=irowtemp;*icoltempp=icoltemp;*jqtempp=jqtemp;
  *f_csp=f_cs;*f_cmp=f_cm;
  *nodeforc2p=nodeforc2;*ndirforc2p=ndirforc2;
  *xforc2p=xforc2;
  
  return;
}
