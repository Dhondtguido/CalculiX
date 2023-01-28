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

void ini_cal(char *jobnamec,char *output,char *fneig,char *kind1,char *kind2,
	     ITG *itempuser,ITG *irobustdesign,ITG *nprint,
	     ITG *neq,ITG *mpcfree,ITG *nbounold,ITG *nforcold,ITG *nloadold,
	     ITG *nbody_,ITG *nbodyold,ITG *network,ITG *nheading_,ITG *nmpc_,
	     ITG *nload_,ITG *nforc_,ITG *nboun_,ITG *nintpoint,ITG *iperturb,
	     ITG *ntmat_,ITG *ithermal,ITG *isolver,ITG *nslavs,ITG *nkon_,
	     ITG *mortar,ITG *jout,ITG *nkon,ITG *nevtot,ITG *ifacecount,
	     ITG *iplas,ITG *npmat_,ITG *mi,ITG *mpcend,ITG *namtot_,
	     ITG *icascade,ITG *ne1d,ITG *ne2d,ITG *infree,ITG *nflow,
	     ITG *irstrt,ITG *nener,ITG *jrstrt,ITG *ntie_,ITG *mcs,ITG *nprop_,
	     ITG *nprop,ITG *itpamp,ITG *nevdamp_,ITG *npt_,ITG *iaxial,
	     ITG *inext,ITG *icontact,ITG *nobject,ITG *nobject_,ITG *iit,
	     ITG *mpcfreeref,ITG *isens,ITG *namtot,ITG *nstam,ITG *ndamp,
	     ITG *nef,ITG *nk_,ITG *ne_,ITG *nalset_,
	     ITG *nmat_,ITG *norien_,ITG *nam_,ITG *ntrans_,ITG *ncs_,
	     ITG *nstate_,ITG *ncmat_,ITG *memmpc_,ITG *nprint_,double *energy,
	     double *ctrl,double *alpha,double *qaold,double *physcon,
	     ITG *istep,ITG *istat,ITG *iprestr,ITG *kode,ITG *nload,
	     ITG *nbody,ITG *nforc,ITG *nboun,ITG *nk,ITG *nmpc,ITG *nam,
	     ITG *nzs_,ITG *nlabel,double *ttime,ITG *iheading,ITG *nfc,
	     ITG *nfc_,ITG *ndc,ITG *ndc_){

  /* used for initialization and re-initialization */
  
  ITG i;

  /* file names */
  
  for(i=132;i<792;i++) strcpy1(&jobnamec[i]," ",1);
  for(i=0;i<132;i++){strcpy1(&fneig[i]," ",1);}

  *iheading=0;
  *nheading_=0;

  /* parameters involved in reading the input deck */
  
  *istat=0;

  /* distritubed loading */
  
  *nloadold=0;
  *nload_=0;
  *nload=0;

  /* body loading */
  
  *nbody_=0;
  *nbodyold=0;
  *nbody=0;

  /* single point constraints */
  
  *nbounold=0;
  *nboun_=0;
  *nboun=0;

  /* concentrated loads */
  
  *nforcold=0;
  *nforc_=0;
  *nforc=0;

  /* temperature loading */
  
  for(i=0;i<2;i++) itempuser[i]=0;
  itempuser[2]=-2;

  /* node numbers */
  
  *nk_=0;
  *nk=0;

  /* information on multiple point constraints */
  
  *nmpc_=0;
  *nmpc=0;
  *memmpc_=0;
  *mpcfree=1;
  *mpcend=-1;
  *mpcfreeref=-1;
  *icascade=0;

  /* information on *coupling combined with *distributing */

  *nfc_=0;
  *nfc=0;
  *ndc_=0;
  *ndc=0;

  /* amplitude information */
  
  *namtot_=0;
  *namtot=0;
  *nam_=0;
  *nam=0;
  *itpamp=0;
  *inext=0;
  *nstam=0;

  /* ties / cyclic symmetry */
  
  *ntie_=0;
  strcpy1(kind1,"T",1);
  strcpy1(kind2,"T",1);
  *ncs_=0;
  *mcs=0;

  /* material information */
  
  *ntmat_=0;
  *npmat_=0;
  *nmat_=0;
  *norien_=0;
  *ntrans_=0;
  *nstate_=0;
  *ncmat_=0;
  *iplas=0;
  *ndamp=0;
  for(i=0;i<14;i++){physcon[i]=0.;}

  /* beam and fluid element properties */
  
  *nprop_=0;
  *nprop=0;

  /* sensitivity analysis and robust design */
  
  *isens=0;
  *nobject_=0;
  *nobject=0;
  for(i=0;i<3;i++) irobustdesign[i]=0;

  /* (modal) dynamic and frequency analysis */
  
  *nevdamp_=0;
  *nevtot=0;
  
  /* element topology */
  
  *ne_=0;
  *nkon_=0;
  *nkon=0;
  *nintpoint=0;
  *nef=0;

  /* 1d and 2d element information */
  
  *ne1d=0;
  *ne2d=0;
  *iaxial=1;
  for(i=0;i<4;i++) infree[i]=0;

  /* info on network elements */
  
  *nflow=0;
  *network=0;

  /* set information */
  
  *nalset_=0;

  /* output information */

  /* caveat: if changing next line:
     - change noelfiles appropriately
     - change nlabel in geomview.f, expand.c, storeresidual.f
     and createmddof.f
     - change the dimension of label in geomview.f
     - change the documentation (tex-file)  */

  *nlabel=55;
  *nprint_=0;
  *nprint=0;
  *kode=0;
  for(i=0;i<2;i++) jout[i]=1;
  strcpy1(output,"asc ",4);

  /* info on the kind of nonlinearity of the step */
  
  for(i=0;i<2;i++) iperturb[i]=0;

  /* ithermal[1] comes from readinput */
  
  ithermal[0]=0;

  /* restart */
  
  *jrstrt=0;
  for(i=0;i<2;i++) irstrt[i]=0;

  /* basic dimensional information */
  
  mi[0]=0;
  mi[1]=3;
  mi[2]=1;

  /* equation solver information */
  
  *isolver=0;
  *nzs_=20000000;
  for(i=0;i<3;i++) neq[i]=0;
  for(i=0;i<2;i++){qaold[i]=0.;}

  /* contact information */
  
  *ifacecount=0;
  *nslavs=0;
  *mortar=-2;
  *icontact=0;

  /* energy requests */
  
  *nener=0;
  for(i=0;i<5;i++){energy[i]=0.;}

  /* control parameters (cf. controlss.f) */

  /* time control parameters */
  
  ctrl[0]=4.5;    // i0
  ctrl[1]=8.5;    // ir
  ctrl[2]=9.5;    // ip
  ctrl[3]=16.5;   // ic
  ctrl[4]=10.5;   // il
  ctrl[5]=4.5;    // ig
  ctrl[6]=0.;     // is
  ctrl[7]=5.5;    // ia
  ctrl[8]=0.;     // ij
  ctrl[9]=0.;     // it
  ctrl[10]=0.25;  // df
  ctrl[11]=0.5;   // dc
  ctrl[12]=0.75;  // db
  ctrl[13]=0.85;  // da
  ctrl[14]=0.;    // ds
  ctrl[15]=0.;    // dh
  ctrl[16]=1.5;   // dd
  ctrl[17]=0.;    // wg

  /* field control parameters */
  
  ctrl[18]=0.005; // ran
  ctrl[19]=0.01;  // can
  ctrl[20]=0.;    // qa0
  ctrl[21]=0.;    // qau
  ctrl[22]=0.02;  // rap
  ctrl[23]=1.e-5; // ea
  ctrl[24]=1.e-3; // cae
  ctrl[25]=1.e-8; // ral

  /* heat transfer and modal dynamics parameter */
  
  ctrl[26]=1.e30; // deltmx

  /* line search parameters */
  
  ctrl[27]=1.5;   // nls
  ctrl[28]=0.25;  // smaxls
  ctrl[29]=1.01;  // sminls
  ctrl[30]=1.;    // fls
  ctrl[31]=1.;    // etls

  /* network parameters */
  
  ctrl[32]=5.e-7; // c1t
  ctrl[33]=5.e-7; // c1f
  ctrl[34]=1.e-4; // c1p
  ctrl[35]=5.e-7; // c2t
  ctrl[36]=5.e-7; // c2f
  ctrl[37]=5.e-7; // c2p
  ctrl[38]=5.e-7; // c2a

  /* viscous parameter */
  
  ctrl[39]=-1.;   // cetol

  /* network parameters */
  
  ctrl[40]=1.e20; // a1t
  ctrl[41]=1.e20; // a1f
  ctrl[42]=1.e20; // a1p
  ctrl[43]=1.e20; // a2t
  ctrl[44]=1.e20; // a2f
  ctrl[45]=1.e20; // a2p
  ctrl[46]=1.e20; // a2a

  /* CFD parameters */
  
  ctrl[47]=1.5;   // scheme
  ctrl[48]=0.5;   // simple or simplec
  ctrl[49]=20.5;  // iitf
  ctrl[50]=0.5;   // iitg
  ctrl[51]=1.5;   // iitp
  ctrl[52]=1.5;   // iitpt

  /* contact parameters */
  
  ctrl[53]=0.001; // delcon
  ctrl[54]=0.1;   // alea
  ctrl[55]=100.5; // kscalemax
  ctrl[56]=60.5;  // itf2f

  /* time integration dynamics */
  
  alpha[0]=0.;
  alpha[1]=0.5;

  /* information on the step, iteration ... */
  
  *istep=0;
  *iit=-1;
  *ttime=0.;

  /* residual stress, pre-stress */
  
  *iprestr=0;
  *npt_=0;
  
  /* default solver */

#if defined(SGI)
  *isolver=4;
#elif defined(PASTIX)
  *isolver=8;
#elif defined(PARDISO)
  *isolver=7;
#elif defined(SPOOLES)
  *isolver=0;
#elif defined(TAUCS)
  *isolver=5;
#else
  *isolver=3;
#endif
  
  return;
}
