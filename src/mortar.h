/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License,or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,  */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not,write to the Free Software       */
/*     Foundation,Inc.,675 Mass Ave,Cambridge,MA 02139,USA.         */

void add_rect(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1,ITG m_1,
	      double *au_2,ITG * irow_2,ITG * jq_2,ITG n_2,ITG m_2,
	      double **au_rp,ITG **irow_rp,ITG * jq_r,ITG *nzs);
       
void bdfill(ITG **irowbdp,ITG *jqbd,double **aubdp,ITG *nzsbd,
	    ITG **irowbdtilp,ITG *jqbdtil,double **aubdtilp,ITG *nzsbdtil,
	    ITG **irowbdtil2p,ITG *jqbdtil2,double **aubdtil2p,ITG *nzsbdtil2,
	    ITG **irowddp,ITG *jqdd,double **auddp,
	    ITG **irowddtilp,ITG *jqddtil,double **auddtilp,
	    ITG **irowddtil2p,ITG *jqddtil2,double **auddtil2p,
	    ITG **irowddinvp,ITG *jqddinv,double **auddinvp,
	    ITG *irowt,ITG *jqt,double *aut,
	    ITG *irowtinv,ITG *jqtinv,double *autinv,
	    ITG *ntie,ITG *ipkon,ITG *kon,
	    char *lakon,ITG *nslavnode,ITG *nmastnode,ITG *imastnode,
	    ITG *islavnode,ITG *islavsurf,ITG *imastsurf,double *pmastsurf,
	    ITG *itiefac,char *tieset,ITG *neq,ITG *nactdof,double *co,
	    double *vold,
	    ITG *iponoels,ITG *inoels,ITG *mi,double *gapmints,double *gap,
	    double* pslavsurf,double* pslavdual,
	    ITG *nintpoint,double *slavnor,ITG *nk,
	    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
	    ITG *ikmpc,ITG *ilmpc,ITG *nslavmpc,ITG *islavmpc,
	    ITG *nmastmpc,ITG *imastmpc,
	    ITG *iit,ITG *iinc,ITG *islavactdof,ITG *islavact,ITG *islavnodeinv,
	    double **Bdp,ITG **irowbp,ITG *jqb,
	    double **Bdhelpp,ITG **irowbhelpp,ITG *jqbhelp,
	    double **Ddp,ITG **irowdp,ITG *jqd,
	    double **Ddtilp,ITG **irowdtilp,ITG *jqdtil,
	    double **Bdtilp,ITG **irowbtilp,ITG *jqbtil,
	    ITG *ithermal);  
    
void buildtquad(ITG *ntie,ITG *ipkon,ITG *kon,ITG *nk,
		char *lakon,ITG *nslavnode,
		ITG *itiefac,char *tieset,
		ITG *islavnode,ITG *islavsurf,
		ITG **irowtp,ITG *jqt,double **autp,
		ITG **irowtinvp,ITG *jqtinv,double **autinvp);
    
void FORTRAN(remlagrangemult,(ITG *ntie,char *tieset,ITG *islavnode,
			      ITG *imastnode,ITG *nslavnode,ITG *nmastnode,
			      ITG *islavact,ITG *nodempc,ITG *nmpc,
			      ITG *ipompc));
    
void contactmortar(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
		   ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
		   char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,ITG *ne,
		   double *cg,double *straight,double *co,
		   double *vold,ITG *ielmat, double *elcon,
		   ITG *istep,ITG *iinc,ITG *iit,ITG *ncmat_,ITG *ntmat_,
		   ITG *ne0,double *vini,
		   ITG *nmethod,ITG *neq,ITG *nzs,ITG *nactdof,ITG *itiefac,
		   ITG *islavsurf,ITG *islavnode,ITG *imastnode,
		   ITG *nslavnode,ITG *nmastnode,double *ad,
		   double **aup,double *b,ITG **irowp,ITG *icol,ITG *jq,
		   ITG *imastop,
		   ITG *iponoels,ITG *inoels,ITG *nzsc,double **aucp,
		   double *adc,ITG **irowcp,ITG *jqc,ITG *islavact,
		   double *gap,
		   double *slavnor,double *slavtan,
		   double *bhat,
		   ITG **irowbdp,ITG *jqbd,double **aubdp,
		   ITG **irowbdtilp,ITG *jqbdtil ,double **aubdtilp,
		   ITG **irowbdtil2p,ITG *jqbdtil2,double **aubdtil2p,
		   ITG **irowddp,ITG *jqdd,double **auddp,
		   ITG **irowddtilp,ITG *jqddtil,double **auddtilp,
		   ITG **irowddtil2p,ITG *jqddtil2,double **auddtil2p,
		   ITG **irowddinvp,ITG *jqddinv,double **auddinvp,
		   ITG *irowt,ITG *jqt,double *aut,  
		   ITG *irowtinv,ITG *jqtinv,double *autinv,   
		   ITG *mi,ITG *ipe,ITG *ime,double *tietol,
		   double *cstress,
		   double *cstressini,double *bp_old,ITG *nk,
		   ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		   ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		   ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		   ITG *nslavspc,ITG *islavspc,ITG *nslavmpc,
		   ITG *islavmpc,ITG *nmastmpc,ITG *imastmpc,
		   double *pslavdual,ITG *islavactdof,ITG *islavtie,
		   double *plicon,ITG *nplicon,ITG *npmat_,ITG *nelcon,
		   double *dtime,
		   ITG *islavnodeinv,
		   double **Bdp,ITG **irowbp,ITG *jqb,
		   double **Bdhelpp,ITG **irowbhelpp,ITG *jqbhelp,
		   double **Ddp,ITG **irowdp,ITG *jqd,
		   double **Ddtilp,ITG **irowdtilp,ITG *jqdtil,
		   double **Bdtilp,ITG **irowbtilp,ITG *jqbtil,
		   double *bet,
		   double *cfsinitil,double *reltime,
		   ITG *ithermal,double *plkcon,ITG *nplkcon);

void stressmortar(double *bhat,double *adc,double *auc,ITG *jqc,
		  ITG *irowc,ITG *neq,double *gap,double *b,ITG *islavact,
		  ITG *irowddinv,ITG *jqddinv,double *auddinv,
		  ITG *irowt,ITG *jqt,double *aut, 
		  ITG *irowtinv,ITG *jqtinv,double *autinv,
		  ITG *ntie,ITG *nslavnode,
		  ITG *islavnode,ITG *nmastnode,ITG *imastnode,double *slavnor,
		  double *slavtan,
		  ITG *nactdof,ITG *iflagact,double *cstress,
		  double *cstressini,ITG *mi,
		  double *cdisp,double *f_cs,double *f_cm,ITG *iit,ITG *iinc,
		  double *vold,double *vini,double* bp,ITG *nk,
		  ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		  ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		  ITG *nslavmpc,
		  ITG *islavmpc,
		  char *tieset,
		  double  *elcon,double *tietol,ITG *ncmat_,ITG *ntmat_,
		  double *plicon,ITG *nplicon,ITG *npmat_,ITG *nelcon,
		  double *dtime,double *cfs,double *cfm,ITG *islavnodeinv,
		  double *Bd,ITG *irowb,ITG *jqb,
		  double *Dd,ITG *irowd,ITG *jqd,
		  double *Ddtil,ITG *irowdtil,ITG *jqdtil,
		  double *Bdtil,ITG *irowbtil,ITG *jqbtil,
		  ITG *nmethod,
		  double *bet,ITG *ithermal,ITG *iperturb,
		  char *labmpc,double *cam,double *veold,
		  double *accold,double *gam,
		  double *cfsini,double *cfstil,double *plkcon,ITG *nplkcon,
		  char *filab,double *f,double *fn,
		  double *qa,ITG *nprint,char *prlab,double *xforc,
		  ITG *nforc,ITG *iponoel);
    
void FORTRAN(spcmpcmortar,(ITG *ntie,ITG *islavnode,ITG *imastnode,
			   ITG *nslavnode,ITG *nmastnode,ITG *nboun,
			   ITG *nmpc,ITG *ipompc,
			   ITG *nodempc,ITG *ikbou,ITG *ilbou,ITG *ikmpc,
			   ITG *ilmpc,ITG *nslavspc,ITG *islavspc,ITG *nsspc,
			   ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,
			   ITG *nmastspc,ITG *imastspc,ITG *nmspc,
			   ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
			   char *jobnamef));

void FORTRAN(create_t,(ITG *ipkon,ITG *kon,char *lakon,ITG *islavsurf,
			 double *dcontr,ITG *idcontr1,ITG *idcontr2,
			 ITG *icounter,ITG *l));

void FORTRAN(create_t_lin,(ITG *ipkon,ITG *kon,char *lakon,ITG *islavsurf,
			     double *dcontr,ITG *idcontr1,ITG *idcontr2,
			     ITG *icounter,ITG *l));

void FORTRAN(create_tinv,(ITG *ipkon,ITG *kon,char *lakon,ITG *islavsurf,
			    double *dcontr,ITG *idcontr1,ITG *idcontr2,
			    ITG *icounter,ITG *l));

void FORTRAN(create_tinv_lin,(ITG *ipkon,ITG *kon,char *lakon,ITG *islavsurf,
				double *dcontr,ITG *idcontr1,ITG *idcontr2,
				ITG *icounter,ITG *l));
       
void FORTRAN(createbd,(ITG *ict,ITG *l,ITG *ipkon,ITG *kon,char *lakon,
		       double *co,double *vold,double* gapmints,ITG *islavsurf,
		       ITG *imastsurf,double *pmastsurf,double *contr,
		       ITG *isconstr,ITG *imcontr,double *dcontr,ITG *idcontr1,
		       ITG *idcontr2,double *gcontr,ITG *igcontr,ITG *mi,
		       double* pslavsurf,double* pslavdual,ITG *nslavnode,
		       ITG *islavnode,ITG *nmastnode,ITG *imastnode,
		       ITG *icounter,ITG *icounter2,ITG *islavact));

void decascade_mortar(ITG *nmpc,ITG *ipompc,ITG **nodempcp,double **coefmpcp,
		      ITG *ikmpc,ITG *ilmpc,ITG *memmpc_,ITG *mpcfree);

void FORTRAN(nortanslav,(char *tieset,ITG *ntie,ITG *ipkon,ITG *kon,
			 char *lakon,char *set,double *co,double *vold,
			 ITG *nset,ITG *islavsurf,ITG *itiefac,ITG *islavnode,
			 ITG *nslavnode,double *slavnor,double *slavtan,
			 ITG *mi));

void FORTRAN(gendualcoeffs,(char *tieset,ITG *ntie,ITG *ipkon,ITG *kon,
			    char *lakon,double *co,double *vold,ITG *islavact,
			    ITG *islavsurf,ITG *itiefac,ITG *islavnode,
			    ITG *nslavnode,ITG *mi,double *pslavsurf,
			    double* pslavdual));

void FORTRAN(genfirstactif,(char *tieset,ITG *ntie,ITG *itietri,
			    double *cg,double *straight,
			    double *co,double *vold,double *xo,double *yo,
			    double *zo,double *x,double *y,double *z,ITG *nx,
			    ITG *ny,ITG *nz,ITG *istep,ITG *iinc,ITG *iit,
			    ITG *mi,ITG *imastop,ITG *nslavnode,ITG *islavnode,
			    char *set,ITG *nset,
			    ITG *istartset,ITG *iendset,ITG *ialset,
			    ITG *islavact,double *tietol));

void FORTRAN(genislavactdof,(ITG *ntie,char *tieset,ITG *nactdof,
			     ITG *nslavnode,ITG *nmastnode,ITG *imastnode,
			     ITG *islavactdof,ITG *islavnode,ITG *mi,
			     ITG *ithermal));
     
void FORTRAN(genislavquadel,(ITG *islavquadel,ITG *jqt,char *lakon,
			     ITG *ipkon,ITG *kon,ITG *ne,ITG *nasym,
			     ITG *nslavquadel));

void FORTRAN(getcontactparams,(double *mu,ITG *regmode,
			       double *fkninv,double *fktauinv,
			       double *p0,double *beta,double *tietol,
			       double *elcon,ITG *itie,
			       ITG *ncmat_,ITG *ntmat_));
     
ITG FORTRAN(getlocno,(ITG *m,ITG *jfaces,ITG *nope));

void FORTRAN(getnumberofnodes,(ITG *nelems,ITG *jfaces,char *lakon,ITG *nope,
			       ITG *nopes,ITG *idummy)); 

void inimortar(double **enerp,ITG *mi,ITG *ne ,ITG *nslavs,ITG *nk,ITG *nener,
	       ITG **ipkonp,char **lakonp,ITG **konp,ITG *nkon,
	       ITG *maxprevcontel,double **xstatep,ITG *nstate_,
	       ITG **islavtiep,double **bpp,ITG **islavactp,
	       double **gapp,
	       double **cdispp,double **cstressp,double **cfsp,
	       double **bpinip,ITG **islavactinip,double **cstressinip,
	       ITG *ntie,char *tieset,ITG *nslavnode,ITG *islavnode,
	       ITG **islavnodeinvp,ITG **islavelinvp,double **pslavdualp,
	       double **autp,ITG **irowtp,ITG **jqtp,
	       double **autinvp,ITG **irowtinvp,ITG **jqtinvp,
	       double **Bdp,ITG **irowbp,ITG **jqbp,
	       double **Bdhelpp,ITG **irowbhelpp,ITG **jqbhelpp,
	       double **Ddp,ITG **irowdp,ITG **jqdp,
	       double **Ddtilp,ITG **irowdtilp,ITG **jqdtilp,
	       double **Bdtilp,ITG **irowbtilp,ITG **jqbtilp,
	       ITG *itiefac,ITG *islavsurf,ITG *nboun,ITG *nmpc,
	       ITG **nslavspcp,ITG **islavspcp,ITG **nslavmpcp,ITG **islavmpcp,
	       ITG **nmastspcp,ITG **imastspcp,ITG **nmastmpcp,ITG **imastmpcp,
	       ITG *imastnode,ITG *nmastnode,ITG *nasym,ITG *mortar,
	       ITG **ielmatp,ITG **ielorienp,ITG *norien,ITG *ipompc,
	       ITG *nodempc,ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       char *jobnamef,char *set,double *co,double *vold,ITG *nset,
	       ITG *nslavquadel);

void insertas(ITG **irowp,ITG **mast1p,ITG *i1,
	      ITG *i2,ITG *ifree,ITG *nzs_,double *contribution,double **bdp);

void insertas_ws(ITG **irowp, ITG *i1,ITG *i2, ITG *ifree, ITG *nzs_,
		 double *contribution, double **bdp);

void matrixsort(double *au,ITG *mast1,ITG *irow,ITG *jq,
		ITG *nzs,ITG *dim);
      
void mortar_prefrd(ITG *ne,ITG *nslavs,ITG *mi,ITG *nk,ITG *nkon,
		   double **stxp,double *cdisp,
		   double *fn,double *cfs,double *cfm);
   
void mortar_postfrd(ITG *ne,ITG *nslavs,ITG *mi,ITG *nk,ITG *nkon,
		    double *fn,double *cfs,double *cfm);
       
void multimortar(double **aup,double *ad,ITG **irowp,ITG *jq,ITG *nzs,
		 double **aucp,double *adc,ITG **irowcp,ITG *jqc,ITG *nzsc,
		 double *aubd,ITG *irowbd,ITG *jqbd,
		 double *aubdtil,ITG *irowbdtil,ITG *jqbdtil,
		 double *aubdtil2,ITG *irowbdtil2,ITG *jqbdtil2,
		 ITG *irowdd,ITG *jqdd,double *audd,
		 ITG *irowddtil2,ITG *jqddtil2,double *auddtil2,
		 ITG *irowddinv,ITG *jqddinv,double *auddinv,
		 double *Bd,ITG *irowb,ITG *jqb,
		 double *Dd,ITG *irowd,ITG *jqd,
		 double *Ddtil,ITG *irowdtil,ITG *jqdtil,
		 ITG *neq,double *b,double *bhat,ITG *islavnode,ITG *imastnode,
		 ITG *nslavnode,ITG *nmastnode,
		 ITG *islavact,ITG *islavactdof,
		 double *gap,
		 double *slavnor,double *slavtan,
		 double *vold,double *vini,double *cstress,double *cstressini,
		 double *bp_old,ITG *nactdof,ITG *ntie,ITG *mi,ITG *nk,
		 ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		 ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		 ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		 ITG *nslavspc,ITG *islavspc,ITG *nslavmpc,
		 ITG *islavmpc,char *tieset,
		 ITG *islavtie,ITG *nelcon,double  *elcon,double *tietol,
		 ITG *ncmat_,ITG *ntmat_,
		 double *plicon,ITG *nplicon,ITG *npmat_,double *dtime,
		 ITG *irowt,ITG *jqt,double *aut,
		 ITG *irowtinv,ITG *jqtinv,double *autinv,
		 ITG *islavnodeinv,ITG *iit,ITG *nmethod,double *bet,
		 ITG *ithermal,
		 double *plkcon,ITG *nplkcon
		 );

void multi_rect(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1,ITG m_1,
		double *au_2,ITG * irow_2,ITG * jq_2,ITG n_2,ITG m_2,
		double **au_rp,ITG **irow_rp,ITG * jq_r,ITG *nzs);

void multi_rectv(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1,ITG m_1,
		 double * b,double ** v_rp);

void multi_scal(double *au_1,ITG * irow_1,ITG * jq_1,
		double *au_2,ITG * irow_2,ITG * jq_2,
		ITG m,ITG n,double*value,ITG *flag);
       
void FORTRAN(opnonsym,(ITG *neq,double *aux,double *b,double *bhat,
		       double *bdd,double*bdu,ITG *jqbd,ITG *irowbd));

void FORTRAN(opnonsymt,(ITG *neq,double *aux,double *b,double *bhat,
			double *bdd,double*bdu,ITG *jqbd,ITG *irowbd));

void premortar(ITG *nzs,ITG *nzsc2,
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
	       ITG *iit,
	       ITG *icol,ITG *irow,ITG *jq,
	       ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       ITG *imastnode,ITG *nmastnode,
	       double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,double *stn,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,double *prestr,
	       ITG *iprestr,char *filab,double *eme,double *emn,
	       double *een,ITG *iperturb,ITG *nactdof,
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
	       double *aut,ITG *irowt,ITG *jqt,
	       double *autinv,ITG *irowtinv,ITG *jqtinv,
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
	       double *alpham,double *betam,
	       double *pslavsurf,double *pmastsurf,
	       double *clearini,ITG *ielprop,double *prop,
	       ITG *islavact,double *cdn,ITG *memmpc_,
	       ITG *idamping,
	       ITG *ilin,ITG *iperturb_sav,
	       ITG *itietri,double *cg,double *straight,ITG *koncont,
	       double *energyini,
	       double *energy,ITG *kscale,ITG *iponoeln,ITG *inoeln,ITG *nener,
	       char *orname,ITG *network,
	       char *typeboun,ITG *num_cpus,double *t0g,double *t1g,
	       double *smscale,ITG *mscalmethod,ITG *nslavquadel,
	       ITG *iponoel);
       
void FORTRAN(regularization_gn_c,(double *lambdap,ITG *divmode,ITG *regmode,
				  double *gnc,double *aninvloc,double *p0,
				  double *beta,double *elcon,ITG *nelcon,
				  ITG *itie,ITG *ntmat_,double *plicon,
				  ITG *nplicon,ITG *npmat_,ITG *ncmat_,
				  double *tietol,double *scal));
     
void FORTRAN(regularization_gt_c,(double *lambdatt,ITG *divmode,
				  double *gtc,double *atauinvloc));
 
void FORTRAN(regularization_slip_lin,(double *utilt,double *bp,double *atauinv,
				      double *resreg,ITG *divmode,
				      ITG *islavact,double *lambdat,
				      double *lambdatilt,double *constantt,
				      ITG *inode,double *n2,
				      double *t,double *that,double *mu,
				      double *rslip,double *ltslip,
				      double *ltu));
     
void FORTRAN(resultsini_mortar,(int *nk,double *v,int *ithermal,int *iperturb,
				int *nactdof,int *iout,double *vold,double *b,
				int *nodeboun,int *ndirboun,double *xboun,
				int *nboun,int *ipompc,int *nodempc,
				double *coefmpc,char *labmpc,int *nmpc,
				int *nmethod,double *cam,double *bet,
				double *gam,double *dtime,int *mi));     

void FORTRAN(slavintmortar,(ITG *ntie,ITG *itietri,ITG *ipkon,ITG *kon,
			    char *lakon,double *straight,ITG *nintpoint,
			    ITG *koncont,double *co,double *vold,double *xo,
			    double *yo,double *zo,double *x,double *y,
			    double *z,ITG *nx,ITG *ny,ITG *nz,ITG *iinc,
			    ITG *islavsurf,ITG *imastsurf,double *pmastsurf,
			    ITG *islavnode,ITG *nslavnode,ITG *imastop,
			    double *gap,ITG *islavact,ITG *mi,ITG *ncont,
			    ITG *ipe,ITG *ime,double *pslavsurf,ITG *i,ITG *l,
			    ITG *ntri,double *tietol,double *reltime,
			    ITG *nmethod));
    
void transpose(double *au,ITG *jq,ITG *irow,ITG *dim,
	       double *au_t,ITG *jq_t,ITG *irow_t);

void trafontmortar2(ITG *neq,ITG *nzs,ITG *islavactdof,ITG *islavact,
		    ITG *nslavnode,
		    double *f_da,double *f_atil,
		    double *au_dan,ITG *irow_dan,ITG *jq_dan,
		    double *au_dam,ITG *irow_dam,ITG *jq_dam,
		    double *au_dai,ITG *irow_dai,ITG *jq_dai,
		    double *au_daa,ITG *irow_daa,ITG *jq_daa,
		    double **au_antilp,ITG **irow_antilp,ITG *jq_antil,
		    double **au_amtilp,ITG **irow_amtilp,ITG *jq_amtil,
		    double **au_aitilp,ITG **irow_aitilp,ITG *jq_aitil,
		    double **au_aatilp,ITG **irow_aatilp,ITG *jq_aatil,
		    double *gap,
		    double *Bd,ITG *irowb,ITG *jqb,
		    double *Dd,ITG *irowd,ITG *jqd,
		    double *Ddtil,ITG *irowdtil,ITG *jqdtil,
		    double *au_bdtil2,ITG *irow_bdtil2,ITG *jq_bdtil2,
		    double *au_ddtil2i,ITG *irow_ddtil2i,ITG *jq_ddtil2i,
		    double *au_ddtil2a,ITG *irow_ddtil2a,ITG *jq_ddtil2a,
		    ITG *m_flagr,ITG *i_flagr,ITG *a_flagr,ITG *a_flag,
		    ITG *i_flag,ITG *m_flag,
		    ITG *row_ln,ITG *row_lm,ITG *row_li,ITG *row_la,
		    double *slavnor,double *slavtan,
		    double *vold,double *vini,double *cstress,
		    double *cstressini,
		    double *bp_old,ITG *nactdof,ITG *islavnode,
		    ITG *ntie,ITG *mi,ITG *nk,
		    ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		    ITG *ipompc,ITG *nodempc,double *coefmpc,
		    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		    ITG *nslavspc,ITG *islavspc,ITG *nslavmpc,
		    ITG *islavmpc,char *tieset,
		    ITG *islavtie,ITG *nelcon,double  *elcon,
		    double *tietol,ITG *ncmat_,ITG *ntmat_,
		    double *plicon,ITG *nplicon,ITG *npmat_,double *dtime,
		    ITG *irowt,ITG *jqt,double *aut,
		    ITG *irowtinv,ITG *jqtinv,double *autinv,
		    ITG *islavnodeinv,
		    ITG *iit,double *beta,ITG *ithermal,
		    double *plkcon,ITG *nplkcon);
    
void trafontspcmpc( double *n,double *t,double *n2,double *that,
		    ITG *islavnodeentry,
		    ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		    ITG *ipompc,ITG *nodempc,double *coefmpc,
		    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		    ITG *nslavspc,ITG *islavspc,ITG *nslavmpc,
		    ITG *islavmpc,ITG *node);

void FORTRAN(updatecont,(ITG *koncont,ITG *ncont,double *co,double *vold,
			 double *cg,double *straight,ITG *mi));
 
void FORTRAN(writematrix,(double *au,double *ad,ITG *irow,ITG *jq,ITG *neq,
			  ITG *number));
 
void FORTRAN(writematrix2,(double *au,double *ad,ITG *irow,ITG *jq,ITG *neq,
			   ITG *number));

void FORTRAN(writevector,(double *ad,ITG *neq,ITG *number));
