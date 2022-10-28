/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2021 Guido Dhondt                     */

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
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif

void massless(ITG *kslav,ITG *lslav,ITG *ktot,ITG *ltot,double *au,double *ad,
	      double *auc,double *adc,
	      ITG *jq,ITG *irow,ITG *neq,ITG *nzs,double *auw,
	      ITG *jqw,ITG *iroww,ITG *nzsw,ITG *islavnode,ITG *nslavnode,
	      ITG *nslavs,ITG *imastnode,ITG *nmastnode,ITG *ntie,ITG *nactdof,
	      ITG *mi,double *vold,double *volddof, double *veold,ITG *nk,
	      double *fext,ITG *isolver,ITG *masslesslinear,double *co,
	      double *springarea,ITG *neqtot,double *qb,
	      double *b,double *dtime,double *aloc,double *fric,ITG *iexpl,
	      ITG *nener,double *ener,ITG *ne,ITG **jqbip,double **aubip,
	      ITG **irowbip,ITG **jqibp,double **auibp,ITG **irowibp,
	      ITG *iclean,ITG *iinc){

  /* determining the RHS of the global system for massless contact */

  ITG *jqwnew=NULL,*irowwnew=NULL,symmetryflag=0,mt=mi[1]+1,
    inputformat=0,*iacti=NULL,nacti=0,itranspose,index,i,j,k,kitermax,
    *jqbb=NULL,*irowbb=NULL,*icolbb=NULL,nzsbb,*jqbi=NULL,*irowbi=NULL,
    nzsbi,*jqib=NULL,*irowib=NULL,nzsib,nrhs=1;

  double *auwnew=NULL,sigma=0.0,*gapdisp=NULL,*gapnorm=NULL,*cvec=NULL,sum,
    *adbbb=NULL,*aubbb=NULL,*gvec=NULL,*gmatrix=NULL,*qi_kbi=NULL,
    *veolddof=NULL,*alglob=NULL,atol,rtol,*aubb=NULL,*adbb=NULL,
    *al=NULL,*alnew=NULL,*eps_al=NULL,*rhs=NULL,*aubi=NULL,
    *auib=NULL,omega,*alocold=NULL,*ddisp=NULL;

  jqbi=*jqbip;aubi=*aubip;irowbi=*irowbip;jqib=*jqibp;auib=*auibp;
  irowib=*irowibp;

  /* clearing memory for the equation solver at the end of the step*/

  if(*iclean==1){
    if(*isolver==0){
#ifdef SPOOLES
      spooles_cleanup_rad();
#endif
    }else if(*isolver==7){
#ifdef PARDISO
      pardiso_cleanup_cp(neqtot,&symmetryflag,&inputformat);
#endif
    }else if(*isolver==8){
#ifdef PASTIX
      pastix_solve_cp(qb,neqtot,&symmetryflag,&nrhs);
#endif
    }
    return;
  }
  
  /* expanding the matrix Wb according to the number of degrees
     of freedom */

  NNEW(jqwnew,ITG,3**nslavs+1);
  NNEW(auwnew,double,*nzsw);
  NNEW(irowwnew,ITG,*nzsw);

  /* Rearrange the row entries in the Wb matrix, column by column
     from the order in islavnode and imastnode to the order as
     dictated by nactdof */
  
  FORTRAN(expand_auw,(auw,jqw,iroww,nslavs,auwnew,jqwnew,irowwnew,
		      nactdof,mi,ktot,neqtot,islavnode,imastnode));

  memcpy(jqw,jqwnew,sizeof(ITG)*(3**nslavs+1));
  memcpy(auw,auwnew,sizeof(double)**nzsw);
  memcpy(iroww,irowwnew,sizeof(ITG)**nzsw);

  SFREE(jqwnew);SFREE(auwnew);SFREE(irowwnew);

  /* extracting Kbb,Kbi,Kib,Kii from the stiffness matrix;
     only for a nonlinear calculation or the first increment
     of a linear calculation */

  if((*masslesslinear==0)||(*iinc==1)){
    NNEW(jqbb,ITG,*neqtot+1);
    NNEW(aubb,double,nzs[0]);
    NNEW(adbb,double,*neqtot);
    NNEW(irowbb,ITG,nzs[0]);
    NNEW(icolbb,ITG,*neqtot);

    NNEW(jqbi,ITG,neq[0]+1);
    NNEW(aubi,double,nzs[0]);
    NNEW(irowbi,ITG,nzs[0]);

    NNEW(jqib,ITG,*neqtot+1);
    NNEW(auib,double,nzs[0]);
    NNEW(irowib,ITG,nzs[0]);

    FORTRAN(extract_matrices,(au,ad,jq,irow,neq,aubb,adbb,jqbb,irowbb,neqtot,
			      &nzsbb,aubi,jqbi,irowbi,&nzsbi,auib,jqib,irowib,
			      &nzsib,ktot,icolbb));

    RENEW(aubb,double,nzsbb);
    RENEW(irowbb,ITG,nzsbb);
    RENEW(aubi,double,nzsbi);
    RENEW(irowbi,ITG,nzsbi);
    RENEW(auib,double,nzsib);
    RENEW(irowib,ITG,nzsib);
  }

  /* calculate the residual force in the contact area and store in gapdisp */
  
  NNEW(gapdisp,double,*neqtot);
  NNEW(qi_kbi,double,*neqtot);
  
  FORTRAN(resforccont,(vold,nk,mi,aubi,irowbi,jqbi,neqtot,ktot,fext,gapdisp,
		       auib,irowib,jqib,nactdof,volddof,neq,qi_kbi));

  /* factorize Kbb and delete Kbb afterwards;
     only for a nonlinear calculation or the first increment
     of a linear calculation */

  if((*masslesslinear==0)||(*iinc==1)){
    if(*isolver==0){
#ifdef SPOOLES
      spooles_factor_rad(adbb,aubb,adbbb,aubbb,&sigma,icolbb,irowbb,
			 neqtot,&nzsbb,&symmetryflag,&inputformat,iexpl);
#else
      printf(" *ERROR in massless: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }else if(*isolver==7){
#ifdef PARDISO
      pardiso_factor_cp(adbb,aubb,adbbb,aubbb,&sigma,icolbb,
			irowbb,neqtot,&nzsbb,&symmetryflag,&inputformat,jqbb,
			&nzsbb,iexpl);
#else
      printf(" *ERROR in massless: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }else if(*isolver==8){
#ifdef PASTIX
      ITG inputformat=1;
      pastix_factor_main_cp(adbb,aubb,adbbb,aubbb,&sigma,icolbb,
			    irowbb,neqtot,&nzsbb,&symmetryflag,&inputformat,
			    jqbb,&nzsbb);
#else
      printf(" *ERROR in massless: the PASTIX library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    SFREE(aubb);SFREE(adbb);SFREE(irowbb);SFREE(icolbb);SFREE(jqbb);
  }

  /* premultiply gapdisp with Kbb^{-1} */

  if(*isolver==0){
#ifdef SPOOLES
    spooles_solve_rad(gapdisp,neqtot);
#endif
  }else if(*isolver==7){
#ifdef PARDISO
    pardiso_solve_cp(gapdisp,neqtot,&symmetryflag,&inputformat,&nrhs);
#endif
  }else if(*isolver==8){
#ifdef PASTIX
    pastix_solve_cp(gapdisp,neqtot,&symmetryflag,&nrhs);
#endif
  }

  NNEW(gapnorm,double,*nslavs);
  NNEW(iacti,ITG,*neqtot);

  /* premultiply g by Wb^T and add g0 => determine active degrees => 
     reduce g to c */

  FORTRAN(detectactivecont,(gapnorm,gapdisp,auw,iroww,jqw,nslavs,springarea,
			    iacti,&nacti));

  /* reduced the gap dimension to the active dofs */
  
  if (nacti>0){

    /* constructing the c-vector of the inclusion equation */
    
    NNEW(cvec,double,nacti);

    for (i=0;i<*neqtot;i++){
      gapdisp[i]=(gapdisp[i]-volddof[ktot[i]-1])/(*dtime); // gapdisp - qb_km1
    }

    // cvec = Wb^T * cvec
    
    for (i=0;i<3**nslavs;++i){
      if (iacti[i]!=0){
	//        index=i;
	//        for (j=jqw[index]-1;j<jqw[index+1]-1;j++){
        for (j=jqw[i]-1;j<jqw[i+1]-1;j++){
          cvec[iacti[i]-1]+=auw[j]*gapdisp[iroww[j]-1];
        }
      }
    }

    /* constructing the g-matrix of the inclusion equation */
    NNEW(gmatrix,double,nacti*nacti);

    /* calculate G = Wb^T.Kbb^(-1).Wb 
       only for active slave degrees of freedom */

    /* loop over the columns of Wb */
    
    for(i=0;i<3**nslavs;i++){
      
      if(iacti[i]!=0) {
	//	index=i;//  contact index

	NNEW(gvec,double,*neqtot);

	/* Filling the vector of Wb column */
	
	for(j=jqw[i]-1;j<jqw[i+1]-1;j++){
	  gvec[iroww[j]-1]=auw[j];
	}
 
	/* Solving the linear system per column */

	if(*isolver==0){
#ifdef SPOOLES
	  spooles_solve_rad(gvec,neqtot);
#endif
	}else if(*isolver==7){
#ifdef PARDISO
	  pardiso_solve_cp(gvec,neqtot,&symmetryflag,&inputformat,&nrhs);
#endif
	}else if(*isolver==8){
#ifdef PASTIX
	  pastix_solve_cp(gvec,neqtot,&symmetryflag,&nrhs);
#endif
	}
	
	/* premultiplying per Wb^T */

	for(j=0;j<3**nslavs;++j){
	  if(iacti[j]!=0){
	    //	    index=j;
	    sum=0.0;
	    for(k=jqw[j]-1;k<jqw[j+1]-1;k++){
	      sum+=auw[k]*gvec[iroww[k]-1];
	    }
	    gmatrix[(iacti[i]-1)*nacti+(iacti[j]-1)]= sum / (*dtime) ;
	  }
	}
	SFREE(gvec);
      }
    }

    /* solve the inclusion problem (augmented Lagrange) */
    
    atol=1.0e-8;
    rtol=1.0e-6;
    kitermax=5000;
    
    //modifier for relaxation: According to observations in Monjaraz 2022,
    // ML contact  may be 2x more aggressive in convergence relaxation
    // than what is shown in Studer 2009.
    // For full FE it should play a big role in performance, unless the total of contact DOF is very high.
    // We control this parameter with OMEGA. Recommended is 100% (no modification) No more than 200% or less than 100%.
    // should not be necesary to decrease it < 100%
    
    omega=1.0;

    NNEW(al,double,nacti);
    NNEW(alnew,double,nacti);
    
    // taking values from aloc for initial guess
    
    for (i=0;i<3**nslavs;++i){
      if (iacti[i]!=0){
	al[iacti[i]-1]=aloc[i];
	alnew[iacti[i]-1]=aloc[i];
      }
    }

    NNEW(eps_al,double,nacti);
    NNEW(alglob,double,*neqtot);
    FORTRAN(auglag_inclusion,
	    (gmatrix,cvec,iacti,&nacti,fric,&atol,&rtol,
	     alglob,&kitermax,auw,jqw,iroww,nslavs,al,
	     alnew,eps_al,&omega));

    if(*nener==1){
      NNEW(alocold,double,3**nslavs);
      memcpy(&alocold[0],&aloc[0],sizeof(double)*3**nslavs);
    }
    
    // storing back values for next initial guess
    
    for(i=0;i<3**nslavs;++i){
      if(iacti[i]!=0){
	aloc[i]=alnew[iacti[i]-1];
      }
    }
    SFREE(al);
    SFREE(alnew);
    SFREE(eps_al);
    SFREE(gmatrix);
    SFREE(cvec);
    
  }else{
    NNEW(alglob,double,*neqtot);
  }

  SFREE(gapdisp); // TODO CMT move this freeing, we need it still
  SFREE(gapnorm);
  
  /* compute  qb = Kbb^{-1}*(Wb*al-qi_kbi+fexb) */

  for(i=0;i<*neqtot;i++){
    qb[i]=alglob[i]-qi_kbi[i]+fext[ktot[i]-1];
  }

  if(*isolver==0){
#ifdef SPOOLES
    spooles_solve_rad(qb,neqtot);
#endif
  }else if(*isolver==7){
#ifdef PARDISO
    pardiso_solve_cp(qb,neqtot,&symmetryflag,&inputformat,&nrhs);
#endif
  }else if(*isolver==8){
#ifdef PASTIX
    pastix_solve_cp(qb,neqtot,&symmetryflag,&nrhs);
#endif
  }

  SFREE(alglob);SFREE(qi_kbi);

  /* compute the energy dissipated by friction */

  if(*nener==1){

    /* relative tangential displacements in this increment
       in a local coordinate system: Wb^T*(qb^{n+1}-qb^n) */

    NNEW(ddisp,double,3**nslavs);
    for(i=0;i<3**nslavs;i++){

      /* treat only tangential direction */

      if(3*(i/3)!=i){

	/* treat only active slave nodes */
	
	if(iacti[i]!=0){
	  for(j=jqw[i]-1;j<jqw[i+1]-1;j++){
	    ddisp[i]+=auw[j]*(qb[iroww[j]-1]-volddof[ktot[iroww[j]-1]-1]);
	  }
	}
      }
    }

    /* friction energy increase = (al^{n+1}+al^n)/2*Wb^T*(qb^{n+1}-qb^n) */
    
    for(i=0;i<3**nslavs;i++){
      if(3*(i/3)!=i){
	if(iacti[i]!=0){
	  ener[2*mi[0]*(*ne+i/3)+1]+=-(alocold[i]+aloc[i])*ddisp[i]/2.;
	}
      }
    }
    SFREE(alocold);SFREE(ddisp);
  }
  
  SFREE(iacti);
  
  /* compute the right hand side of the global equation system */
  
  /* calculate Kii*qi */
  
  NNEW(rhs,double,neq[0]); 
  FORTRAN(op,(&neq[0],volddof,rhs,ad,au,jq,irow));

  /* calculate Kii*qi+Kib*qb */

  itranspose=0;
  FORTRAN(mulmatvec_asym,(auib,jqib,irowib,neqtot,qb,rhs,&itranspose));
  itranspose=1;
  FORTRAN(mulmatvec_asym,(aubi,jqbi,irowbi,&neq[0],qb,rhs,&itranspose));

  /* deleting the matrices Kbi and Kib only in the nonlinear case */
  
  if(*masslesslinear==0){
    SFREE(jqbi);SFREE(aubi);SFREE(irowbi);SFREE(jqib);SFREE(auib);SFREE(irowib);
  }
    
  /* calculate (Mii/Delta_t-Dii/2)*u_i^{k-1/2} */

  /* switch for the velocity from a nodal to a dof representation */
  
  NNEW(veolddof,double,neq[0]);
  for(k=0;k<*nk;++k){
    for(j=0;j<mt;j++){
      if(nactdof[mt*k+j]>0){
        veolddof[nactdof[mt*k+j]-1]=veold[mt*k+j];
      }
    }
  }

  FORTRAN(op,(&neq[0],veolddof,b,adc,auc,jq,irow));

  for(i=0;i<neq[0];++i){b[i]=fext[i]-rhs[i]+b[i];}
  SFREE(rhs);
  SFREE(veolddof);

  /* clearing memory for the equation solver;
     only in the nonlinear case */

  if(*masslesslinear==0){
    if(*isolver==0){
#ifdef SPOOLES
      spooles_cleanup_rad();
#endif
    }else if(*isolver==7){
#ifdef PARDISO
      pardiso_cleanup_cp(neqtot,&symmetryflag,&inputformat);
#endif
    }else if(*isolver==8){
#ifdef PASTIX
      pastix_solve_cp(qb,neqtot,&symmetryflag,&nrhs);
#endif
    }
  }

  *jqbip=jqbi;*aubip=aubi;*irowbip=irowbi;*jqibp=jqib;*auibp=auib;
  *irowibp=irowib;

  return;
}





