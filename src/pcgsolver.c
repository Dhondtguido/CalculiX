/*

  pcgsolver.c	--	The preconditioned conjugate gradient solver(s)
  
  				Implemented solver/preconditioner:
  				    preconditioned cg-solver with
  				       * diagonal scaling only
  				       * incomplete Cholesky on full matrix
  				
  				Most functions are based on:
  				H.R. Schwarz: FORTRAN-Programme zur Methode 
                                             der finiten Elemente,
  				              Teubner, 1981 
 
                               The present version is based on the c-code by
                               Martin Ruecker and Ernst Rank of the Technical
                               University of Munich 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

#define GOOD 0
#define BAD 1
#define FALSE 0
#define TRUE 1

/* Prototyping	*/

ITG cgsolver (double *A, double *x, double *b, ITG neq, ITG len, ITG *ia, 
	       ITG *iz,double *eps, ITG *niter, ITG precFlg);
void PCG (double *A, double *x, double *b, ITG neq, ITG len, ITG *ia, 
	  ITG *iz,double *eps, ITG *niter, ITG precFlg,
	  double *rho, double *r, double *g, double *C, double *z);
void CG (double *A, double *x, double *b, ITG neq, ITG len, ITG *ia, 
	 ITG *iz,double *eps, ITG *niter,double *r, double *p, double *z);
void Scaling (double *A, double *b, ITG neq, ITG *ia, ITG *iz, double *d);
void MatVecProduct (double *A, double *p, ITG neq, ITG *ia, ITG *iz, 
		    double *z);
void PreConditioning (double *A, ITG neq, ITG len, ITG *ia, ITG *iz, 
		      double alpha, ITG precFlg,double *C, ITG *ier);
void Mrhor (double *C, ITG neq, ITG *ia, ITG *iz, double *r, double *rho);
void InnerProduct (double *a, double *b, ITG n, double *Sum);
	
/* **********************************************************************

The (preconditioned) conjugate gradient solver

  parameter:                                                                
           A       compact row oriented storage of lower left of matrix A  
           neq     order of A, number of equations                         
           len     number of non zero entries in Matrix A                  
           ia      column indices of corresponding elements in Matrix A    
           iz      row indices (diagonal elements indices)               
           x       solution vector                                       
           b       right hand side                                       
           eps     required accuracy -> residual                         
           niter   maximum number of iterations -> number of iterations  
           precFlg preconditioning flag                                  

The compact row oriented storage of sparse quadratic matrices is decsribed in
H.R. Schwarz: FORTRAN-Programme zur Methode der finiten Elemente, pp.66-67, 
Teubner, 1981
 
********************************************************************** 
*/

ITG cgsolver (double *A, double *x, double *b, ITG neq, ITG len, 
	       ITG *ia, ITG *iz, 
				double *eps, ITG *niter, ITG precFlg)
{
  ITG i=0;
  double *Factor=NULL,*r=NULL,*p=NULL,*z=NULL,*C=NULL,*g=NULL,*rho=NULL;

  /*  reduce row and column indices by 1 (FORTRAN->C)   */

  for (i=0; i<neq; i++)	--iz[i];
  for (i=0; i<len; i++)	--ia[i];

  /*  Scaling the equation system A x + b = 0  */

  NNEW(Factor,double,neq);
  Scaling(A,b,neq,ia,iz,Factor);

  /*  SOLVER/PRECONDITIONING TYPE  */

  /*  Conjugate gradient solver without preconditioning  */

  if (!precFlg)
    {
      NNEW(r,double,neq);
      NNEW(p,double,neq);
      NNEW(z,double,neq);
      CG(A,x,b,neq,len,ia,iz,eps,niter,r,p,z);
      SFREE(r);SFREE(p);SFREE(z);
    }
  
  /* Conjugate gradient solver with incomplete Cholesky preconditioning on
     full matrix */
  
  else if (precFlg==3)
    {
      NNEW(rho,double,neq);
      NNEW(r,double,neq);
      NNEW(g,double,neq);
      NNEW(C,double,len);
      NNEW(z,double,neq);
      PCG(A,x,b,neq,len,ia,iz,eps,niter,precFlg,rho,r,g,C,z);
      SFREE(rho);SFREE(r);SFREE(g);SFREE(C);SFREE(z);
    }
  
  /*  Backscaling of the solution vector  */
  
  for (i=0; i<neq; i++)	x[i] *= Factor[i];
  
  /*  That's it  */
  
  SFREE(Factor);
  return GOOD;
}

/* **********************************************************************

  Conjugate gradient solver with partial Cholesky preconditioning
  parameter:

           A       compact row oriented storage of lower left of the scaled
                   matrix A  
           x       solution vector 
           b       right hand side
           neq     order of A, number of equations 
           len     number of non zero entries in Matrix A 
           ia      column indices 
           iz      row indices (diagonal elements indices)
           eps     required accuracy -> residual
           niter   maximum number of iterations -> number of iterations
           precFlg preconditioning flag 

  The function corresponds to function PACHCG() in H.R. Schwarz: FORTRAN-Pro-
  gramme zur Methode der finiten Elemente, p.115, Teubner, 1981

********************************************************************** 
*/

void PCG (double *A, double *x, double *b, ITG neq, ITG len, ITG *ia, 
	  ITG *iz,double *eps, ITG *niter, ITG precFlg,
	  double *rho, double *r, double *g, double *C, double *z)
{
  ITG			i=0, k=0, ncg=0,iam,ier=0;
  double		alpha=0.0, ekm1=0.0, rrho=0.0;
  double		rrho1=0.0, gz=0.0, qk=0.0;
  double          c1=0.005,qam,err,ram=0;
  
  
  /*  initialize result and residual vectors  */
  
  qam=0.;iam=0;
  for (i=0; i<neq; i++)
    {
      x[i] = 0.0;	  
      r[i] = b[i];
      err=fabs(r[i]);
      if(err>1.e-20){qam+=err;iam++;}
    }
  if(iam>0) qam=qam/iam;
  else {*niter=0;return;}

  /*  preconditioning  */
  
  printf("Cholesky preconditioning\n\n");
  
  printf("alpha=%f\n",alpha);
  PreConditioning(A,neq,len,ia,iz,alpha,precFlg,C,&ier);
  while (ier==0)
    {
      if (alpha<=0.0)	alpha=0.005;
      alpha += alpha;
      printf("alpha=%f\n",alpha);
      PreConditioning(A,neq,len,ia,iz,alpha,precFlg,C,&ier);
    }
  
  /* solving the system of equations using the iterative solver */
  
  printf("Solving the system of equations using the iterative solver\n\n");
  
  /*  main iteration loop  */
  
  for (k=1; k<=*niter; k++, ncg++)
    {
      
      /*  solve M rho = r, M=C CT  */
      
      Mrhor(C,neq,ia,iz,r,rho);
      
      /*  inner product (r,rho)  */
      
      InnerProduct(r,rho,neq,&rrho);
      
      /*  If working with Penalty-terms for boundary conditions you can get 
	  numerical troubles because RRHO becomes a very large value. 
	  With at least two iterations the result may be better !!! */
      
      /*  convergency check */
      
      printf("iteration= %" ITGFORMAT ", error= %e, limit=%e\n",ncg,ram,c1*qam);
      if (k!=1 && (ram<=c1*qam))	break;
      if (k!=1)
	{
	  ekm1 = rrho/rrho1;
	  for (i=0; i<neq; i++)	g[i] = ekm1*g[i]-rho[i];
	}
      else
	{
	  
	  /*  g1 = -rho0  */
	  
	  for (i=0; i<neq; i++)	g[i] = -rho[i];
	}
      
      /*  solve equation system z = A g_k  */
      
      MatVecProduct(A,g,neq,ia,iz,z);
      
      /*  inner product (g,z)  */
      
      InnerProduct(g,z,neq,&gz);
      qk = rrho/gz;
      ram=0.;
      for (i=0; i<neq; i++)
	{
	  x[i] += qk*g[i];
	  r[i] += qk*z[i];
	  err=fabs(r[i]);
	  if(err>ram) ram=err;
	}
      rrho1 = rrho;
    }
  if(k==*niter){
    printf(" *ERROR in PCG: no convergence;");
    FORTRAN(stop,());
  } 
  *eps = rrho;
  *niter = ncg;
  
  return;
}


/* **********************************************************************

  Scaling the equation system A x + b = 0

  The equation system is scaled in consideration of keeping the symmetry in
  such a way that the diagonal elements of matrix A are 1. This is performed
  by the diagonal matrix Ds with the diagonal elements d_i = 1/sqrt(a_ii).
  The given equation system Ax+b=0 is transformed into
               -1                     - -   -
     Ds A Ds  Ds x + Ds b = 0   or    A x + b = 0
                         _                                             _ 
  with the scaled Matrix A= Ds A Ds  and  the scaled right hand side b= Ds b.
  The scaling factor Ds is needed later for backscaling of the solution 
  vector
     _    -1                       _ 
     x = Ds x      resp.    x = Ds x

  parameter:
           A       compact row oriented storage of lower left of matrix A
           b       right hand side
           neq     order of A, number of equations
           ia      column indices
           iz      row indices (diagonal elements indices)

  The function corresponds to function SCALKO() in H.R. Schwarz: FORTRAN-Pro-
  gramme zur Methode der finiten Elemente, p.105, Teubner, 1981

********************************************************************** 
*/

void Scaling (double *A, double *b, ITG neq, ITG *ia, ITG *iz, double *d)
{
  ITG			i=0, j=0, jlo=0, jup=0;
  
  /*  extract diagonal vector from matrix A  */
  
  for (i=0; i<neq; i++)		d[i] =  1.0/sqrt(A[iz[i]]);
  
  /*  scale right hand side (Ax=b -> Ax+b=0: negative sign)  */
  
  for (i=0; i<neq; i++)	b[i] *= -d[i];
  
  /*  scale matrix A  */
  
  A[iz[0]] *= d[0]*d[0];
  for (i=1; i<neq; i++)
    {
      jlo = iz[i-1]+1;
      jup = iz[i];
      for (j=jlo; j<=jup; j++)	A[j] *= d[i]*d[ia[j]];
    }
  
  return;
}

/* **********************************************************************

Computation of matrix vector product z = A p

  parameter:
             A       compact row oriented storage of lower left of matrix A
             p       vector
             neq     order of A, number of equations
             ia      column indices
             iz      row indices (diagonal elements indices)

  The function corresponds to function APZ() in H.R. Schwarz: FORTRAN-
  Programme zur Methode der finiten Elemente, p.111, Teubner, 1981

********************************************************************** 
*/

void MatVecProduct (double *A, double *p, ITG neq, ITG *ia, ITG *iz, 
		    double *z)
{
  ITG				i=0, j=0, jlo=0, jup=0, k=0;
  
  /*  matrix vector product  */
  
  z[0] = A[iz[0]]*p[0];
  for (i=1; i<neq; i++)
    {
      z[i] = A[iz[i]]*p[i];
      
      /*  first non-zero element in current row  */
      
      jlo = iz[i-1]+1;	
      
      /*  last non-zero off-diagonal element  */
      
      jup = iz[i]-1;				
      for (j=jlo; j<=jup; j++)
	{
	  k = ia[j];
	  z[i] += A[j]*p[k];
	  z[k] += A[j]*p[i];
	}
    }
  return;
}

/*--Preconditioning of the equation system Ax+b=0----------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--Partial Cholesky decomposition of scaled matrix A. The off-diagonal matrix   --	*/
/*--elements are divided by 1/(1+alpha) until a decomposition of the positive    --	*/
/*--definit matrix A exists. C is obtained by ignoring the fill-in during the    --	*/
/*--decomposition. In case of successfull decomposition the preconditioning      --	*/
/*--matrix C is returned. Otherwise the function is called again with new alpha. --	*/
/*--alpha has to be chosen as small as possible, because preconditioning effect  --	*/
/*--decreases with increasing alpha.-----------------------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--parameter:                                                                   --	*/
/*--           A       compact row oriented storage of lower left of matrix A    --	*/
/*--           neq     order of A, number of equations                           --	*/
/*--           ia      column indices                                            --	*/
/*--           iz      row indices (diagonal elements indices)                   --	*/
/*--           alpha   see above                                                 --	*/
/*--           precFlg ==3: partial Cholesky partition on full matrix A          --	*/
/*---------------------------------------------------------------------------------	*/
/*--The function corresponds to function PARTCH() in H.R. Schwarz: FORTRAN-Pro-  --	*/
/*--gramme zur Methode der finiten Elemente, p.117, Teubner, 1981                --	*/
/*---------------------------------------------------------------------------------	*/
void PreConditioning (double *A, ITG neq, ITG len, ITG *ia, ITG *iz, 
							double alpha, ITG precFlg,
		      double *C, ITG *ier)
{
	ITG				i=0, j=0, jlo=0, jup=0, k=0, klo=0, kup=0, l=0, llo=0, lup=0;
	ITG				id=0, nILU=0, m=0;
	double			factor;
		factor = 1.0/(1.0+alpha);
/*..division of the off-diagonal elements of A by factor (1.0+alpha)...............	*/
	C[0] = A[0];
	for (i=1; i<neq; i++)
	{
		jlo = iz[i-1]+1;			/*..first non-zero element in current row......	*/
		jup = iz[i];				/*..diagonal element...........................	*/
		C[jup] = A[jup];			/*..copy of diagonal element...................	*/
		for (j=jlo; j<jup; j++)		/*..all non-zero off-diagonal elements.........	*/
			C[j] = A[j]*factor;
	}
	nILU = neq;		/*..ILU on full matrix.................	*/
/*..partial Cholesky decomposition of scaled matrix A..............................	*/
	for (i=1; i<nILU; i++)
	{
		jlo = iz[i-1]+1;			/*..first non-zero element in current row......	*/
		jup = iz[i];				/*..diagonal element...........................	*/
		for (j=jlo; j<jup; j++)		/*..all non-zero off-diagonal elements.........	*/
		{
			C[j] /= C[iz[ia[j]]];
			klo = j+1;				/*..next non-zero element in current row.......	*/
			kup = jup;				/*..diagonal element in current row............	*/
			for (k=klo; k<=kup; k++)
			{
				m = ia[k];
				llo = iz[m-1]+1;
				lup = iz[m];
				for (l=llo; l<=lup; l++)
				{
					if (ia[l]>ia[j])	break;
					if (ia[l]<ia[j])	continue;
					C[k] -= C[j]*C[l];
					break;
				}
			}
		}
		id = iz[i];
		if (C[id]<1.0e-6){
		  return;
		}
		C[id] = sqrt(C[id]);
	}
	*ier=1;
	return;
}
/*---------------------------------------------------------------------------------	*/



/*--Solution of the equation system  M rho = r-------------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--The solution of the equation system M rho = r, with M = C*CT, where C is the --	*/
/*--partial Cholesky decomposition of matrix A, represent a preconditioning step.--	*/
/*--The equation  system is solved by forward/backward substitution.---------------	*/
/*---------------------------------------------------------------------------------	*/
/*--parameter:                                                                   --	*/
/*--           C       compact row oriented storage of preconditioned matrix A   --	*/
/*--           neq     order of A, number of equations                           --	*/
/*--           ia      column indices                                            --	*/
/*--           iz      row indices (diagonal elements indices)                   --	*/
/*--           r       residual vector                                           --	*/
/*---------------------------------------------------------------------------------	*/
/*--The function corresponds to function CRHOER() in H.R. Schwarz: FORTRAN-Pro-  --	*/
/*--gramme zur Methode der finiten Elemente, p.118, Teubner, 1981                --	*/
/*---------------------------------------------------------------------------------	*/
void Mrhor (double *C, ITG neq, ITG *ia, ITG *iz, double *r, double *rho)
{
	ITG				i=0, j=0, jlo=0, jup=0;
	double			s=0.0;
/*..solve equation sytem by forward/backward substitution..........................	*/
	rho[0] = r[0];
	for (i=1; i<neq; i++)
	{
		s = 0.0;
		jlo = iz[i-1]+1;			/*..first non-zero element in current row......	*/
		jup = iz[i];				/*..diagonal element in current row............	*/
		for (j=jlo; j<jup; j++)		/*..all non-zero off-diagonal element..........	*/
			s += C[j]*rho[ia[j]];
		rho[i] = (r[i]-s)/C[jup];
	}
	for (i=neq-1; i>0; i--)
	{
		rho[i] /= C[iz[i]];
		jlo = iz[i-1]+1;			/*..first non-zero element in current row......	*/
		jup = iz[i]-1;				/*..diagonal element in current row............	*/
		for (j=jlo; j<=jup; j++)	/*..all non-zero off-diagonal element..........	*/
			rho[ia[j]] -= C[j]*rho[i];
	}
	return;
}
/*---------------------------------------------------------------------------------	*/




/*--Calculation of the inner product of two (distributed) vectors------------------	*/
/*---------------------------------------------------------------------------------	*/
void InnerProduct (double *a, double *b, ITG n, double *Sum)
{
	ITG 		i=0;
/*..local vectors..................................................................	*/
	*Sum=0.;
		for (i=0; i<n; i++)	       *Sum += a[i]*b[i];
	return;
}
/*---------------------------------------------------------------------------------	*/




/*--(parallel) conjugate gradient solver-------------------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--parameter:                                                                   --	*/
/*--           A       compact row oriented storage of lower left of the scaled  --	*/
/*--                   matrix A                                                  --	*/
/*--           x       solution vector                                           --	*/
/*--           b       right hand side                                           --	*/
/*--           neq     order of A, number of equations                           --	*/
/*--           len     number of non zero entries in Matrix A                    --	*/
/*--           ia      column indices                                            --	*/
/*--           iz      row indices (diagonal elements indices)                   --	*/
/*--           eps     required accuracy -> residual                             --	*/
/*--           niter   maximum number of iterations -> number of iterations      --	*/
/*---------------------------------------------------------------------------------	*/
void CG (double *A, double *x, double *b, ITG neq, ITG len, ITG *ia, ITG *iz,
		double *eps, ITG *niter, double *r, double *p, double *z)
{
	ITG			i=0, k=0, ncg=0,iam;
	double		ekm1=0.0,c1=0.005,qam,ram=0.,err;
	double		rr=0.0, pz=0.0, qk=0.0, rro=0.0;


  /* solving the system of equations using the iterative solver */

  printf("Solving the system of equations using the iterative solver\n\n");

/*..initialize result, search and residual vectors.................................	*/
	qam=0.;iam=0;
	for (i=0; i<neq; i++)
	{
		x[i] =  0.0;					/*..start vector x0 = 0....................	*/
		r[i] =  b[i];					/*..residual vector r0 = Ax+b = b..........	*/
		p[i] = -r[i];					/*..initial search vector..................	*/
		err=fabs(r[i]);
		if(err>1.e-20){qam+=err;iam++;}
	}
	if(iam>0) qam=qam/neq;
	else {*niter=0;return;}
	/*else qam=0.01;*/
/*..main iteration loop............................................................	*/
	for (k=1; k<=(*niter); k++, ncg++)
	{
/*......inner product rT r......................................................... */
		InnerProduct(r,r,neq,&rr);
		printf("iteration= %" ITGFORMAT ", error= %e, limit=%e\n",ncg,ram,c1*qam);
/*......If working with Penalty-terms for boundary conditions you can get nume-....	*/
/*......rical troubles because RR becomes a very large value. With at least two....	*/
/*......iterations the result may be better !!!....................................	*/
/*......convergency check..........................................................	*/
		if (k!=1 && (ram<=c1*qam)) break;
/*......new search vector..........................................................	*/
		if (k!=1)
		{
			ekm1 = rr/rro;
			for (i=0; i<neq; i++)	p[i] = ekm1*p[i]-r[i];
		}
/*......matrix vector product A p = z..............................................	*/
		MatVecProduct(A,p,neq,ia,iz,z);
/*......inner product pT z.........................................................	*/
		InnerProduct(p,z,neq,&pz);
/*......step size..................................................................	*/
		qk = rr/pz;
/*......approximated solution vector, residual vector..............................	*/
		ram=0.;
		for (i=0; i<neq; i++)
		{
			x[i] = x[i] + qk*p[i];
			r[i] = r[i] + qk*z[i];
			err=fabs(r[i]);
			if(err>ram) ram=err;
		}
/*......store previous residual....................................................	*/
		rro = rr;

	}
	if(k==*niter){
	  printf(" *ERROR in PCG: no convergence;");
	  FORTRAN(stop,());
	} 
	*eps = rr;						/*..return residual............................	*/
	*niter = ncg;					/*..return number of iterations................	*/
/*..That's it......................................................................	*/
	return;
}
/*---------------------------------------------------------------------------------	*/
