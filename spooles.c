
/*     CalculiX - A 3-dimensional finite element program                   */
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
/*
 * The implementation is derived from the SPOOLES sample described in
 * AllInOne.ps
 *  created -- 98jun04, cca
 *
 * Converted to something that resembles C and
 * support for multithreaded solving added.
 * (C) 2003 Manfred Spraul
 */

/* spooles_factor and spooles_solve occur twice in this routine: once
   with their plane names and once with _rad appended to the name. This is
   necessary since the factorized stiffness matrices (plain names) and the
   factorized radiation matrices (_rad appended) are kept at the same time
   in the program */

#ifdef SPOOLES

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "CalculiX.h"
#include "spooles.h"

#if USE_MT
int num_cpus = -1;
#endif

#define TUNE_MAXZEROS		1000
#define TUNE_MAXDOMAINSIZE	800
#define TUNE_MAXSIZE		64

#define RNDSEED		7892713
#define MAGIC_DTOL		0.0
#define MAGIC_TAU		100.0

/*
 * Substeps for solving A X = B:
 *
 *  (1) form Graph object
 *  (2) order matrix and form front tree
 *  (3) get the permutation, permute the matrix and 
 *      front tree and get the symbolic factorization
 *  (4) compute the numeric factorization
 *  (5) read in right hand side entries
 *  (6) compute the solution
 *
 * The ssolve_main functions free the input matrices internally
 */

static void ssolve_creategraph(Graph ** graph, ETree ** frontETree,
			 InpMtx * mtxA, int size, FILE * msgFile)
{
	IVL *adjIVL;
	int nedges;

	*graph = Graph_new();
	adjIVL = InpMtx_fullAdjacency(mtxA);
	nedges = IVL_tsize(adjIVL);
	Graph_init2(*graph, 0, size, 0, nedges, size, nedges, adjIVL,
		    NULL, NULL);
	if (DEBUG_LVL > 1) {
		fprintf(msgFile, "\n\n graph of the input matrix");
		Graph_writeForHumanEye(*graph, msgFile);
		fflush(msgFile);
	}
	/* (2) order the graph using multiple minimum degree */

	/*maxdomainsize=neqns/100; */
	/*if (maxdomainsize==0) maxdomainsize=1; */
	/* *frontETree = orderViaMMD(*graph, RNDSEED, DEBUG_LVL, msgFile) ; */
	/* *frontETree = orderViaND(*graph,maxdomainsize,RNDSEED,DEBUG_LVL,msgFile); */
	/* *frontETree = orderViaMS(*graph,maxdomainsize,RNDSEED,DEBUG_LVL,msgFile); */

	*frontETree =
	    orderViaBestOfNDandMS(*graph, TUNE_MAXDOMAINSIZE,
				  TUNE_MAXZEROS, TUNE_MAXSIZE, RNDSEED,
				  DEBUG_LVL, msgFile);
	if (DEBUG_LVL > 1) {
		fprintf(msgFile, "\n\n front tree from ordering");
		ETree_writeForHumanEye(*frontETree, msgFile);
		fflush(msgFile);
	}
}

static void ssolve_permuteA(IV ** oldToNewIV, IV ** newToOldIV,
			 IVL ** symbfacIVL, ETree * frontETree,
			 InpMtx * mtxA, FILE * msgFile, int *symmetryflagi4)
{
	int *oldToNew;

	*oldToNewIV = ETree_oldToNewVtxPerm(frontETree);
	oldToNew = IV_entries(*oldToNewIV);
	*newToOldIV = ETree_newToOldVtxPerm(frontETree);
	ETree_permuteVertices(frontETree, *oldToNewIV);
	InpMtx_permute(mtxA, oldToNew, oldToNew);
	if(*symmetryflagi4!=2) InpMtx_mapToUpperTriangle(mtxA);
	InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS);
	InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);
	*symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA);
	if (DEBUG_LVL > 1) {
		fprintf(msgFile, "\n\n old-to-new permutation vector");
		IV_writeForHumanEye(*oldToNewIV, msgFile);
		fprintf(msgFile, "\n\n new-to-old permutation vector");
		IV_writeForHumanEye(*newToOldIV, msgFile);
		fprintf(msgFile, "\n\n front tree after permutation");
		ETree_writeForHumanEye(frontETree, msgFile);
		fprintf(msgFile, "\n\n input matrix after permutation");
		InpMtx_writeForHumanEye(mtxA, msgFile);
		fprintf(msgFile, "\n\n symbolic factorization");
		IVL_writeForHumanEye(*symbfacIVL, msgFile);
		fflush(msgFile);
	}
}

static void ssolve_postfactor(FrontMtx *frontmtx, FILE *msgFile)
{
	FrontMtx_postProcess(frontmtx, DEBUG_LVL, msgFile);
	if (DEBUG_LVL > 1) {
		fprintf(msgFile, "\n\n factor matrix after post-processing");
		FrontMtx_writeForHumanEye(frontmtx, msgFile);
		fflush(msgFile);
	}
}

static void ssolve_permuteB(DenseMtx *mtxB, IV *oldToNewIV, FILE* msgFile)
{
	DenseMtx_permuteRows(mtxB, oldToNewIV);
	if (DEBUG_LVL > 1) {
		fprintf(msgFile,
			"\n\n right hand side matrix in new ordering");
		DenseMtx_writeForHumanEye(mtxB, msgFile);
		fflush(msgFile);
	}
}

static void ssolve_permuteout(DenseMtx *mtxX, IV *newToOldIV, FILE *msgFile)
{
	DenseMtx_permuteRows(mtxX, newToOldIV);
	if (DEBUG_LVL > 1) {
		fprintf(msgFile, "\n\n solution matrix in original ordering");
		DenseMtx_writeForHumanEye(mtxX, msgFile);
		fflush(msgFile);
	}
}

 void factor(struct factorinfo *pfi, InpMtx *mtxA, int size, FILE *msgFile,
             int *symmetryflagi4)
{
	Graph *graph;
	IVL *symbfacIVL;
	Chv *rootchv;

	/* Initialize pfi: */
	pfi->size = size;
	pfi->msgFile = msgFile;
	pfi->solvemap = NULL;
	DVfill(10, pfi->cpus, 0.0);

	/*
	 * STEP 1 : find a low-fill ordering
	 * (1) create the Graph object
	 */
	ssolve_creategraph(&graph, &pfi->frontETree, mtxA, size, pfi->msgFile);

	/*
	 * STEP 2: get the permutation, permute the matrix and 
	 *      front tree and get the symbolic factorization
	 */
	ssolve_permuteA(&pfi->oldToNewIV, &pfi->newToOldIV, &symbfacIVL, pfi->frontETree,
		     mtxA, pfi->msgFile, symmetryflagi4);

	/*
	 * STEP 3: initialize the front matrix object
	 */
	{
		pfi->frontmtx = FrontMtx_new();
		pfi->mtxmanager = SubMtxManager_new();
		SubMtxManager_init(pfi->mtxmanager, NO_LOCK, 0);
		FrontMtx_init(pfi->frontmtx, pfi->frontETree, symbfacIVL, SPOOLES_REAL,
			      *symmetryflagi4, FRONTMTX_DENSE_FRONTS,
			      SPOOLES_PIVOTING, NO_LOCK, 0, NULL,
			      pfi->mtxmanager, DEBUG_LVL, pfi->msgFile);
	}

	/* 
	 * STEP 4: compute the numeric factorization
	 */
	{
		ChvManager *chvmanager;
		int stats[20];
		int error;

		chvmanager = ChvManager_new();
		ChvManager_init(chvmanager, NO_LOCK, 1);
		IVfill(20, stats, 0);
		rootchv = FrontMtx_factorInpMtx(pfi->frontmtx, mtxA, MAGIC_TAU, MAGIC_DTOL,
						chvmanager, &error, pfi->cpus,
						stats, DEBUG_LVL, pfi->msgFile);
		ChvManager_free(chvmanager);
		if (DEBUG_LVL > 1) {
			fprintf(msgFile, "\n\n factor matrix");
			FrontMtx_writeForHumanEye(pfi->frontmtx, pfi->msgFile);
			fflush(msgFile);
		}
		if (rootchv != NULL) {
			fprintf(pfi->msgFile, "\n\n matrix found to be singular\n");
			exit(-1);
		}
		if (error >= 0) {
			fprintf(pfi->msgFile, "\n\nerror encountered at front %d",
				error);
			exit(-1);
		}
	}
	/*
	 * STEP 5: post-process the factorization
	 */
	ssolve_postfactor(pfi->frontmtx, pfi->msgFile);

	/* cleanup: */
	IVL_free(symbfacIVL);
	InpMtx_free(mtxA);
	Graph_free(graph);
}

DenseMtx *fsolve(struct factorinfo *pfi, DenseMtx *mtxB)
{
	DenseMtx *mtxX;
	/*
	 * STEP 6: permute the right hand side into the new ordering
	 */
	{
		DenseMtx_permuteRows(mtxB, pfi->oldToNewIV);
		if (DEBUG_LVL > 1) {
			fprintf(pfi->msgFile,
				"\n\n right hand side matrix in new ordering");
			DenseMtx_writeForHumanEye(mtxB, pfi->msgFile);
			fflush(pfi->msgFile);
		}
	}
	/*
	 * STEP 7: solve the linear system
	 */
	{
		mtxX = DenseMtx_new();
		DenseMtx_init(mtxX, SPOOLES_REAL, 0, 0, pfi->size, 1, 1, pfi->size);
		DenseMtx_zero(mtxX);
		FrontMtx_solve(pfi->frontmtx, mtxX, mtxB, pfi->mtxmanager, pfi->cpus,
			       DEBUG_LVL, pfi->msgFile);
		if (DEBUG_LVL > 1) {
			fprintf(pfi->msgFile, "\n\n solution matrix in new ordering");
			DenseMtx_writeForHumanEye(mtxX, pfi->msgFile);
			fflush(pfi->msgFile);
		}
	}
	/*
	 * STEP 8:  permute the solution into the original ordering
	 */
	ssolve_permuteout(mtxX, pfi->newToOldIV, pfi->msgFile);

	/* cleanup: */
	DenseMtx_free(mtxB);

	return mtxX;
}

#ifdef USE_MT 
void factor_MT(struct factorinfo *pfi, InpMtx *mtxA, int size, FILE *msgFile, int *symmetryflagi4)
{
	Graph *graph;
	IV *ownersIV;
	IVL *symbfacIVL;
	Chv *rootchv;

	/* Initialize pfi: */
	pfi->size = size;
	pfi->msgFile = msgFile;
	DVfill(10, pfi->cpus, 0.0);

	/*
	 * STEP 1 : find a low-fill ordering
	 * (1) create the Graph object
	 */
	ssolve_creategraph(&graph, &pfi->frontETree, mtxA, size, msgFile);

	/*
	 * STEP 2: get the permutation, permute the matrix and 
	 *      front tree and get the symbolic factorization
	 */
	ssolve_permuteA(&pfi->oldToNewIV, &pfi->newToOldIV, &symbfacIVL, pfi->frontETree,
		     mtxA, msgFile, symmetryflagi4);

	/*
	 * STEP 3: Prepare distribution to multiple threads/cpus
	 */
	{
		DV *cumopsDV;
		int nfront;

		nfront = ETree_nfront(pfi->frontETree);

		pfi->nthread = num_cpus;
		if (pfi->nthread > nfront)
			pfi->nthread = nfront;

		cumopsDV = DV_new();
		DV_init(cumopsDV, pfi->nthread, NULL);
		ownersIV = ETree_ddMap(pfi->frontETree, SPOOLES_REAL, *symmetryflagi4,
				       cumopsDV, 1. / (2. * pfi->nthread));
		if (DEBUG_LVL > 1) {
			fprintf(msgFile,
				"\n\n map from fronts to threads");
			IV_writeForHumanEye(ownersIV, msgFile);
			fprintf(msgFile,
				"\n\n factor operations for each front");
			DV_writeForHumanEye(cumopsDV, msgFile);
			fflush(msgFile);
		} else {
			fprintf(msgFile, "\n\n Using %d threads\n",
				pfi->nthread);
		}
		DV_free(cumopsDV);
	}

	/*
	 * STEP 4: initialize the front matrix object
	 */
	{
		pfi->frontmtx = FrontMtx_new();
		pfi->mtxmanager = SubMtxManager_new();
		SubMtxManager_init(pfi->mtxmanager, LOCK_IN_PROCESS, 0);
		FrontMtx_init(pfi->frontmtx, pfi->frontETree, symbfacIVL, SPOOLES_REAL,
			      *symmetryflagi4, FRONTMTX_DENSE_FRONTS,
			      SPOOLES_PIVOTING, LOCK_IN_PROCESS, 0, NULL,
			      pfi->mtxmanager, DEBUG_LVL, pfi->msgFile);
	}

	/*
	 * STEP 5: compute the numeric factorization in parallel
	 */
	{
		ChvManager *chvmanager;
		int stats[20];
		int error;

		chvmanager = ChvManager_new();
		ChvManager_init(chvmanager, LOCK_IN_PROCESS, 1);
		IVfill(20, stats, 0);
		rootchv = FrontMtx_MT_factorInpMtx(pfi->frontmtx, mtxA, MAGIC_TAU, MAGIC_DTOL,
						   chvmanager, ownersIV, 0,
						   &error, pfi->cpus, stats, DEBUG_LVL,
						   pfi->msgFile);
		ChvManager_free(chvmanager);
		if (DEBUG_LVL > 1) {
			fprintf(msgFile, "\n\n factor matrix");
			FrontMtx_writeForHumanEye(pfi->frontmtx, pfi->msgFile);
			fflush(pfi->msgFile);
		}
		if (rootchv != NULL) {
			fprintf(pfi->msgFile, "\n\n matrix found to be singular\n");
			exit(-1);
		}
		if (error >= 0) {
			fprintf(pfi->msgFile, "\n\n fatal error at front %d", error);
			exit(-1);
		}
	}

	/*
	 * STEP 6: post-process the factorization
	 */
	ssolve_postfactor(pfi->frontmtx, pfi->msgFile);

	/*
	 * STEP 7: get the solve map object for the parallel solve
	 */
	{
		pfi->solvemap = SolveMap_new();
		SolveMap_ddMap(pfi->solvemap, *symmetryflagi4,
			       FrontMtx_upperBlockIVL(pfi->frontmtx),
			       FrontMtx_lowerBlockIVL(pfi->frontmtx), pfi->nthread, ownersIV,
			       FrontMtx_frontTree(pfi->frontmtx), RNDSEED, DEBUG_LVL,
			       pfi->msgFile);
	}

	/* cleanup: */
	InpMtx_free(mtxA);
	IVL_free(symbfacIVL);
	Graph_free(graph);
	IV_free(ownersIV);
}

DenseMtx *fsolve_MT(struct factorinfo *pfi, DenseMtx *mtxB)
{
	DenseMtx *mtxX;
	/*
	 * STEP 8: permute the right hand side into the new ordering
	 */
	ssolve_permuteB(mtxB, pfi->oldToNewIV, pfi->msgFile);


	/*
	 * STEP 9: solve the linear system in parallel
	 */
	{
		mtxX = DenseMtx_new();
		DenseMtx_init(mtxX, SPOOLES_REAL, 0, 0, pfi->size, 1, 1, pfi->size);
		DenseMtx_zero(mtxX);
		FrontMtx_MT_solve(pfi->frontmtx, mtxX, mtxB, pfi->mtxmanager,
					pfi->solvemap, pfi->cpus, DEBUG_LVL,
					pfi->msgFile);
		if (DEBUG_LVL > 1) {
			fprintf(pfi->msgFile, "\n\n solution matrix in new ordering");
			DenseMtx_writeForHumanEye(mtxX, pfi->msgFile);
			fflush(pfi->msgFile);
		}
	}

	/*
	 * STEP 10: permute the solution into the original ordering
	 */
	ssolve_permuteout(mtxX, pfi->newToOldIV, pfi->msgFile);

	/* Cleanup */
	DenseMtx_free(mtxB);

	return mtxX;
}

#endif

/** 
 * factor a system of the form (au - sigma * aub)
 * 
*/

FILE *msgFile;
struct factorinfo pfi;

void spooles_factor(double *ad, double *au,  double *adb, double *aub, 
             double *sigma,ITG *icol, ITG *irow, ITG *neq, ITG *nzs, 
	     ITG *symmetryflag, ITG *inputformat, ITG *nzs3)
{
	int size = *neq;
	int symmetryflagi4=*symmetryflag;
	InpMtx *mtxA;

	if(symmetryflagi4==0){
	    printf(" Factoring the system of equations using the symmetric spooles solver\n");
	}else if((symmetryflagi4==2)&&(*inputformat!=4)){
	    printf(" Factoring the system of equations using the unsymmetric spooles solver\n");
	}

/*	if(*neq==0) return;*/
 
	if ((msgFile = fopen("spooles.out", "a")) == NULL) {
		fprintf(stderr, "\n fatal error in spooles.c"
			"\n unable to open file spooles.out\n");
	}

	/*
	 * Create the InpMtx object from the CalculiX matrix
	 *      representation
	 */

	{
	    int col, ipoint, ipo;
	    int nent,i,j;
	    
	    mtxA = InpMtx_new();
	    
	    if((*inputformat==0)||(*inputformat==3)){
		nent = *nzs + *neq;	/* estimated # of nonzero entries */
	    }else if(*inputformat==1){
		nent=2**nzs+*neq;
	    }else if(*inputformat==2){
		nent=0;
		for(i=0;i<*neq;i++){
		    for(j=0;j<*neq;j++){
			if(fabs(ad[i*(int)*nzs+j])>1.e-20) nent++;
		    }
		}
	    }else if(*inputformat==4){
	      nent = *nzs;
	    }
	    
	    InpMtx_init(mtxA, INPMTX_BY_ROWS, SPOOLES_REAL, nent, size);
	    
	    /* inputformat:
	       0: sparse lower triangular matrix in ad (diagonal)
	          and au (lower triangle)
	       1: sparse full matrix in ad (diagonal)
	          and au (first lower triangle, then upper triangle; lower
	          and upper triangle have nonzero's at symmetric positions)
	       2: full matrix in field ad
	       3: sparse full matrix in ad (diagonal)
	          and au (off-diagonal positions)
               4: sparse full matrix in au, row after row
                  (CFD applications, call from compfluid.c)  */
	    
	    if(*inputformat==0){
		ipoint = 0;
		
		if(*sigma==0.){
		    for (col = 0; col < size; col++) {
//			printf("row=%d,col=%d,value=%e\n",col,col,ad[col]);
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
//			printf("row=%d,col=%d,value=%e\n",col,row,au[ipo]);
			    InpMtx_inputRealEntry(mtxA, col, row,
						  au[ipo]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
		else{
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]-*sigma*adb[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    InpMtx_inputRealEntry(mtxA, col, row,
						  au[ipo]-*sigma*aub[ipo]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
	    }else if(*inputformat==1){
		ipoint = 0;
		
		if(*sigma==0.){
		    for (col = 0; col < size; col++) {
//			printf("row=%d,col=%d,value=%e\n",col,col,ad[col]);
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
//			printf("row=%d,col=%d,value=%e\n",row,col,au[ipo]);
			    InpMtx_inputRealEntry(mtxA, row,col,
						  au[ipo]);
//			printf("row=%d,col=%d,value=%e\n",col,row,au[ipo+*nzs]);
			    InpMtx_inputRealEntry(mtxA, col,row,
						  au[ipo+(int)*nzs3]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
		else{
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]-*sigma*adb[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    InpMtx_inputRealEntry(mtxA, row,col,
						  au[ipo]-*sigma*aub[ipo]);
			    InpMtx_inputRealEntry(mtxA, col,row,
						  au[ipo+(int)*nzs3]-*sigma*aub[ipo+(int)*nzs3]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
	    }else if(*inputformat==2){
		for(i=0;i<*neq;i++){
		    for(j=0;j<*neq;j++){
			if(fabs(ad[i*(int)*nzs+j])>1.e-20){
			    InpMtx_inputRealEntry(mtxA,j,i,
						  ad[i*(int)*nzs+j]);
			}
		    }
		}
	    }else if(*inputformat==3){
		ipoint = 0;
		
		if(*sigma==0.){
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    /* next change on 27.06.2020 */
			    /*InpMtx_inputRealEntry(mtxA, col, row,
						    au[ipo]);*/
			    InpMtx_inputRealEntry(mtxA, row, col,
						  au[ipo]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
		else{
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]-*sigma*adb[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    /* next change on 27.06.2020 */
			    /* InpMtx_inputRealEntry(mtxA, col, row, 
			       au[ipo]-*sigma*aub[ipo]);*/
			    InpMtx_inputRealEntry(mtxA, row, col, 
						  au[ipo]-*sigma*aub[ipo]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
	    }else if(*inputformat==4){
	        int row;
	        for(row = 0; row < size; row++){
		    for(ipo = irow[row]; ipo < irow[row+1]; ipo++){
		        col = icol[ipo]-1;
			//			printf("spooles %d %d %d %f\n",ipo,row,col,au[ipo]);
		        InpMtx_inputRealEntry(mtxA, col, row, 
					au[ipo]);
		     }
	        }
		  
	    }		  
	    
	    InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);
            
	    if (DEBUG_LVL > 1) {
		fprintf(msgFile, "\n\n input matrix");
		InpMtx_writeForHumanEye(mtxA, msgFile);
		fflush(msgFile);
	    }
	}
	
	/* solve it! */


#ifdef USE_MT
	/* Rules for parallel solve:
           a. determining the maximum number of cpus:
              - if NUMBER_OF_CPUS>0 this is taken as the number of
                cpus in the system
              - else it is taken from _SC_NPROCESSORS_CONF, if strictly
                positive
              - else 1 cpu is assumed (default)
           b. determining the number of cpus to use
              - if CCX_NPROC_EQUATION_SOLVER>0 then use
                CCX_NPROC_EQUATION_SOLVER cpus
              - else if OMP_NUM_THREADS>0 use OMP_NUM_THREADS cpus
              - else use the maximum number of cpus
	 */
	if (num_cpus < 0) {
	    int sys_cpus;
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
	    
	    envloc = getenv("CCX_NPROC_EQUATION_SOLVER");
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
	    
	}
	printf(" Using up to %d cpu(s) for spooles.\n\n", num_cpus);
	if (num_cpus > 1) {
	    /* do not use the multithreaded solver unless
	     * we have multiple threads - avoid the
		 * locking overhead
		 */
		factor_MT(&pfi, mtxA, size, msgFile,&symmetryflagi4);
	} else {
		factor(&pfi, mtxA, size, msgFile,&symmetryflagi4);
	}
#else
	printf(" Using 1 cpu for spooles.\n\n");
	factor(&pfi, mtxA, size, msgFile,&symmetryflagi4);
#endif
}

/** 
 * solve a system of equations with rhs b
 * factorization must have been performed before 
 * using spooles_factor
 * 
*/

void spooles_solve(double *b, ITG *neq)
{
	/* rhs vector B
	 * Note that there is only one rhs vector, thus
	 * a bit simpler that the AllInOne example
	 */
	int size = *neq;
	DenseMtx *mtxB,*mtxX;

//	printf(" Solving the system of equations using the symmetric spooles solver\n");

	{
		int i;
		mtxB = DenseMtx_new();
		DenseMtx_init(mtxB, SPOOLES_REAL, 0, 0, size, 1, 1, size);
		DenseMtx_zero(mtxB);
		for (i = 0; i < size; i++) {
			DenseMtx_setRealEntry(mtxB, i, 0, b[i]);
		}
		if (DEBUG_LVL > 1) {
			fprintf(msgFile, "\n\n rhs matrix in original ordering");
			DenseMtx_writeForHumanEye(mtxB, msgFile);
			fflush(msgFile);
		}
	}

#ifdef USE_MT
	//	printf(" Using up to %d cpu(s) for spooles.\n\n", num_cpus);
	if (num_cpus > 1) {
		/* do not use the multithreaded solver unless
		 * we have multiple threads - avoid the
		 * locking overhead
		 */
		mtxX=fsolve_MT(&pfi, mtxB);
	} else {
		mtxX=fsolve(&pfi, mtxB);
	}
#else
	//	printf(" Using 1 cpu for spooles.\n\n");
	mtxX=fsolve(&pfi, mtxB);
#endif

	/* convert the result back to Calculix representation */
	{
		int i;
		for (i = 0; i < size; i++) {
			b[i] = DenseMtx_entries(mtxX)[i];
		}
	}
	/* cleanup */
	DenseMtx_free(mtxX);
}

void spooles_cleanup()
{
        
  	FrontMtx_free(pfi.frontmtx);
	IV_free(pfi.newToOldIV);
	IV_free(pfi.oldToNewIV);
	SubMtxManager_free(pfi.mtxmanager);
	if (pfi.solvemap)
		SolveMap_free(pfi.solvemap);
	ETree_free(pfi.frontETree);
	fclose(msgFile);
}


/** 
 * factor a system of the form (au - sigma * aub)
 * 
*/

FILE *msgFilf;
struct factorinfo pfj;

void spooles_factor_rad(double *ad, double *au,  double *adb, double *aub, 
			double *sigma,ITG *icol, ITG *irow,
			ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
			ITG *iexpl)
{
        int symmetryflagi4=*symmetryflag;
	int size = *neq;
	InpMtx *mtxA;

	/* messages only for implicit calculations (for explicit calculations
	   the output should be kept to a minimum because of performance
	   reasons*/
	
	if(*iexpl<=1){
	  if(symmetryflagi4==0){
	    printf(" Factoring the system of equations using the symmetric spooles solver\n");
	  }else if((symmetryflagi4==2)&&(*inputformat!=4)){
	    printf(" Factoring the system of radiation equations using the unsymmetric spooles solver\n");
	  }
	}

/*	if(*neq==0) return;*/
 
	if ((msgFilf = fopen("spooles.out", "a")) == NULL) {
		fprintf(stderr, "\n fatal error in spooles.c"
			"\n unable to open file spooles.out\n");
	}

	/*
	 * Create the InpMtx object from the Calculix matrix
	 *      representation
	 */

	{
	    int col, ipoint, ipo;
	    int nent,i,j;
	    
	    mtxA = InpMtx_new();
	    
	    if((*inputformat==0)||(*inputformat==3)){
		nent = *nzs + *neq;	/* estimated # of nonzero entries */
	    }else if(*inputformat==1){
		nent=2**nzs+*neq;
	    }else if(*inputformat==2){
		nent=0;
		for(i=0;i<*neq;i++){
		    for(j=0;j<*neq;j++){
			if(fabs(ad[i*(int)*nzs+j])>1.e-20) nent++;
		    }
		}
	    }else if(*inputformat==4){
	        nent = *nzs;
	    }
	    
	    InpMtx_init(mtxA, INPMTX_BY_ROWS, SPOOLES_REAL, nent, size);
	    
	    /* inputformat:
	       0: sparse lower triangular matrix in ad (diagonal)
	          and au (lower triangle)
	       1: sparse lower + upper triangular matrix in ad (diagonal)
	          and au (first lower triangle, then upper triangle; lower
	          and upper triangle have nonzero's at symmetric positions)
	       2: full matrix in field ad
	       3: sparse upper triangular matrix in ad (diagonal)
	          and au (upper triangle)  */
	    
	    if(*inputformat==0){
		ipoint = 0;
		
		if(*sigma==0.){
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    InpMtx_inputRealEntry(mtxA, col, row,
						  au[ipo]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
		else{
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]-*sigma*adb[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    InpMtx_inputRealEntry(mtxA, col, row,
						  au[ipo]-*sigma*aub[ipo]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
	    }else if(*inputformat==1){
		ipoint = 0;
		
		if(*sigma==0.){
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    InpMtx_inputRealEntry(mtxA, row,col,
						  au[ipo]);
			    InpMtx_inputRealEntry(mtxA, col,row,
						  au[ipo+(int)*nzs]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
		else{
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]-*sigma*adb[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    InpMtx_inputRealEntry(mtxA, row,col,
						  au[ipo]-*sigma*aub[ipo]);
			    InpMtx_inputRealEntry(mtxA, col,row,
						  au[ipo+(int)*nzs]-*sigma*aub[ipo+(int)*nzs]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
	    }else if(*inputformat==2){
		for(i=0;i<*neq;i++){
		    for(j=0;j<*neq;j++){
			if(fabs(ad[i*(int)*nzs+j])>1.e-20){
			    InpMtx_inputRealEntry(mtxA,j,i,
						  ad[i*(int)*nzs+j]);
			}
		    }
		}
	    }else if(*inputformat==3){
		ipoint = 0;
		
		if(*sigma==0.){
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    InpMtx_inputRealEntry(mtxA, row, col,
						  au[ipo]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
		else{
		    for (col = 0; col < size; col++) {
			InpMtx_inputRealEntry(mtxA, col, col, ad[col]-*sigma*adb[col]);
			for (ipo = ipoint; ipo < ipoint + icol[col]; ipo++) {
			    int row = irow[ipo] - 1;
			    InpMtx_inputRealEntry(mtxA, row, col, 
						  au[ipo]-*sigma*aub[ipo]);
			}
			ipoint = ipoint + icol[col];
		    }
		}
	    }else if(*inputformat==4){
	        int row;
	        for(row = 0; row < size; row++){
		    for(ipo = irow[row]; ipo < irow[row+1]; ipo++){
		        col = icol[ipo] - 1;
		        InpMtx_inputRealEntry(mtxA, col, row, 
					au[ipo]);
		     }
	        }
		  
	    }		  
	    
	    InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);
	    
	    if (DEBUG_LVL > 1) {
		fprintf(msgFilf, "\n\n input matrix");
		InpMtx_writeForHumanEye(mtxA, msgFilf);
		fflush(msgFilf);
	    }
	}

	/* solve it! */


#ifdef USE_MT
	/* Rules for parallel solve:
           a. determining the maximum number of cpus:
              - if NUMBER_OF_CPUS>0 this is taken as the number of
                cpus in the system
              - else it is taken from _SC_NPROCESSORS_CONF, if strictly
                positive
              - else 1 cpu is assumed (default)
           b. determining the number of cpus to use
              - if CCX_NPROC_EQUATION_SOLVER>0 then use
                CCX_NPROC_EQUATION_SOLVER cpus
              - else if OMP_NUM_THREADS>0 use OMP_NUM_THREADS cpus
              - else use the maximum number of cpus
	 */
	if (num_cpus < 0) {
	    int sys_cpus;
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
	    
	    envloc = getenv("CCX_NPROC_EQUATION_SOLVER");
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
	    
	}
	if(*iexpl<=1){
	  printf(" Using up to %d cpu(s) for spooles.\n\n", num_cpus);
	}
	if (num_cpus > 1) {
		/* do not use the multithreaded solver unless
		 * we have multiple threads - avoid the
		 * locking overhead
		 */
		factor_MT(&pfj, mtxA, size, msgFilf,&symmetryflagi4);
	} else {
		factor(&pfj, mtxA, size, msgFilf,&symmetryflagi4);
	}
#else
	if(*iexpl<=1){
	  printf(" Using 1 cpu for spooles.\n\n");
	}
	factor(&pfj, mtxA, size, msgFilf,&symmetryflagi4);
#endif
}

/** 
 * solve a system of equations with rhs b
 * factorization must have been performed before 
 * using spooles_factor
 * 
*/

void spooles_solve_rad(double *b, ITG *neq)
{
	/* rhs vector B
	 * Note that there is only one rhs vector, thus
	 * a bit simpler that the AllInOne example
	 */
	int size = *neq;
	DenseMtx *mtxB,*mtxX;

	//	printf(" Solving the system of radiation equations using the unsymmetric spooles solver\n");

	{
		int i;
		mtxB = DenseMtx_new();
		DenseMtx_init(mtxB, SPOOLES_REAL, 0, 0, size, 1, 1, size);
		DenseMtx_zero(mtxB);
		for (i = 0; i < size; i++) {
			DenseMtx_setRealEntry(mtxB, i, 0, b[i]);
		}
		if (DEBUG_LVL > 1) {
			fprintf(msgFilf, "\n\n rhs matrix in original ordering");
			DenseMtx_writeForHumanEye(mtxB, msgFilf);
			fflush(msgFilf);
		}
	}

#ifdef USE_MT
	//	printf(" Using up to %d cpu(s) for spooles.\n\n", num_cpus);
	if (num_cpus > 1) {
		/* do not use the multithreaded solver unless
		 * we have multiple threads - avoid the
		 * locking overhead
		 */
		mtxX=fsolve_MT(&pfj, mtxB);
	} else {
		mtxX=fsolve(&pfj, mtxB);
	}
#else
	//	printf(" Using 1 cpu for spooles.\n\n");
	mtxX=fsolve(&pfj, mtxB);
#endif

	/* convert the result back to Calculix representation */
	{
		int i;
		for (i = 0; i < size; i++) {
			b[i] = DenseMtx_entries(mtxX)[i];
		}
	}
	/* cleanup */
	DenseMtx_free(mtxX);
}

void spooles_cleanup_rad()
{
        
  	FrontMtx_free(pfj.frontmtx);
	IV_free(pfj.newToOldIV);
	IV_free(pfj.oldToNewIV);
	SubMtxManager_free(pfj.mtxmanager);
	if (pfj.solvemap)
		SolveMap_free(pfj.solvemap);
	ETree_free(pfj.frontETree);
	fclose(msgFilf);
}


/** 
 * spooles: Main interface between Calculix and spooles:
 *
 * Perform 3 operations:
 * 1) factor
 * 2) solve
 * 3) cleanup
 *   
 */

void spooles(double *ad, double *au, double *adb, double *aub, double *sigma,
             double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs, 
             ITG *symmetryflag, ITG *inputformat, ITG *nzs3)
{

    if(*neq==0){
	return;
    }else if(*neq==1){
	b[0]/=ad[0];
	return;
    }

    /*    FORTRAN(spooles_write,(ad,au,adb,aub,sigma,b,icol,irow,neq,nzs,
	  symmetryflag,inputformat,nzs3));*/
    
  spooles_factor(ad,au,adb,aub,sigma,icol,irow,neq,nzs,symmetryflag,
		 inputformat,nzs3);
  
  spooles_solve(b,neq);
  
  spooles_cleanup();

}

#endif
