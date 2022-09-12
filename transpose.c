/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

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
#include <time.h>
#include <string.h> 
#include "CalculiX.h"
#include "mortar.h"
/**
 *  function to transpose sparse matrix

 * Author: Saskia Sitzmann

 * [in] au 	matrix A values
 * [in] jq	matrix A column pointer to irow
 * [in] irow	matrix A rows
 * [in] dim	dimension of matrix A (dim x dim)
 * [out] au_t 	matrix A^T values
 * [out] jq_t	matrix A^T column pointer to irow
 * [out] irow_t	matrix A^T rows
 **/

void transpose(double *au,ITG *jq,ITG *irow,ITG *dim,
	       double *au_t,ITG *jq_t,ITG *irow_t){
  
  ITG  i,j,nzs_t=jq[*dim]-1,*mast1=NULL;
  
  /*transpose*/	
  
  NNEW(mast1,ITG,nzs_t);
  
  for(j=0;j<nzs_t;j++){au_t[j]=au[j];}
  
  /* mast1 contains the rows of the original matrix,
     irowbt the columns */
  
  for(j=0;j<*dim;j++){
    for(i=jq[j]-1;i<jq[j+1]-1;i++){
      mast1[i]=irow[i];
      irow_t[i]=j+1;
    }
  }  
  
  matrixsort(au_t,mast1,irow_t,jq_t,&nzs_t,dim); 
  SFREE(mast1);
  return;
}
