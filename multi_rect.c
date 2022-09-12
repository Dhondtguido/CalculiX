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
 *  sparse matrix multiplication a_r=a_1^T*a_2
 * 
 *  [in] au_1 		values of matrix 1
 *  [in] irow_1		rows of matrix 1
 *  [in] jq_1		column pointer to iorw_1
 *  [in]  n_1		number of rows matrix 1  
 *  [in]  m_1		number of columns matrix 1
  *  [in] au_2 		values of matrix 2
 *  [in] irow_2		rows of matrix 2
 *  [in] jq_2		column pointer to iorw_2
 *  [in]  n_2		number of rows matrix 2  
 *  [in]  m_2		number of columns matrix 2
 *  [out] au_rp		values of result matrix
 *  [out] irow_rp		rows of result matrix
 *  [out] jq_r		column pointer to iorw_r
 *  [out] nzs		number of non-zero entries in au_r
**/

void multi_rect(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1, ITG m_1,
		double *au_2,ITG * irow_2,ITG * jq_2,ITG n_2, ITG m_2,
		double **au_rp,ITG **irow_rp,ITG * jq_r,ITG *nzs){ 
  
  /*Result fields*/
  
  ITG *irow=NULL,ifree=1,i,j,l,m,flag=0;
  
  double *autmp=NULL,value;
  
  /*Perform a_1T*a_2*/
  
  if(n_1!=n_2) {
    printf(" *ERROR in mutli_rec : Matrix sizes are not compatible\n");
    return;
  } 
  
  irow=*irow_rp;
  autmp=*au_rp; 
  
  jq_r[0]=1;
  for(j=0;j<m_2;j++){
    m=j+1;
    if(jq_2[j+1]-jq_2[j]>0){
      for (i=0;i<m_1;i++){
	l=i+1;
	flag=0;
	value=0.0;
	multi_scal(au_1,irow_1,jq_1,au_2,irow_2,jq_2,i,j,&value,&flag);
	if(flag!=0) insertas_ws(&irow,&l,&m,&ifree,nzs,&value,&autmp);
      }
    }
    jq_r[m]=ifree;
  }
  
  /* Sort the column and compute jq*/
  
  *nzs=ifree-1;
  RENEW(autmp,double,*nzs);
  RENEW(irow,ITG,*nzs);
  *irow_rp=irow;*au_rp=autmp;   

  return;
}
