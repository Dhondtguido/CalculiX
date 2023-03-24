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
 *Calculate C=A+B for sparse matrices A,B
 * 
 *  [in] au_1		sparse matrix A
 *  [in] irow_1           row numbers for matrix A
 *  [in] jq_1             colum numbers for matrix A
 *  [in] n_1              # rows  for matrix A 
 *  [in] m_1              # columns  for matrix A
 *  [in] au_2		sparse matrix B
 *  [in] irow_2           row numbers for matrix B
 *  [in] jq_2             colum numbers for matrix B
 *  [in] n_2              # rows  for matrix B
 *  [in] m_2              # columns  for matrix B
 *  [out] au_rp		sparse matrix C
 *  [out] irow_rp          row numbers for matrix C
 *  [out] jq_r            colum numbers for matrix C
 *  [out] nzs             # nonzero entries in C 
 *
*/

void add_rect(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1, ITG m_1,
	      double *au_2,ITG * irow_2,ITG * jq_2,ITG n_2, ITG m_2,
	      double **au_rp,ITG **irow_rp,ITG * jq_r,ITG *nzs){
  
  /*Result fields*/
    
  ITG *irow=NULL,ifree=1,j,m,pt1,pt2,row1,row2;
  
  double *autmp=NULL,value;
  
  if((m_1!=m_2)||(n_1!=n_2)){
    printf(" *ERROR in mutli_rec : Matrix sizes are not compatible\n");
    return;
  } 
  
  irow=*irow_rp;
  autmp=*au_rp; 
  
  jq_r[0]=1;
  
  for(j=0;j<m_2;j++){
    pt1=jq_1[j]-1;
    pt2=jq_2[j]-1;
    m=j+1; 
    while((pt1<jq_1[m]-1)||(pt2<jq_2[m]-1)){
      if ((pt1<jq_1[m]-1)&&(pt2<jq_2[m]-1)){ //normal case
	row1=irow_1[pt1];
	row2=irow_2[pt2];
	if (row1==row2){
	  value=au_1[pt1]+au_2[pt2];
	  insertas_ws(&irow,&row1,&m,&ifree,nzs,&value,&autmp);
	  pt1++;
	  pt2++;
	}else{
	  if (row1<row2) {
	    insertas_ws(&irow,&row1,&m,&ifree,nzs,&au_1[pt1],&autmp);
	    pt1++;
	  }else{
	    insertas_ws(&irow,&row2,&m,&ifree,nzs,&au_2[pt2],&autmp);
	    pt2++;
	  }
	}
      }
      else{ 
	if (pt1<jq_1[m]-1){ //column 2 finished
	  insertas_ws(&irow,&irow_1[pt1],&m,&ifree,nzs,&au_1[pt1],&autmp);
	  pt1++;				
	}
	else{ //column 1 finished
	  insertas_ws(&irow,&irow_2[pt2],&m,&ifree,nzs,&au_2[pt2],&autmp);
	  pt2++;						
	}
      }
    }	   
    jq_r[m]=ifree;
  }
  
  *nzs=ifree-1;
  RENEW(irow,ITG,*nzs);
  RENEW(autmp,double,*nzs);
  *irow_rp=irow;*au_rp=autmp;
  
  return;
}
