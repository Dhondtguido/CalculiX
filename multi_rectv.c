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
#include <string.h> 
#include "CalculiX.h"
#include "mortar.h"
/**  
 *  sparse matrix-vector multiplication \f$ v_r=a_1^Tb \f$
 * 
 *  [in] au_1 		values of matrix 1
 *  [in] irow_1		rows of matrix 1
 *  [in] jq_1		column pointer to iorw_1
 *  [in]  n_1		number of rows matrix 1  
 *  [in]  m_1		number of columns matrix 1
 *  [in]  b 		vector b
 *  [out] v_rp		result vector
**/

void multi_rectv(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1, ITG m_1,
		 double * b, double ** v_rp){ 
  
  /*Result fields*/
  
  ITG i,j;
  double *v_r=NULL;
  
  
  NNEW(v_r,double,m_1);   
  for(j=0;j<m_1;j++){
    for (i=jq_1[j]-1;i<jq_1[j+1]-1;i++){
      v_r[j]+=au_1[i]*b[irow_1[i]-1];
    }
  }
  
  *v_rp=v_r;
  
  return;
}
