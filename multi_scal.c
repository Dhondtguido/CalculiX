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
 *  Performs the scalar product of the mth column of au_1 and the nth column of au_2
 *
 *  [in] au_1 		values of matrix 1
 *  [in] irow_1		rows of matrix 1
 *  [in] jq_1		column pointer to iorw_1
 *  [in] au_2 		values of matrix 2
 *  [in] irow_2		rows of matrix 2
 *  [in] jq_2		column pointer to iorw_2
 *  [in] m			row number
 *  [in] n			column number
 *  [out] value		result
 *  [in] flag		not used
**/

void multi_scal(double *au_1,ITG * irow_1,ITG * jq_1,
		double *au_2,ITG * irow_2,ITG * jq_2,
		ITG m,ITG n,double*value,ITG *flag){
  
  /*Performs the scalar product of the mth column of au_1 and 
    the nth column of au_2*/
    
  ITG pt1,pt2;
  double val=0.0;
  
  pt1=jq_1[m]-1;
  pt2=jq_2[n]-1; 
  while((pt1<jq_1[m+1]-1)&&(pt2<jq_2[n+1]-1)){
    if ((pt1<jq_1[m+1]-1)&&(pt2<jq_2[n+1]-1)){ //normal case
      if (irow_1[pt1]==irow_2[pt2]){
	val+=au_1[pt1]*au_2[pt2];
	pt1++;
	pt2++;
	*flag=1;	
      }else{
	if (irow_1[pt1]<irow_2[pt2]) {
	  pt1++;
	}else{
	  pt2++;
	}
      }
    }
  }	   	
  
  *value=val;
  return;
}
