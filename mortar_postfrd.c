/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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
#include <time.h>
#include "CalculiX.h"
#include "mortar.h"
/** 
*  function undoing preparations for mortar output in frd-file
* Author: Saskia Sitzmann
 *  [in,out] ne		number of elements
 *  [in] nslavs		number of slave nodes
 *  [in] mi		(1) max # of integration points per element (2) max degree of freedom per element
 *  [in] nk		number of nodes 
 *  [in,out] nkon		size of kon
 *  [in,out] fn		internal forces
 *  [in] cfs		contat forces
 *  [in] cfm		not used any more
**/
void mortar_postfrd(ITG *ne,ITG *nslavs, ITG *mi, ITG *nk, ITG *nkon,
		   double *fn, double *cfs, double *cfm){
  
  ITG i,k,l,mt=mi[1]+1;
  
  /* substract slave nodes from ne and nkon again */    
  *ne-=*nslavs;
  *nkon-=*nslavs; 
  /* substract contact forces from reaction rorces Fn */ 
  for(i=0;i<mt**nk;i++){     
    fn[i]=fn[i]-cfs[i];
  }
  
  return;
}
