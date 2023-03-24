/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "CalculiX.h"

void elementcpuload(ITG *neapar,ITG *nebpar,ITG *ne,ITG *ipkon,ITG *num_cpus){

 /*  divides the elements into ranges with an equal number of
     active elements (element numbering may have gaps) for
     parallel processing on different cpus */

    ITG i,nepar,*ipar=NULL,idelta,isum;
    
    NNEW(ipar,ITG,*ne);

    nepar=0;
    for(i=0;i<*ne;i++){
	if(ipkon[i]>-1){

	    /* active element */
	    
	    ipar[nepar]=i;
	    nepar++;
	}
    }
    if(nepar<*num_cpus) *num_cpus=nepar;

    /* dividing the element number range into num_cpus equal numbers of 
       active elements */

    idelta=nepar/(*num_cpus);
    isum=0;
    for(i=0;i<*num_cpus;i++){
	neapar[i]=ipar[isum];
	if(i!=*num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=nepar;
	}
	nebpar[i]=ipar[isum-1];
    }
    
    SFREE(ipar);
    
    return;
}
