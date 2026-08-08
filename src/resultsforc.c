/*     CalculiX - A 3-dimensional finite element program                 */
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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

void resultsforc(ITG *nk,double *f,double *fn,ITG *nactdof,ITG *ipompc,
		 ITG *nodempc,double *coefmpc,char *labmpc,ITG *nmpc,
		 ITG *mi,double *fmpc,ITG *calcul_fn,ITG *calcul_f,
                 ITG *num_cpus,ITG *iponoel){

  ITG i,j,ist,node,ndir,mt=mi[1]+1,index,index2;

  double forcempc;
    
  /*     subtracting the mpc force (for each linear mpc there is one
	 force; the actual force in a node belonging to the mpc is
	 obtained by multiplying this force with the nodal coefficient.
	 The force has to be subtracted from f, since it does not
	 appear on the rhs of the equation system */
    
  if(*calcul_fn==1){
    for(i=0;i<*nmpc;i++){
      ist=ipompc[i]-1;
      node=nodempc[3*ist]-1;
      ndir=nodempc[3*ist+1];
      if(ndir>3) continue;
      forcempc=fn[mt*node+ndir]/coefmpc[ist];
      fmpc[i]=forcempc;
      fn[mt*node+ndir]=0.;
      index=nodempc[3*ist+2]-1;
      if(index==-1) continue;
      do{
	node=nodempc[3*index]-1;
	ndir=nodempc[3*index+1];
	fn[mt*node+ndir]-=coefmpc[index]*forcempc;
	index=nodempc[3*index+2]-1;
	if(index==-1) break;
      }while(1);
    }
  }
    
  //     calculating the system force vector

  if(*calcul_f==1){
    forparll(&mt,nactdof,f,fn,nk,num_cpus);
  }
    
  /* adding the mpc force again to fn */

  if(*calcul_fn==1){
    for(i=0;i<*nmpc;i++){
      ist=ipompc[i]-1;
      node=nodempc[3*ist]-1;
      ndir=nodempc[3*ist+1];
      if(ndir>3) continue;
      forcempc=fmpc[i];
      if(iponoel[node]!=0){
	fn[mt*node+ndir]=forcempc*coefmpc[ist];
      }
      index=nodempc[3*ist+2]-1;
      if(index==-1) continue;

      do{
	node=nodempc[3*index]-1;
	ndir=nodempc[3*index+1];
	if(iponoel[node]!=0){
	  fn[mt*node+ndir]+=coefmpc[index]*forcempc;
	}
	index=nodempc[3*index+2]-1;
	if(index==-1) break;

      }while(1);
    }
  }
    
  return;

}
