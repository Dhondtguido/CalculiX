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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <locale.h>
#include "CalculiX.h"
#include "mortar.h"

void printenergy(ITG *iexpl,double *ttime,double *theta,double *tper,
		 double *energy,ITG *ne,ITG *nslavs,double *ener,
		 double *energyref,double *allwk,double *dampwk,
		 double *ea,double *energym,double *energymold,ITG *jnz,
		 ITG *mscalmethod,ITG *mortar,ITG *mi){
  
  ITG i;

  double denergymax;

  setlocale(LC_NUMERIC, "C");

  if(*iexpl>1){
    printf(" actual total time=%e\n\n",*ttime+*theta**tper);
    if(*mortar==-1){
      energy[3]=0.;
      for(i=*ne;i<*ne+*nslavs;i++){
	energy[3]+=ener[2*mi[0]*i+1];
      }
    }
  }
	
  printf(" initial energy (at start of step) = %e\n\n",*energyref);

  printf(" since start of the step: \n");
  printf(" external work = %e\n",*allwk);
  printf(" work performed by the damping forces = %e\n",*dampwk);
  printf(" netto work = %e\n\n",*allwk+*dampwk);

  printf(" actual energy: \n");
  printf(" internal energy = %e\n",energy[0]);
  printf(" kinetic energy = %e\n",energy[1]);
  printf(" elastic contact energy = %e\n",energy[2]);
  printf(" energy lost due to friction = %e\n",energy[3]);
  printf(" total energy  = %e\n\n",energy[0]+energy[1]+energy[2]+energy[3]);

  printf(" energy increase = %e\n\n",energy[0]+energy[1]+energy[2]
	 +energy[3]-*energyref);

  printf(" energy balance (absolute) = %e \n",energy[0]+energy[1]
	 +energy[2]+energy[3]-*energyref-*allwk-*dampwk);

  /* Belytschko criterion */

  denergymax=energy[0];
  if(denergymax<energy[1]){
    denergymax=energy[1];}
  if(denergymax<fabs(*allwk)) denergymax=fabs(*allwk);

  if(denergymax>*ea**energym){
    *energym=(*energymold**jnz+denergymax)/(*jnz+1);}
  else {
    *energym=*energymold;}
  *energymold=*energym;   

  if(*energym>1.e-30){
    printf(" energy balance (relative) = %f %% \n\n",
	   fabs((energy[0]+energy[1]+energy[2]+energy[3]-*energyref-*allwk
		 -*dampwk)/(*energym)*100.));
  }else{
    printf(" energy balance (relative) =0 %% \n\n");
  }
	
  /*Energy balance to evaluate mass scaling*/
      
  if((*mscalmethod==1)||(*mscalmethod==3)){
    printf(" artificial energy due to selective mass scaling = %e\n",
	   energy[4]);
	    
    printf(" energy balance with mass scaling(relative) = %f %% \n\n",
	   fabs((energy[0]+energy[1]+energy[2]+energy[3]+energy[4]
		 -*energyref-*allwk-*dampwk)/(*energym)*100.));
	    
  }
  
return;
}
