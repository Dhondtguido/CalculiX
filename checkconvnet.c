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
#include "CalculiX.h"
#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif

void checkconvnet(ITG *icutb, ITG *iin,
		  double *cam1t, double *cam1f, double *cam1p,
		  double *cam2t, double *cam2f, double *cam2p,
		  double *camt, double *camf, double *camp,
		  ITG *icntrl, double *dtheta, double *ctrl,
                  double *cam1a,double *cam2a,double *cama,
                  double *vamt, double *vamf, double *vamp, double *vama,
                  double *qa, double *qamt, double *qamf,
                  double *ramt, double *ramf, double *ramp, ITG *iplausi,
                  ITG *iaxial){
  
  ITG i0,ir,ip,ic,il,ig,ia,idivergence,iin_dyn,dyna_flag_1,dyna_flag_2;
  
  double c2t,c2f,c2p,c2a,c1t,c1f,c1p,a2t,a2f,a2p,a2a,a1t,a1f,a1p,qamp=1.,
         df,dc,db,dd,ran,can,rap,ea,cae,ral;

  i0=ctrl[0];ir=ctrl[1];ip=ctrl[2];ic=ctrl[3];il=ctrl[4];ig=ctrl[5];ia=ctrl[7];
  df=ctrl[10];dc=ctrl[11];db=ctrl[12];dd=ctrl[16];ran=ctrl[18];can=ctrl[19];
  rap=ctrl[22];ea=ctrl[23];cae=ctrl[24];ral=ctrl[25];c1t=ctrl[32];c1f=ctrl[33];
  c1p=ctrl[34];c2t=ctrl[35];c2f=ctrl[36];c2p=ctrl[37];c2a=ctrl[38];
  a1t=ctrl[40];a1f=ctrl[41];a1p=ctrl[42];a2t=ctrl[43];a2f=ctrl[44],a2p=ctrl[45];
  a2a=ctrl[46];
  
  /* checks for dynamic convergence
     - to avoid oscillations the check is done after a special number
       inn_dyn of iterations */

  iin_dyn=50;
  if(*iin%iin_dyn==0){
  
  /* criteria 1: change of sign */
     dyna_flag_1=0;
     if(*camt**cam1t>=0 && *camt**cam2t>=0){dyna_flag_1=dyna_flag_1;}
     else{dyna_flag_1=dyna_flag_1+1;}
     if(*camf**cam1f>=0 && *camf**cam2f>=0){dyna_flag_1=dyna_flag_1;}
     else{dyna_flag_1=dyna_flag_1+1;}
     if(*camp**cam1p>=0 && *camp**cam2p>=0){dyna_flag_1=dyna_flag_1;}
     else{dyna_flag_1=dyna_flag_1+1;}
     if(*cama**cam1a>=0 && *cama**cam2a>=0){dyna_flag_1=dyna_flag_1;}
     else{dyna_flag_1=dyna_flag_1+1;}  
  
  /* criteria 2: progression check */
     dyna_flag_2=0;
     if(fabs(*camt)<=fabs(*cam1t) && fabs(*cam1t)<=fabs(*cam2t)){dyna_flag_2=dyna_flag_2;}
     else{dyna_flag_2=dyna_flag_2+1;}
     if(fabs(*camf)<=fabs(*cam1f) && fabs(*cam1f)<=fabs(*cam2f)){dyna_flag_2=dyna_flag_2;}
     else{dyna_flag_2=dyna_flag_2+1;}  
     if(fabs(*camp)<=fabs(*cam1p) && fabs(*cam1p)<=fabs(*cam2p)){dyna_flag_2=dyna_flag_2;}
     else{dyna_flag_2=dyna_flag_2+1;}  
     if(fabs(*cama)<=fabs(*cam1a) && fabs(*cam1a)<=fabs(*cam2a)){dyna_flag_2=dyna_flag_2;}
     else{dyna_flag_2=dyna_flag_2+1;}    
  
  }
  
  
/* Das wird aus meiner Sicht nicht mehr benoetigt 
  if(*cam1t<*cam2t) {*cam2t=*cam1t;}
  if(*cam1f<*cam2f) {*cam2f=*cam1f;}
  if(*cam1p<*cam2p) {*cam2p=*cam1p;}
  if(*cam1a<*cam2a) {*cam2a=*cam1a;}
  */
  
    
  /* check for convergence or divergence; 
     the convergence check consists of 
     - a comparison of the correction in
       the latest network iteration with the change since the 
       start of the network calculations 
     - a comparison of the residual in the latest network
       iteration with mean typical values of the equation terms */

  if((fabs(*camt)<=c2t**vamt)&&(*ramt<c1t**qamt)&&(fabs(*camt)<=a2t)&&
     (*ramt<a1t/(*iaxial))&&(fabs(*camf)<=c2f**vamf)&&(*ramf<c1f**qamf)&&
     (fabs(*camf)<=a2f/(*iaxial))&&(*ramf<a1f/(*iaxial))&&
     (fabs(*camp)<=c2p**vamp)&&(*ramp<c1p*qamp)&&(fabs(*camp)<a2p)&&
     (*ramp<a1p)&&
     (fabs(*cama)<=c2a**vama)&&(fabs(*cama)<=a2a)&&(*iplausi==1)&&
     (*iin>3)){
      
      /* increment convergence reached */
      
      printf("      flow network: convergence in gas iteration %" ITGFORMAT " \n\n",*iin);
      *icntrl=1;
      *icutb=0;
  }
  
  else {

      idivergence=0;

      /* divergence based on temperatures */
      
      if((*iin>=20*i0)||(fabs(*camt)>1.e20)){
	  if((fabs(*cam1t)>=fabs(*cam2t))&&(fabs(*camt)>=fabs(*cam2t))&&(fabs(*camt)>c2t**vamt)){
	      idivergence=1;
	  }
      }

      /* divergence based on the mass flux */
      
      if((*iin>=20*i0)||(fabs(*camf)>1.e20)){
	  if((fabs(*cam1f)>=fabs(*cam2f))&&(fabs(*camf)>=fabs(*cam2f))&&(fabs(*camf)>c2f**vamf)){
	      idivergence=1;
	  }
      }

      /* divergence based on pressures */
      
      if((*iin>=20*i0)||(fabs(*camp)>1.e20)){
	  if((fabs(*cam1p)>=fabs(*cam2p))&&(fabs(*camp)>=fabs(*cam2p))&&(fabs(*camp)>c2p**vamp)){
	      idivergence=1;
	  }
      }

      /* divergence based on geometry */
      
      if((*iin>=20*i0)||(fabs(*cama)>1.e20)){
	  if((fabs(*cam1a)>=fabs(*cam2a))&&(fabs(*cama)>=fabs(*cam2a))&&(fabs(*cama)>c2a**vama)){
	      idivergence=1;
	  }
      }

      /* divergence based on the number of iterations */

      if(*iin>20*ic) idivergence=1;

      /* divergence based on singular matrix or negative pressures */

      if(*iin==0) idivergence=1;
      
      if(idivergence==1){
	  *dtheta=*dtheta*df;
	  printf("\n network divergence; the under-relaxation parameter is decreased to %e\n",*dtheta);
	  printf(" the network iteration for the increment is reattempted\n\n");
	  *iin=0;
	  (*icutb)++;
	  if(*icutb>ia){
	      qa[2]=0.25;
	      *icntrl=1;
//	    printf("\n *ERROR: too many cutbacks\n");
//	    FORTRAN(stop,());
	  }
      }else{
	  if(*iin%iin_dyn==0){
	     if(dyna_flag_1==0 && dyna_flag_2==0 && *iplausi==1){
		printf("      good convergence --> *dtheta is increased %" ITGFORMAT "\n",*iin);
		*dtheta=*dtheta*(1.2);
		if(*dtheta>=1){
	           *dtheta=1.;
		}	        
	     }
	     else if(dyna_flag_1!=0 && dyna_flag_2!=0 && *iplausi!=1){
		printf("      bad convergence progression --> *dtheta is decreased %" ITGFORMAT "\n",*iin);
		*dtheta=*dtheta*(0.8);	        
	     }
	     else{
	        printf("      no convergence\n\n"); 
	     }
	  }
	  else{
	     printf("      no convergence\n\n"); 
	  }
      }
  }
  return;
}
