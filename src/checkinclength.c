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

void checkinclength(double *time,double *ttime,double *theta, double *dtheta,
          ITG *idrct, double *tper,double *tmax, double *tmin, double *ctrl, 
          double *amta,ITG *namta, ITG *itpamp, ITG *inext, double *dthetaref, 
          ITG *itp,ITG *jprint, ITG *jout){

    ITG id,istart,iend,inew,ireduceincrement;
    double reftime;

    ITG i0,ir,ip,ic,il,ig,ia;
    double df,dc,db,dd,ran,can,rap,ea,cae,ral,da;
    i0=ctrl[0];ir=ctrl[1];ip=ctrl[2];ic=ctrl[3];il=ctrl[4];ig=ctrl[5];ia=ctrl[7];
    df=ctrl[10];dc=ctrl[11];db=ctrl[12];da=ctrl[13];dd=ctrl[16];
    ran=ctrl[18];can=ctrl[19];rap=ctrl[22];
    ea=ctrl[23];cae=ctrl[24];ral=ctrl[25];

    /* check whether the new increment size is not too big */
    
    if(*dtheta>*tmax){
	*dtheta=*tmax;
//	printf(" the increment size exceeds thetamax and is decreased to %e\n\n",*dtheta**tper);
    } 
    
    /* if itp=1 the increment just finished ends at a time point */
    
    if((*itpamp>0)&&(*idrct==0)){
	if(namta[3**itpamp-1]<0){
	    reftime=*ttime+*time+(*dtheta)**tper;
	}else{
	    reftime=*time+(*dtheta)**tper;
	}
	istart=namta[3**itpamp-3];
	iend=namta[3**itpamp-2];
	FORTRAN(identamta,(amta,&reftime,&istart,&iend,&id));
	if(id<istart){
	    inew=istart;
	}else{
	    inew=id+1;
	}
//	printf("istart=%" ITGFORMAT ",iend=%" ITGFORMAT ",inext=%" ITGFORMAT ",inew=%" ITGFORMAT "\n",istart,iend,*inext,inew);

            /* inew: smallest time point exceeding time+dtheta*tper
               inext: smallest time point exceeding time */

	/* the check with *tmin
           was introduced to circumvent the following problem: if the new data point is
           smaller than the next data point, but the distance is very small, the next *dheta
           calculated a few lines below may be zero due to the subtraction of two nearly equal
           numbers; a zero *dtheta leads to a fatal error */

	ireduceincrement=0;
	if(*inext<iend){
	    if(fabs((amta[2**inext-2]-reftime)/(*tper))<*tmin){
		ireduceincrement=1;
	    }
	}

	if((*inext<inew)||(ireduceincrement==1)){
//	if((*inext<inew)||(fabs((amta[2**inext-2]-reftime)/(*tper))<*tmin)){
	  //	if((*inext<inew)||(fabs((amta[2**inext-2]-reftime))<1.e-10)){
	  
	    if(namta[3**itpamp-1]<0){
		*dtheta=(amta[2**inext-2]-*ttime-*time)/(*tper);
	    }else{
		*dtheta=(amta[2**inext-2]-*time)/(*tper);
	    }
	    (*inext)++;
	    *itp=1;
//	    printf(" the increment size exceeds a time point and is decreased to %e\n\n",*dtheta**tper);
	}else{*itp=0;}
    }

    /* check whether the step length is not exceeded */

    if(*dtheta>1.-*theta){
	*dtheta=1.-*theta;
	*dthetaref=*dtheta;
	printf(" the increment size exceeds the remainder of the step and is decreased to %e\n\n",*dtheta**tper);
//	if(*dtheta<=1.e-6){(*ttime)+=(*dtheta**tper);}
    }
    
    return;
}
