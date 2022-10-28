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
#include <string.h>
#include <stdlib.h>
#include "CalculiX.h"

void crackpropdata(char *jobnamec,ITG *nelcon,double *elcon,double **crconp,
		   ITG *ncrconst,ITG *ncrtem,ITG *imat,char *matname,
		   ITG *ntmat_,ITG *ncmat_,char **paramp,ITG *nparam,ITG *law){

  /* user routine to read the crack propagation data


      INPUT:

      jobnamec           name of the input deck (without .inp)
      nelcon(1,i)        number of temperature data points for material i
      nelcon(2,i)        material code number for material i
      elcon(0..ncmat_,1..ntmat_,i)
                         material constants for material i
      imat               crack material number
      matname(i)         name of material i
      ntmat_             max. number of material temperature data points
      ncmat_             max. number of material constants


      OUTPUT (required)

      crcon(0..ncrconst,1..nrctem)
                         crack propagation constants
      ncrconst           number of crack propagation constants per 
                         temperature
      ncrtem             number of temperature data points for the crack 
                         material


      OUTPUT (optional; is subsequently available for further use in
              routine crackrate.f)

      param(1...nparam)  parameters returned for use in the crack propagation
                         law (character*132)
      nparam             number of parameters (integer)
      law                number of the crack propagation law (integer) */

  
  ITG i,j,index;

  double *crcon=NULL;

  crcon=*crconp;

  *ncrconst=-nelcon[(*imat-1)*2]-100;
  *ncrtem=nelcon[(*imat-1)*2+1];

  index=(*ncmat_+1)**ntmat_*(*imat-1);

  NNEW(crcon,double,(*ncrconst+1)**ncrtem);

  for(j=0;j<*ncrtem;j++){
    for(i=0;i<=*ncrconst;i++){
      crcon[j*(*ncrconst+1)+i]=elcon[index+j*(*ncmat_+1)+i];
    }
  }

  *crconp=crcon;
  
  return;
}
