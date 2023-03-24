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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

ITG strcpy2(char *s1, const char *s2, ITG length)
{

  /* copies up to 'length' characters from s2 into s1;
     copying stops after '\0' was copied */

  char b;
  ITG i,blank=0;
  
  for(i=0;i<length;i++) {
      if(blank==0){
	  b=*s2;
	  if(b=='\0')blank=1;
      }
      *s1=*s2;s2++;
      if(blank==1) return 0;
      s1++;
  }
  return 0;
}
	  




