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

ITG strcmp1(const char *s1, const char *s2)
{
  ITG a,b;

  do {
    a=*s1++;
    b=*s2++;

/* the statement if((a=='\0')||(b=='\0')) has been treated separately
   in order to avoid the first field (s1) to be defined one longer
   than required; s1 is assumed to be a variable field, s2 is
   assumed to be a fixed string */

    if(b=='\0'){
      a='\0';
      b='\0';
      break;
    }
    if(a=='\0'){
      a='\0';
      b='\0';
      break;
    }
  }while(a==b);
  return(a-b);
}
	  
