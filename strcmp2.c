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

ITG strcmp2(const char *s1, const char *s2, ITG length)
{

/* comparison of the first "length" characters unless s1
   and/or s2 has less characters */

  ITG a,b,i;

  i=0;
  do {
    a=*s1++;
    b=*s2++;

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
    i++;
  }while((a==b)&&(i<length));
  return(a-b);
}
	  
