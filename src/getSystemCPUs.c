/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License, or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */
/*
	Added for a consistent way to determine the number of available cores
	for both *NIX and Windows systems. Convergent Mechanical, Dec 5, 2012
*/

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

#ifdef __WIN32
#include <windows.h>

ITG getSystemCPUs(){
	SYSTEM_INFO sysinfo;
	GetSystemInfo( &sysinfo );
	return sysinfo.dwNumberOfProcessors;
}

#else
#include <unistd.h>

ITG getSystemCPUs(){
	return sysconf(_SC_NPROCESSORS_CONF);;
}
#endif
