!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine fsuper(time,t,a,b,h1,h2,h3,h4,h5,h6,func,funcp)
!
      implicit none
!
      real*8 time,t,a,b,h1,h2,h3,h4,h5,h6,fexm,fexp,func,funcp
!
      fexm=dexp(h1*t)
      fexp=dexp(-h2*t)
!
!     function
!
      func=(a+b*time)*(fexm*h3+fexp*h4)
     &    -b*(fexm*(t*h3-h5)+fexp*(t*h4+h6))
!
!     derivative of the function
!
      funcp=(a+b*time)*(fexm-fexp)-b*(fexm*(t-h3)-fexp*(t+h4))

!
      return
      end
