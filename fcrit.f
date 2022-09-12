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
      subroutine fcrit(time,t,a,b,ze,d,dd,h1,h2,h3,h4,func,funcp)
!
      implicit none
!
      real*8 time,t,a,b,ze,d,dd,h1,h2,h3,h4,fexp,func,funcp
!
      fexp=dexp(-h1*t)
!
!     function
!
      func=((a+b*time)*(-t*h2-h3)-b*(-t*t*h2-2.d0*t*h3-2.d0*h4))*fexp
!
!     derivative of the function
!
      funcp=((a+b*time)*t-b*(h3+t*h2+t*t))*fexp
!
      return
      end
