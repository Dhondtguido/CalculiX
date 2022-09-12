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
      subroutine fsub(time,t,a,b,dd,h1,h2,h3,h4,func,funcp)
!
      implicit none
!
      real*8 time,t,a,b,dd,h1,h2,h3,h4,fexp,fsin,fcos,func,funcp,
     &  h8,h9,h10,h11,h12,h13
!
      fexp=dexp(-h1*t)
      fsin=dsin(dd*t)
      fcos=dcos(dd*t)
      h8=(a+b*time)*fexp/h2
      h9=-b*fexp/h2
      h10=-h8*h1
      h11=h8*dd
      h12=h9*(-h1*t-h3/h2)
      h13=h9*(dd*t+h4)
!
!     function
!
c      fsub=(a+b*time)*fexp*(-h1*fsin-dd*fcos)/h2-b*fexp/h2*((-h1*t-h3/h2)*
c     &     fsin-(dd*t+h4)*fcos)
      func=h10*fsin-h11*fcos+h12*fsin-h13*fcos
!
!     derivative of the function
!
      funcp=-h1*func+dd*(h10*fcos+h11*fsin+h12*fcos+h13*fsin)
!
      return
      end
