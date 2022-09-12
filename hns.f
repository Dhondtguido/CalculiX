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
      subroutine hns(xflow,rho,b,theta,dg,sqrts0,h1,h2)
!
!     determine the flow depth h2 downstream of a hydraulic jump, 
!     corresponding to a upstream flow depth of h1
!     
!
!     INPUT:
!
!     xflow:     fluid mass flow (may be negative)
!     rho:       fluid density
!     b:         width at bottom
!     theta:     angle describing the width increase with fluid depth
!     dg:        earth acceleration
!     sqrts0:    sqrt(1.d0-s0*s0), where s0 is the change of bottom
!                height with channel length
!     h1:        depth upstream of jump
!     
!     OUTPUT:
!     
!     h2:        depth downstream of jump
!
      implicit none
!      
      real*8 b,rho,dg,sqrts0,xflow,h1,h2,hk,ygA1,A1,zeta1,
     &     A2,ygA2,h2new,theta,tth
!
      call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
!
      h2=h1*(-1.d0+dsqrt(1.d0+8.d0*(hk/h1)**3))/2.d0
!
      if(dabs(theta).lt.1.d-10) return
!
!     hns for a trapezoid, non-rectangular cross section
!
      tth=dtan(theta)
      ygA1=h1*h1*(b/2.d0+h1*tth/3.d0)
      A1=h1*(b+h1*tth)
      zeta1=(xflow/rho)**2/(dg*A1)+sqrts0*ygA1
!
      do
        A2=h2*(b+h2*tth)
        ygA2=h2*h2*(b/2.d0+h2*tth/3.d0)
        h2new=h2*dsqrt((zeta1-(xflow/rho)**2/(dg*A2))/
     &       (sqrts0*ygA2))
        if(dabs(h2new-h2).lt.1.d-3*h2) exit
        h2=h2new
      enddo
!
      write(*,*) 'hns ','h1= ',h1,'h2= ',h2,'hk= ',hk
!
      return
      end
      

