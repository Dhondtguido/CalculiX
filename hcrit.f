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
      subroutine hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
!
!     determine the critical depth for a trapezoidal channel cross
!     section:
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
!     
!     OUTPUT:
!     
!     hk:        critical depth
!
      implicit none
!      
      real*8 xflow,rho,b,dg,sqrts0,hk,theta,tth,c1,third,hknew
!
      if(dabs(xflow).lt.1.d-20) then
        hk=0.d0
        return
      elseif(b.lt.1.d-20) then
!
!       triangular cross section
!
        hk=(2*xflow*xflow/(dg*sqrts0*tth*tth))**(1.d0/5.d0)
        return
      endif
!
      third=1.d0/3.d0
!
      hk=((xflow/(rho*b))**2/(dg*sqrts0))**third
!
      if(dabs(theta).lt.1.d-10) return
!
!     critical depth for trapezoid, non-rectangular cross section
!
      tth=dtan(theta)
      c1=((xflow/rho)**2/(dg*sqrts0))**third
!     
      do
        hknew=c1*(b+2.d0*hk*tth)**third/(b+hk*tth)
        if(dabs(hknew-hk).lt.1.d-3*hk) exit
        hk=hknew
      enddo
!
      return
      end
      

