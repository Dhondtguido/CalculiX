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
      subroutine hnorm(xflow,rho,b,theta,dg,s0,friction,xks,he)
!
!     determine the normal depth for a trapezoidal channel cross
!     section:
!
!     INPUT:
!
!     xflow:     fluid mass flow (may be negative)
!     rho:       fluid density
!     b:         width at bottom
!     theta:     angle describing the width increase with fluid depth
!     dg:        earth acceleration
!     s0:        change of bottom height with channel length
!     friction:  friction coefficient (if White-Coolebrook)
!     xks:       if > 0: White-Coolebrook law is to be taken      
!                if <= 0: Manning law is to be taken; dabs(xks) is 
!                         the Manning coefficient
!     
!     OUTPUT:
!     
!     he:        normal depth
!
      implicit none
!      
      real*8 xflow,rho,b,dg,s0,he,theta,tth,friction,xks,henew,
     &     third,threetenth,fourtenth,cthi,c1
!
!     bottom is steepening: he is infinite
!
      if(s0.le.0.d0) then
        he=1.d30
        return
!
!     no flow
!
      elseif(dabs(xflow).lt.1.d-20) then
        he=0.d0
        return
      endif
!
      third=1.d0/3.d0
!
      tth=dtan(theta)
      cthi=1.d0/dcos(theta)
!
      if(xks.gt.0.d0) then
!
!       White-Coolebrook; initial value: 
!
        c1=friction*(xflow/rho)**2/(8.d0*dg*s0)
!
!       triangular cross section
!
        if(b.lt.1.d-20) then
          he=(c1*2.d0*cthi/tth**3)**(1.d0/5.d0)
          return
        endif
!        
!       initial value        
!        
        he=(c1/(b*b))**third
!        
        do
          henew=(c1*(b+2.d0*he*cthi))**third/(b+he*tth)
          if(dabs(henew-he).lt.1.d-3*he) exit
          he=henew
        enddo
      else
!
!       Manning; initial value: 
!
        threetenth=3.d0/10.d0
        fourtenth=4.d0/10.d0
!     
        c1=((xks*xflow/rho)**2/s0)**threetenth
!
!       triangular cross section
!
        if(b.lt.1.d-20) then
          he=(c1*(2d0*cthi)**fourtenth/tth)**(1.d0/1.6d0)
          return
        endif
!        
!       initial value        
!        
        he=c1*b**(fourtenth-1.d0)
!        
        do
          henew=c1*(b+2d0*he*cthi)**fourtenth/(b+he*tth)
          if(dabs(henew-he).lt.1.d-3*he) exit
          he=henew
        enddo
!        
      endif
!
      return
      end
      

