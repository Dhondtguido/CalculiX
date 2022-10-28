!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
!     Solve the Bresse equation for the turbulent stationary flow
!     in channels with a non-erosive bottom
!     
      subroutine calcdhds(xflow,b,tth,cthi,s0,sqrts0,friction,xks,h,
     &     dg,rho,dhds)
!
      implicit none
!
      real*8 xflow,b,tth,cthi,s0,sqrts0,friction,xks,h,dg,area,bb,p,sf,
     &     rho,dhds
!     
      area=h*(b+h*tth)
      bb=b+2.d0*h*tth
      p=b+2.d0*h*cthi
!
      if(xks.gt.0.d0) then
!
!       White-Colebrook
!
        sf=friction*p*(xflow/rho)**2/(8.d0*dg*area**3)
      else
!
!       Manning
!
        sf=(xks*xflow/rho)**2*p**(4.d0/3.d0)/(area**(10.d0/3.d0))
      endif
!
      dhds=(s0-sf)/(sqrts0-(xflow/rho)**2*bb/(dg*area**3))
!
      return
      end
      
