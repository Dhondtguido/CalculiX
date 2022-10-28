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
      subroutine pk_y0_yg(p2p1,beta,kappa,y0,yg)
!
      implicit none
!
      real*8 p2p1,beta,kappa,y0,yg,pcrit
!
!     adiabatic expansion factor y0 measured (eq.15-17) 
!
!     author: Yannick Muller
!     
      pcrit=(2.d0/(kappa+1.d0))**(kappa/(kappa-1.d0))

      if(p2p1.ge.0.63d0) then
         y0=1d0-(0.41d0+0.35d0*beta**4.d0)/kappa*(1.d0-p2p1)
      else
         y0=1d0-(0.41d0+0.35d0*beta**4.d0)/kappa*(1.d0-0.63d0)
     &        -(0.3475d0+0.1207d0*beta**2.d0-0.3177d0*beta**4.d0)
     &        *(0.63d0-p2p1)
!         
      endif
!    
!     adiabatic expension factor yg isentropic eq 18
!
      if(p2p1.ge.1d0) then
         yg=1.d0
!         
      elseif (p2p1.ge.pcrit) then
         yg=p2p1**(1.d0/kappa)*dsqrt(kappa/(kappa-1.d0)
     &        *(1.d0-p2p1**((kappa-1.d0)/kappa)))/dsqrt(1.d0-p2p1)
!      
      else
!     critical pressure ratio
         yg=(2.d0/(kappa+1.d0))**(1.d0/(kappa-1.d0))
     &       *dsqrt(kappa/(kappa+1.d0))/dsqrt(1.d0-p2p1)
      endif
!     
      return
!      
      end
