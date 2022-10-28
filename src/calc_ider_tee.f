!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
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
!     author: Yannick Muller
!
      subroutine calc_ider_tee(df,pt1,Tt1,xflow1,xflow2,pt2,
     &Tt2,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)
!
      implicit none
!
      integer ider,iflag
!
      real*8 
     &df(6),
     &pt1,
     &pt2,
     &Tt1,
     &Tt2,
     &xflow1,
     &xflow2,
     &A1,
     &A2,
     &kappa,
     &R,
     &calc_residual_tee,
!     Accuracy of the numerical derivatives
     &eps,
!     Step for the derivative
     &h,
     &f0,
     &zeta_fac,
     &zeta
!
      eps = 1.0e-8
!
      f0 = calc_residual_tee(pt1,Tt1,xflow1,xflow2,pt2,
     &Tt2,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)
!
      h = eps*dabs(pt1)
      if(h.eq.0)then
         h = eps
      endif
      df(1) = (calc_residual_tee(pt1+h,Tt1,xflow1,xflow2,pt2,
     &Tt2,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)-f0)/h
!
      h = eps*dabs(Tt1)
      if(h.eq.0)then
         h = eps
      endif
      df(2) = (calc_residual_tee(pt1,Tt1+h,xflow1,xflow2,pt2,
     &Tt2,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)-f0)/h
!
      h = eps*dabs(xflow1)
      if(h.eq.0)then
         h = eps
      endif
      df(3) = (calc_residual_tee(pt1,Tt1,xflow1+h,xflow2,pt2,
     &Tt2,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)-f0)/h
!
      h = eps*dabs(xflow2)
      if(h.eq.0)then
         h = eps
      endif
      df(4) = (calc_residual_tee(pt1,Tt1,xflow1,xflow2+h,pt2,
     &Tt2,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)-f0)/h
!
      h = eps*dabs(pt2)
      if(h.eq.0)then
         h = eps
      endif
      df(5) = (calc_residual_tee(pt1,Tt1,xflow1,xflow2,pt2+h,
     &Tt2,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)-f0)/h
!
      h = eps*dabs(Tt2)
      if(h.eq.0)then
         h = eps
      endif
      df(6) = (calc_residual_tee(pt1,Tt1,xflow1,xflow2,pt2,
     &Tt2+h,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)-f0)/h
!
      return
      end
