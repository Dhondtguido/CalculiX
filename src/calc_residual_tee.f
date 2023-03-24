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
      real*8 function calc_residual_tee(pt1,Tt1,xflow1,xflow2,
     &pt2,Tt2,A1,A2,zeta_fac,kappa,R,ider,iflag,zeta)
!
!     Function for calculating the residual of both branches of a tee
!
!     author: Yannick Muller
!
      implicit none
!
      integer 
!     Isothermal/adiabatic
     &icase,
!     Residual/Derivatives
     &ider,
!     Critical conditions at inlet or outlet
     &icrit1,
     &icrit2,
!     Where we are in the program
     &iflag
!
      real*8 
!     Residual
     &f,
!     Kappa stuff
     &kappa,
     &R,
!     State variables
     &pt1,
     &pt2,
     &Tt1,
     &Tt2,
     &Ts0,
     &Ts1,
     &Ts2,
     &xflow1,
     &xflow2,
!
     &pt2_lim,
!
     &zeta,
     &zeta_fac,
!
!     Areas
     &A1,
     &A2,
!     Reduced mass flows
     &Q0,
     &Q1,
     &Q2,
     &Q_crit,
!     Pressure ratios
     &pspt_crit,
     &pspt0,
     &pspt1,
     &pspt2,
!     Flow velocities
     &w1,
     &w2,
     &w1w2,
     &w2w1,
!    Mach numbers,
     &M1,
     &M2
!
      icrit1 = 0
      icrit2 = 0
!     setting icase (always adiabatic)
      icase=0;
!
!     Critical values
      Q_crit = dsqrt(kappa/R)*
     &   (1+0.5d0*(kappa-1))**(-0.5d0*(kappa+1)/(kappa-1))
      pspt_crit = (2/(kappa+1)) ** (kappa/(kappa-1))
!
!     These reduced mass flows are equivalent
!     (reduced mass flow at inlet)
      Q0 = xflow1*dsqrt(Tt1)/pt1/A1
      Q1 = xflow2*dsqrt(Tt1)/pt1/A2
      if(Q1.ge.Q_crit) then
         Q1 = Q_crit
         icrit1 = 1
         write(*,*)'*WARNING in Tee:'
         write(*,*)'Critical conditions at 1'
      endif
!     Reduced mass flow at outlet
      Q2 = xflow2*dsqrt(Tt1)/pt2/A2
      if(Q2.ge.Q_crit) then
         Q2 = Q_crit
         icrit2 = 1
         write(*,*)'*WARNING in Tee:'
         write(*,*)'Critical conditions at 2'
      endif
!
!     Flow velocity at inlet
!     Static temperature
      Ts0=Tt1
      call ts_calc(xflow1,Tt1,pt1,kappa,r,A1,Ts0,icase)
!     Pressure ratio
      pspt0 = (Ts0/Tt1)**(kappa/(kappa-1))
!     Velocity
      call wpi(w1, pspt0, Q0, 
     &      dsqrt(Tt1),kappa,R) 
!
!     Flow velocity at outlet
!     Static temperature
      call ts_calc(xflow2,Tt1,pt1,kappa,r,A2,Ts1,icase)
!     Pressure ratio
      pspt1 = (Ts1/Tt1)**(kappa/(kappa-1))
!     Velocity      
      call wpi(w2, pspt1, Q1, 
     &      dsqrt(Tt1),kappa,R) 
!
!     Velocity ratio
c      w2w1=w2/w1
c      w1w2=w1/w2
      if(w2.eq.0.d0)then
         w1w2=1d30
      else
         w1w2=w1/w2
      endif
      if(w1.eq.0.d0)then
         w2w1=1d30
      else
         w2w1=w2/w1
      endif
!
!     Zeta calculation
      zeta=1.d0+0.3d0*W2W1**2
      zeta=zeta*(W1W2)**2
!
      zeta = zeta_fac*zeta
!
!     Residual calculation
      if(icrit2.ne.1) then
         if(icrit1.ne.1) then
            f = pt2 - pt1*pspt1**zeta
         else
            f = xflow2*dsqrt(Tt1)/pt1/A2-Q_crit
         endif
      else
         f = xflow2*dsqrt(Tt1)/pt2/A2-Q_crit
      endif
!
      if(iflag.eq.3) then
!
         write(1,57)'             zeta= ',zeta
 57      format(1x,a,f9.4)
!
      else if (iflag.eq.4) then

!        Calculate Mach numbers
         call machpi(M1,pspt0,kappa,R)
         call ts_calc(xflow2,Tt2,pt2,kappa,r,A2,Ts2,icase)
!        Pressure ratio
         pspt2 = (Ts2/Tt2)**(kappa/(kappa-1))
         call machpi(M2,pspt2,kappa,R)
      
!         write(1,80)'Inlet: Tt1= ',Tt1,
!     &              ', pt1= ',pt1,', M1= ',M1
!     
!         write(1,77)'mass flow = ',xflow2,', kappa = ',kappa,
!     &              ', zeta= ',zeta
!
!         write(1,80)'Outlet: Tt2= ',Tt2,
!     &              ', pt2= ',pt2,', M2= ',M2
!     
! 80   format(3x,a,f10.6,a,f10.2,a,f10.6)
! 77   format(3x,a,f10.6,a,f10.2,a,f10.6)
      endif
!
      calc_residual_tee=f
!
      return
      end
