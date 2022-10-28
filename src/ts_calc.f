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
      subroutine ts_calc(xflow,Tt,pt,kappa,r,A,Ts,icase)
!     
!     Calculation of the static temperature Ts from the total 
!     temperature Tt, the total pressure Pt and the mass flow xflow
!
!     Use of the General Gas Equation
!
!     this subroutine solves the implicit equation
!     xflow*dsqrt(Tt)/(A*Pt)-C*(Tt/Ts)**expon*(Tt/Ts-1)**0.5d0=0.d0
!
!          expon=-0.5d0*(kappa+1.d0)/(kappa-1.d0)
!          C=dsqrt(2.d0/r*kappa/(kappa-1.d0))
!
!     xflow is the mass flow for 360 degrees, sign is irrelevant
!     for this routine
!
!     author: Yannick Muller
!     
      implicit none
!
      integer icase,i,rf
!     
      real*8 xflow,Tt,pt,Ts,kappa,r,f,df,A,expon,Ts_old,C,TtzTs,
     &     deltaTs,TtzTs_crit, Qred_crit,Qred,h1,h2,h3,Ts_min,
     &     Ts_max,f_min,f_max
!
!
!
!     regula falsi flag (rf=0 newton-raphson is active, 
!     rf=1 regula falsi is active)
      rf=0
!
      expon=-0.5d0*(kappa+1.d0)/(kappa-1.d0)
!     
      C=dsqrt(2.d0/r*kappa/(kappa-1.d0))
!     
!     f=xflow*dsqrt(Tt)/(A*pt)-C*(Tt/Ts)**expon*(Tt/Ts-1)**0.5d0
!
!     df=-C*Tt/Ts**expon*(expon/Ts*(Tt/Ts-1)**0.5d0
!     &     -0.5d0*Tt/Ts/Ts*(Tt/Ts-1.d0)**(-0.5d0))
!     
      Ts_old=Tt
!    
      if(dabs(xflow).le.1d-10) then
         Ts=Tt
         return
      endif
!
      Qred=abs(xflow)*dsqrt(Tt)/(A*pt)
!
!     optimised estimate of T static
!
      Ts=Tt/(1.d0+(Qred**2.d0/C**2.d0))   
!     
!     adiabatic
!     
      if(icase.eq.0) then
         TtzTs_crit=(kappa+1.d0)/2.d0
!         
!     isothermal
!
      else
         TtzTs_crit=(1d0+(kappa-1.d0)/(2.d0*kappa))
      endif
!
      Qred_crit=C*(TtzTs_crit)**expon*(Ttzts_crit-1.d0)**0.5d0
!     
      if(Qred.ge.Qred_crit) then
         Ts=Tt/TtzTs_crit
         return
      endif
!  
      i=0
!     
!     start of the Newton-Raphson-Procedure to solve the nonlinear
!     equation
!
      do
!     
         if((Ts.ge.Tt).or.(Ts.le.Tt/TtzTs_crit))then
!
!           to avoid estimated values, in that case the slower 
!           regula falsi is used
!
            rf=1
            exit
         endif
!     
         i=i+1
         Ttzts=Tt/Ts
         h1=Ttzts-1.d0
         h2=dsqrt(h1)
         h3=Ttzts**expon
!     
         f=C*h2*h3
!     
         df=f*(expon+0.5d0*Ttzts/h1)/Ts
!
         f=Qred-f
         deltaTs=-f/df
!     
         Ts=Ts+deltaTs
!     
         if((((dabs(Ts-Ts_old)/ts_old).le.1.d-8))
     &        .or.((dabs(Ts-Ts_old)).le.1.d-10)) then
            exit
         else if(i.gt.20) then
            rf=1
            exit
         endif
         Ts_old=Ts
      enddo
!     
      if(rf.eq.0)then
         return
      endif
!     
!     regula falsi procedure
!
      i=1
      Ts_min=Tt/TtzTs_crit
      Ts_max=Tt
      Ts_old=Tt+1.d0
!
      f_min=Qred-C*dsqrt(Tt/Ts_min-1)*(Tt/Ts_min)**expon 
      f_max=Qred
!     
      do
         i=i+1
         Ts=(Ts_min+Ts_max)/2.d0
         f=Qred-C*dsqrt(Tt/Ts-1)*(Tt/Ts)**expon
!     
         if((((dabs(Ts-Ts_old)/Ts_old).le.1.d-8))
     &        .or.((dabs(Ts-Ts_old)).le.1.d-10)) then
            exit 
         endif
!     
         if(i.gt.10000000)then
            write(*,*)'*ERROR in ts_calc.f'
            write(*,*)'       max. iteration number exceeded'
           call exit(201)
         endif
!     
         if(f_min*f.le.0.d0)then
            Ts_max=Ts
            f_max=f
         else
            Ts_min=Ts
            f_min=f
         endif 
!     
         Ts_old=Ts
! 
      enddo  
!
      return     
      end
      
      
      
