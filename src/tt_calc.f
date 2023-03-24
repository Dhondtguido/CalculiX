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
      subroutine tt_calc(xflow,Tt,pt,kappa,r,A,Ts,icase)
!     
!     this subroutine solves the implicit equation
!     f=xflow*dsqrt(Tt)/(a*Pt)-C*(TtdT)**expon*(Ttdt-1)**0.5d0
!     in order to find Tt when Ts , xflow, pt and a are given
!
!     author: Yannick Muller
!     
      implicit none
!
      integer icase,i
!     
      real*8 xflow,Tt,pt,Ts,kappa,r,f,df,A,expon,Tt_old,C,TtzTs,
     &     deltaTt,TtzTs_crit,Qred,h1,h2,h3
!     
      expon=-0.5d0*(kappa+1.d0)/(kappa-1.d0)
!     
      C=dsqrt(2.d0/r*kappa/(kappa-1.d0))
!     
!     f=xflow*dsqrt(Tt)/(A*pt)-C*(Tt/Ts)**expon*(Tt/Ts-1)**0.5d0
!
!     df=-C*Ttdt**expon*(expon/Ts*(Tt/Ts-1)**0.5d0
!     &     -0.5d0*Tt/Ts/Ts*(Tt/Ts-1.d0)**(-0.5d0))
!     
      if(dabs(xflow).le.1e-9) then
         Tt=Ts
         return
      endif
!
!     initial guess
!
      Qred=dabs(xflow)*dsqrt(Ts)/(A*pt)
      Tt=Ts*(1+(Qred**2/C**2))
!     
!     adiabatic
!     
      if(icase.eq.0) then
!
         TtzTs_crit=(kappa+1.d0)/2.d0
!         
!     isothermal
!
      else
!     
         TtzTs_crit=(1d0+(kappa-1.d0)/(2.d0*kappa))
!
      endif
!
      if(Tt/Ts.gt.TtzTs_crit) then
         Tt=Ts*(TtzTs_crit+1.d0)/2.d0
      endif
!
      i=0
      Tt_old=Tt
!
      do 
         i=i+1
         Ttzts=Tt/Ts
         h1=Ttzts-1.d0
         h2=dsqrt(h1)
         h3=Ttzts**expon
!     
         f=C*h2*h3
!     
         df=0.5d0*dabs(xflow)/(A*pt*dsqrt(Tt))
     &        -C*h2*h3*(expon/Tt+1.d0/(2.d0*h1*Ts))
!
         Qred=abs(xflow)*dsqrt(Tt)/(A*pt)
!
         f=Qred-f
         deltaTt=-f/df
!     
         Tt=Tt+deltaTt
!     
         if((((dabs(Tt-Tt_old)/Tt_old).le.1.d-8))
     &        .or.((dabs(Tt-Tt_old)).le.1.d-10) 
     &        .or.((dabs(f).le.1d-5).and.(deltaTt.lt.1d-3))) then
!     
            if(Tt/Ts.gt.TtzTs_crit) then
               Tt=Ts*TtzTs_crit
            endif
            exit
         else if((i.gt.40)) then
            Tt=0.99d0*Ts*TtzTs_crit
            exit
         endif
         Tt_old=Tt
      enddo
!     
      return
      end
      
      
      
