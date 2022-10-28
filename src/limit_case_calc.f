!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2005 Guido Dhondt
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
      subroutine limit_case_calc(a2,pt1,Tt2,xflow,zeta,r,kappa,
     &     pt2_lim,M2)
!     
!     For restrictor elements A1<A2
!     if A1 is critical, pt2 can be as low as possible without any effect on the flow
!     it is necessary to compute pt2_lim which satisfies the element equation
!
!     xflow*dsqrt(R*Tt1)/(A1*Pt1*dsqrt(kappa)-dsqrt(2/(kappa-1)*(Pt1/Pt2)
!     **((kappa-1)-1))/(zeta*kappa)))/(Pt1/Pt2)**((kappa+1)/(2*zeta*kappa))=0
!
!     **Changed 15.11.2007**
!     since we are in the case A1<=A2 
!     Pt1/Pt2=(1+0.5*(k-1)M1**2)**(zeta*kappa/(kappa-1))
!     A1 is critical M1=1 which simplifies the previous equation
!     (Pt1/Pt2)_crit=(1+0.5*(k-1))**(zeta*kappa/(kappa-1))
!     It is then possible to determine pt2_lim knowing zeta and Pt1
!
!     Once Pt2_lim is calculated it is possible in turn to calculate
!     the corresponding Mach number M2 satisfying the flow equation
!     for section A2 in terms of M2
!     xflow*dsqrt(R*Tt2)/(A2*Pt2*dsqrt(kappa))-M2/(1+(kappa-1)/2*M2**2)
!     **(0.5*(kappa+1)/(kappa-1))
!
!     author: Yannick Muller
!
      implicit none
!
      integer i
!
      real*8 pt1,Tt2,xflow,zeta,r,kappa,pt2_lim,M2,Qred,expon1,
     &     expon2,root,a2,f,df,pt2,km1,kp1,pt1pt2
!
!     **changed 15.11.2006**
      pt2=0.99d0*pt1
      pt1pt2=pt1/pt2
      km1=kappa-1.d0
      kp1=kappa+1.d0
      expon1=-0.5d0*(kp1)/(zeta*kappa)
      expon2=(km1)/(zeta*kappa)
      root=2.d0/(km1)*((pt1pt2)**expon2-1.d0)
      Qred=dsqrt(kappa/R)*(pt1pt2)**(-0.5d0*kp1/(kappa*zeta))
     &     *dsqrt(2.d0/km1*((pt1pt2)**(km1/(kappa*zeta))-1.d0))
!
      pt2_lim=pt1/(1.d0+0.5d0*(km1))**(zeta*kappa/(km1))

!
!     M2_lim calculation
!
      M2=0.5d0
      expon1=-0.5d0*(kp1)/(km1)
      Qred=dabs(xflow)*dsqrt(R*Tt2)/(A2*pt2_lim*dsqrt(kappa))
      if(Qred.gt.((1.d0+0.5d0*(km1))**expon1)) then
         Qred=(1.d0+0.5d0*(km1))**expon1
      endif
!
      do 
         root=(1.d0+0.5d0*(km1)*M2**2.d0)
         f=Qred-M2*root**(expon1)
!
         df=root**expon1*(-1d0+0.5d0*(kp1)*M2**2.d0*root**(-1.d0))
!
         if(dabs(-f/df).le.1d-6) then
            M2=M2-f/df
            
            exit
         endif
!
         M2=M2-f/df
      enddo
!
      return
      end
         
