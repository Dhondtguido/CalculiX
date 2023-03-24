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
!     this subroutine solves iteratively the following equation
!     to determine the pressure for which section A2 is critical
!
!     author: Yannick Muller
!
      subroutine pt2_lim_calc (pt1,a2,a1,kappa,zeta,pt2_lim)
!
      implicit none
!
      integer i 
!
      real*8 pt1,a2,a1,kappa,pt2_lim,x,zeta,f,df,expon1,
     &     expon2,expon3,cte,a2a1,kp1,km1,delta_x,fact1,fact2,term
!
      x=0.999d0
!
!     x belongs to interval [0;1]
!
!     modified 25.11.2007 
!     since Pt1/Pt2=(1+0.5(kappa)-M)**(zeta*kappa)/(kappa-1)
!     and for zeta1 elements type M_crit=M1=1
!     and for zeta2 elements type M_crit=M2 =1
!     it is not necessary to iteratively solve the flow equation.
!     Instead the previous equation is solved to find pt2_crit
!
      if(zeta.ge.0d0) then
         kp1=kappa+1.d0
         km1=kappa-1.d0
         a2a1=a2/a1
         expon1=-0.5d0*kp1/(zeta*kappa)
         expon2=-0.5d0*kp1/km1
         cte=a2a1*(0.5d0*kp1)**expon2
         expon3=-km1/(zeta*kappa)
         i=0
!
!        
         do
            i=i+1
!     
            f=x**(-1.d0)-cte*x**(expon1)
     &           *(2.d0/km1*(x**expon3-1.d0))**(-0.5d0)
!     
            df=-1.d0/X**2-cte*(x**expon1
     &           *(2.d0/km1*(x**expon3-1.d0))**(-0.5d0))
     &           *(expon1/X-1.d0/km1*expon3*x**(expon3-1.d0)
     &           *(2.d0/km1*(x**expon3-1.d0))**(-1.d0))
            
            delta_x=-f/df
!     
            if(( dabs(delta_x/x).le.1.d-8)
     &           .or.(dabs(delta_x/1d0).le.1.d-10)) then
!
               pt2_lim=pt1*X
!
               exit
            endif
            if(i.gt.25)then
                pt2_lim=pt1/(1.d0+0.5d0*km1)**(zeta*kappa/km1)
                exit
             endif
!     
            x=delta_x+x
!     
         enddo
!
      else
!
         do 
            kp1=kappa+1.d0
            km1=kappa-1.d0
            a2a1=a2/a1
            expon1=kp1/(zeta*kappa)
            expon2=km1/(zeta*kappa)
            expon3=kp1/km1
            cte=a2a1**2*(0.5d0*kp1)**(-expon3)*(2.d0/km1)**(-1.d0)
            fact1=x**(-expon1)
            fact2=x**(-expon2)
            term=fact2-1.d0
!     
            f=x**(-2.d0)-cte*fact1*term**(-1.d0)
!     
            df=-2.d0*x**(-3.d0)-cte*(x**(-expon1-1.d0)*term**(-1.d0))
     &           *(-expon1+expon2*(X**(-expon2))*fact2*term**(-1.d0))
!     
            delta_x=-f/df
!     
            if(( dabs(delta_x/x).le.1.d-8)
     &           .or.(dabs(delta_x/1d0).le.1.d-10)) then
               pt2_lim=pt1*X
               exit
            endif
!     
            x=delta_x+x
!     
         enddo
!     
      endif
      
      return
      end
