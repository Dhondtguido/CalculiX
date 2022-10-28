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
!     calculate the maximal admissible pressure ratio pt2/pt1
!
!     1) for icase=1 (sonic conditions are possible at 1):
!                    assuming M1=1 M2 is calculated iteratively
!                    using a dichotomy scheme
!        for icase=2 (sonic conditions are possible at 2):      
!                    assuming M2=1 M1 is calculated iteratively
!                    using a dichotomy scheme
!
!     2)the ratio of the critical pressure ratio  
!       Qred_1/Qred_2crit=Pt2/Pt1=D(M1)/D(M2_crit)
!       is computed 
!       [D(M)=M*(1+0.5*(kappa-1)*M)**(-0.5*(kappa+1)/(kappa-1))]
!       (general gas equation)
!
!     author: Yannick Muller/Guido Dhondt
!
!     Attention: this routine is only used for rotating gas pipes
!                with rotational speed 0, i.e. static pipes with a possibly
!                varying cross section. In such pipes the total pressure
!                always decreases versus length in flow direction
!                beta > 0: flow from node1 to node2
!                beta < 0: flow from node2 to node1
!      
      subroutine pt2zpt1_rot(pt2,pt1,kappa,r,l,pt2zpt1_c,crit,icase,
     &     M1,M2,ca,cb,cc,alpha,beta,Qred1_crit,Qred2_crit,A1,A2)
!     
      implicit none
!
      logical crit
!
      integer icase,i
!     
      real*8 pt2,pt1,kappa,l,M1,pt2zpt1,pt2zpt1_c,km1,kp1,r,Qred1_crit,
     &     fmin,Qred2_crit,f,fmax,M1_min,M1_max,Z1,Z1_min,Z1_max,Z2,
     &     Z2_min,Z2_max,M2,M2_min,M2_max,tcdkm1,km1d2,ca,cb,cc,alpha,
     &     beta,A1,A2,M1_root,M2_root
!
!
!
      crit=.false.
!
!     useful variables and constants
!
      km1=kappa-1.d0
      kp1=kappa+1.d0
      km1d2=km1/2.d0
      tcdkm1=2.d0*cc/(kappa-1.d0)
!
      if(icase.eq.1) then
!
!        M1=1         
!        computing M2 using dichotomy method (dividing the interval
!        with the funciton root iteratively by 2)
!     
         i=1
!
!        because of
!        (alpha+beta*Z2_min)/(alpha+beta))**(cb/beta)
!        we have to avoid a root in the initial dichotomy interval
!
!        for icase.eq.1 we have:
!     alpha+beta*Z is somewhere in the interval [0,1] negative;
!         
!        if beta >= 0, alpha<0 must be true;
!        since alpha+beta*Z is monotonically increasing, a root in the
!        interval is only possible for alpha+beta>0
!         
!        if beta <= 0, alpha+beta<0 must be true;
!        since alpha+beta*Z is monotonically decreasing, a root in the
!        interval is only possible for alpha>0
!
         if(beta.ge.0.d0) then
            if(alpha+beta.gt.0.d0) then
               M2_root=dsqrt(-alpha/beta)
               M2_min=max(0.001d0,M2_root+0.001d0)
               M2_max=0.999d0
            else
               M2_min=0.001d0
               M2_max=0.999d0
            endif
         else
            if(alpha.gt.0.d0) then
               M2_root=dsqrt(-alpha/beta)
               M2_min=max(0.001d0,M2_root+0.001d0)
               M2_max=0.999d0
            else
               M2_min=0.001d0
               M2_max=0.999d0
            endif
         endif
               
c         M2_min=0.001d0
c         M2_max=1
         Z2_min=M2_min**2
         Z2_max=M2_max**2
!
         fmin=dlog(Z2_min**ca*
     &        ((alpha+beta*Z2_min)/(alpha+beta))**(cb/beta)*
     &        ((1.d0+km1d2*Z2_min)/(1.d0+km1d2))**(tcdkm1))-l
!     
         fmax=dlog(Z2_max**ca*
     &        ((alpha+beta*Z2_max)/(alpha+beta))**(cb/beta)*
     &        ((1.d0+km1d2*Z2_max)/(1.d0+km1d2))**(tcdkm1))-l
!
         if(fmin*fmax.gt.0.d0) then
            pt2zpt1_c=1.d30
            Qred2_crit=1.d30
            crit=.false.
            return
         endif
!
         do
            i=i+1
            M2=(M2_min+M2_max)*0.5d0
            Z2=M2**2
!     
            f=dlog(Z2**ca*
     &           ((alpha+beta*Z2)/(alpha+beta))**(cb/beta)*
     &           ((1.d0+km1d2*Z2)/(1.d0+km1d2))**(tcdkm1))-l
!     
            if(abs(f).le.1d-6) then
               exit
            endif
            if(i.gt.50) then
               exit
            endif
!     
            if(fmin*f.le.0.d0) then
               M2_max=M2
               fmax=f
            else
               M2_min=M2
               fmin=f
            endif
         enddo
         write(*,*) 'pt2zpt1_rot M2 ',M2
!     
         pt2zpt1_c=(0.5d0*kp1)**(-0.5d0*kp1/km1)
     &        *(1.d0+0.5d0*km1*Z2)**(0.5d0*kp1/km1)*A1/(A2*M2)
!     
         Qred2_crit=M2*dsqrt(kappa/r)
     &        *(1.d0+0.5d0*km1*Z2)**(-0.5d0*kp1/km1)
!     
         pt2zpt1=pt2/pt1
         write(*,*) 'pt2zpt1_rot pt2/pt1 ',pt2zpt1_c,pt2zpt1
         if(beta.ge.0.d0) then
            if(pt2zpt1.le.pt2zpt1_c) then
               crit=.true.
            endif
         else
            if(pt2zpt1.ge.pt2zpt1_c) then
               crit=.true.
            endif
         endif
!     
      elseif (icase.eq.2) then
!
!        M2=1
!        computing M1 using dichotomy method (dividing the interval
!        with the function root iteratively by 2)
!
         i=1
!         
!        because of
!        (alpha+beta*Z2_min)/(alpha+beta))**(cb/beta)
!        we have to avoid a root in the initial dichotomy interval
!
!        for icase.eq.2 we have:
!     alpha+beta*Z is somewhere in the interval [0,1] positive;
         
!        if beta >= 0 alpha+beta*Z is monotonically increasing, and a root
!        in the interval is only possible for alpha<0
         
!        if beta <= 0 alpha+beta*Z is monotonically decreasing, and a root
!        in the interval is only possible for alpha+beta<0
!
!
         if(beta.ge.0.d0) then
            if(alpha.lt.0.d0) then
               M1_root=dsqrt(-alpha/beta)
               M1_min=max(0.001d0,M1_root+0.001d0)
               M1_max=0.999d0
            else
               M1_min=0.001d0
               M1_max=0.999d0
            endif
         else
            if(alpha+beta.lt.0.d0) then
               M1_root=dsqrt(-alpha/beta)
               M1_min=max(0.001d0,M1_root+0.001d0)
               M1_max=0.999d0
            else
               M1_min=0.001d0
               M1_max=0.999d0
            endif
         endif
c         M1_min=0.001d0
c         M1_max=1
         Z1_min=M1_min**2
         Z1_max=M1_max**2
!
         fmin=dlog((1.d0/Z1_min)**ca*
     &        ((alpha+beta)/(alpha+beta*Z1_min))**(cb/beta)*
     &        ((1.d0+km1d2)/(1.d0+km1d2*Z1_min))**(tcdkm1))-l
!     
         fmax=dlog((1.d0/Z1_max)**ca*
     &        ((alpha+beta)/(alpha+beta*Z1_max))**(cb/beta)*
     &        ((1.d0+km1d2)/(1.d0+km1d2*Z1_max))**(tcdkm1))-l
!
         if(fmin*fmax.gt.0.d0) then
            pt2zpt1_c=1.d30
            Qred1_crit=1.d30
            crit=.false.
            return
         endif
!
         do
            i=i+1
            M1=(M1_min+M1_max)*0.5d0
            Z1=M1**2
!     
            f=dlog((1.d0/Z1)**ca*
     &           ((alpha+beta)/(alpha+beta*Z1))**(cb/beta)*
     &           ((1.d0+km1d2)/(1.d0+km1d2*Z1))**(tcdkm1))-l
!     
            if(abs(f).le.1d-6) then
               exit
            endif
            if(i.gt.50) then
               exit
            endif
!     
            if(fmin*f.le.0.d0) then
               M1_max=M1
               fmax=f
            else
               M1_min=M1
               fmin=f
            endif
         enddo
!     
         pt2zpt1_c=M1*(0.5d0*kp1)**(0.5d0*kp1/km1)
     &        *(1.d0+0.5d0*km1*Z1)**(-0.5d0*kp1/km1)*A1/A2
!     
         Qred1_crit=M1*dsqrt(kappa/r)
     &        *(1.d0+0.5d0*km1*Z1)**(-0.5d0*kp1/km1)
!     
         pt2zpt1=pt2/pt1
         if(beta.ge.0.d0) then
            if(pt2zpt1.le.pt2zpt1_c) then
               crit=.true.
            endif
         else
            if(pt2zpt1.ge.pt2zpt1_c) then
               crit=.true.
            endif
         endif
!     
      endif
!     
      return
      end      
      
      
