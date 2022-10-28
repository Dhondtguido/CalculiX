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
      subroutine cubic(a0,a1,a2,solreal,solimag,n)
!
!     solves the cubic equation x**3+a2*x**2+a1*x+a0=0
!     (analytical solution)
!     there are three solutions in C (complex plane)
!     the real part of the solutions is stored in solreal(1..3),
!     the imaginary part in solimag(1..3). The real solutions have
!     the lower indices, their number is n
!
!     Reference: Abramowitz, M. and Stegun, I., Handbook of 
!     Mathematical Functions (10th printing, 1972 p 17)
!     
      implicit none
!      
      integer n
!
      real*8 a0,a1,a2,solreal(3),solimag(3),q,r,d,s1,s2,
     &  a,phi,s1r,s1i
!
      write(30,*) 'a2,a1,a0 ',a2,a1,a0
      q=a1/3.d0-a2*a2/9.d0
      r=(a1*a2-3.d0*a0)/6.d0-(a2**3)/27.d0
!
      d=q**3+r*r
      write(30,*) 'q,r,d ',q,r,d
!
      if(d.gt.0) then
!
!        one real solution, two complex conjugate complex
!        solutions
!
         n=1
         s1=(r+dsqrt(d))**(1.d0/3.d0)
         s2=(r-dsqrt(d))
         if(s2.gt.0.d0) then
            s2=s2**(1.d0/3.d0)
         else
            s2=-(-s2)**(1.d0/3.d0)
         endif
        write(30,*) 'd>0 s1,s2 ',s1,s2
!
         solreal(1)=(s1+s2)-a2/3.d0
         solreal(2)=-(s1+s2)/2.d0-a2/3.d0
         solreal(3)=solreal(2)
!
         solimag(1)=0.d0
         solimag(2)=(s1-s2)*dsqrt(3.d0)/2.d0
         solimag(3)=-solimag(2)
      else
!
!        three real solutions
!
         n=3
!
!        amplitude and phase of s1
!         
         a=(r*r-d)**(1.d0/6.d0)
         phi=(datan2(dsqrt(-d),r))/3.d0
c         phi=(datan(dsqrt(-d)/r))/3.d0
        write(30,*) 'd <=0 a,phi ',a,phi
!
!        real and imaginary part of s1
!
         s1r=a*dcos(phi)
         s1i=a*dsin(phi)
        write(30,*) 'd >=0 s1r,s1i ',s1r,s1i
!
         solreal(1)=2.d0*s1r-a2/3.d0
         solreal(2)=-s1r-a2/3.d0-s1i*dsqrt(3.d0)
         solreal(3)=-s1r-a2/3.d0+s1i*dsqrt(3.d0)
!
         solimag(1)=0.d0
         solimag(2)=0.d0
         solimag(3)=0.d0
      endif
!
      return
      end
      

