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
!     d{K(X)}/dTeint
!
!     author: Yannick Muller
!     
      subroutine dKdt(x,u,uprime,rpar,ipar)
!     
      implicit none
      integer ipar
      real*8 x,u(1),uprime(1),rpar(*),zk0,phi,Tup,
     &     xflow,Pup,f1_x,K_x,lambda1,df1dk,Rurd,f_k,kup
!
      external f_k
!
!     defining the parameters
      phi=rpar(1)
      lambda1=rpar(2)
      zk0=rpar(3)
      Pup=rpar(4)
      Tup=rpar(5)
      rurd=rpar(6)
      xflow=rpar(7)
      kup=rpar(8)
!     
!     find K(X) for the given x

      k_x=f_k(x,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
!
      k_x=dsqrt(K_x/x)
!
!     f1_x
      f1_x= (zk0*K_x)**(7.d0/4.d0)
     &     -(1-K_x)/dabs(1-K_x)*dabs(1-K_x)**(7d0*4d0)
!     
!     df1dK
      df1dK=7d0/4d0*zk0**(7d0/4d0)*K_x**(3.d0/4.d0)
     &     +7d0/4d0*dabs(1-K_x)**(3.d0/4.d0)
!     
!     
      uprime(1)=-0.8d0*x**1.6d0*lambda1*Pup**(0.8d0)
     &     /(xflow*Tup**1.8d0)*f1_x+u(1)
     &     *(lambda1*x**1.6d0*Pup**0.8d0/(xflow*Tup**0.8d0)
     &     *df1dK-2/x)
!     
      return
!     
      end
