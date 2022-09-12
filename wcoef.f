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
      subroutine wcoef(v,vo,al,um)
!
!     computation of the coefficients of w in the derivation of the
!     second order element stiffness matrix
!
      implicit none
!
      real*8 v(3,3,3,3),vo(3,3)
!
      real*8 a2u,al,um,au,p1,p2,p3
!
!
!
      a2u=al+2.d0*um
      au=al+um
!
      p1=vo(1,1)+1.d0
      p2=vo(2,2)+1.d0
      p3=vo(3,3)+1.d0
!
      v(1,1,1,1)=a2u*p1*p1+um*(vo(1,2)**2+vo(1,3)**2)
      v(2,1,1,1)=au*vo(1,2)*p1
      v(3,1,1,1)=au*vo(1,3)*p1
      v(1,2,1,1)=v(2,1,1,1)
      v(2,2,1,1)=a2u*vo(1,2)**2+um*(p1*p1+vo(1,3)**2)
      v(3,2,1,1)=au*vo(1,2)*vo(1,3)
      v(1,3,1,1)=v(3,1,1,1)
      v(2,3,1,1)=v(3,2,1,1)
      v(3,3,1,1)=a2u*vo(1,3)**2+um*(p1*p1+vo(1,2)**2)
!
      v(1,1,2,1)=al*vo(2,1)*p1+
     &  um*(2.d0*vo(2,1)*p1+vo(1,2)*p2+vo(2,3)*vo(1,3))
      v(2,1,2,1)=al*p1*p2+um*vo(2,1)*vo(1,2)
      v(3,1,2,1)=al*vo(2,3)*p1+um*vo(2,1)*vo(1,3)
      v(1,2,2,1)=al*vo(2,1)*vo(1,2)+um*p1*p2
      v(2,2,2,1)=al*vo(1,2)*p2+
     &  um*(vo(2,1)*p1+2.d0*vo(1,2)*p2+vo(2,3)*vo(1,3))
      v(3,2,2,1)=al*vo(2,3)*vo(1,2)+um*vo(1,3)*p2
      v(1,3,2,1)=al*vo(2,1)*vo(1,3)+um*vo(2,3)*p1
      v(2,3,2,1)=al*vo(1,3)*p2+um*vo(2,3)*vo(1,2)
      v(3,3,2,1)=a2u*vo(2,3)*vo(1,3)+
     &  um*(vo(2,1)*p1+vo(1,2)*p2)
!
      v(1,1,3,1)=al*vo(3,1)*p1+
     &  um*(vo(1,3)*p3+2.d0*vo(3,1)*p1+vo(3,2)*vo(1,2))
      v(2,1,3,1)=al*vo(3,2)*p1+um*vo(3,1)*vo(1,2)
      v(3,1,3,1)=al*p1*p3+um*vo(3,1)*vo(1,3)
      v(1,2,3,1)=al*vo(3,1)*vo(1,2)+um*vo(3,2)*p1
      v(2,2,3,1)=a2u*vo(3,2)*vo(1,2)+
     &  um*(vo(1,3)*p3+vo(3,1)*p1)
      v(3,2,3,1)=al*vo(1,2)*p3+um*vo(3,2)*vo(1,3)
      v(1,3,3,1)=al*vo(3,1)*vo(1,3)+um*p1*p3
      v(2,3,3,1)=al*vo(3,2)*vo(1,3)+um*vo(1,2)*p3
      v(3,3,3,1)=al*vo(1,3)*p3+
     &  um*(2.d0*vo(1,3)*p3+vo(3,1)*p1+vo(3,2)*vo(1,2))
!
      v(1,1,1,2)=al*vo(2,1)*p1+
     &  um*(vo(1,2)*p2+2.d0*vo(2,1)*p1+vo(1,3)*vo(2,3))
      v(2,1,1,2)=al*vo(1,2)*vo(2,1)+um*p1*p2
      v(3,1,1,2)=al*vo(1,3)*vo(2,1)+um*vo(2,3)*p1
      v(1,2,1,2)=al*p1*p2+um*vo(1,2)*vo(2,1)
      v(2,2,1,2)=al*vo(1,2)*p2+
     &  um*(2.d0*vo(1,2)*p2+vo(2,1)*p1+vo(1,3)*vo(2,3))
      v(3,2,1,2)=al*vo(1,3)*p2+um*vo(1,2)*vo(2,3)
      v(1,3,1,2)=al*vo(2,3)*p1+um*vo(1,3)*vo(2,1)
      v(2,3,1,2)=al*vo(1,2)*vo(2,3)+um*vo(1,3)*p2
      v(3,3,1,2)=a2u*vo(1,3)*vo(2,3)+
     &  um*(vo(1,2)*p2+vo(2,1)*p1)
!
      v(1,1,2,2)=a2u*vo(2,1)**2+um*(p2*p2+vo(2,3)**2)
      v(2,1,2,2)=au*vo(2,1)*p2
      v(3,1,2,2)=au*vo(2,3)*vo(2,1)
      v(1,2,2,2)=v(2,1,2,2)
      v(2,2,2,2)=a2u*p2*p2+um*(vo(2,1)**2+vo(2,3)**2)
      v(3,2,2,2)=au*vo(2,3)*p2
      v(1,3,2,2)=v(3,1,2,2)
      v(2,3,2,2)=v(3,2,2,2)
      v(3,3,2,2)=a2u*vo(2,3)**2+um*(p2*p2+vo(2,1)**2)
!
      v(1,1,3,2)=a2u*vo(3,1)*vo(2,1)+
     &  um*(vo(3,2)*p2+vo(2,3)*p3)
      v(2,1,3,2)=al*vo(3,2)*vo(2,1)+um*vo(3,1)*p2
      v(3,1,3,2)=al*vo(2,1)*p3+um*vo(3,1)*vo(2,3)
      v(1,2,3,2)=al*vo(3,1)*p2+um*vo(3,2)*vo(2,1)
      v(2,2,3,2)=al*vo(3,2)*p2+
     &  um*(2.d0*vo(3,2)*p2+vo(2,3)*p3+vo(3,1)*vo(2,1))
      v(3,2,3,2)=al*p2*p3+um*vo(3,2)*vo(2,3)
      v(1,3,3,2)=al*vo(3,1)*vo(2,3)+um*vo(2,1)*p3
      v(2,3,3,2)=al*vo(3,2)*vo(2,3)+um*p2*p3
      v(3,3,3,2)=al*vo(2,3)*p3+
     &  um*(vo(3,2)*p2+2.d0*vo(2,3)*p3+vo(3,1)*vo(2,1))
!
      v(1,1,1,3)=al*vo(3,1)*p1+
     &  um*(vo(1,3)*p3+2.d0*vo(3,1)*p1+vo(1,2)*vo(3,2))
      v(2,1,1,3)=al*vo(1,2)*vo(3,1)+um*vo(3,2)*p1
      v(3,1,1,3)=al*vo(1,3)*vo(3,1)+um*p1*p3
      v(1,2,1,3)=al*vo(3,2)*p1+um*vo(1,2)*vo(3,1)
      v(2,2,1,3)=a2u*vo(1,2)*vo(3,2)+
     &  um*(vo(1,3)*p3+vo(3,1)*p1)
      v(3,2,1,3)=al*vo(1,3)*vo(3,2)+um*vo(1,2)*p3
      v(1,3,1,3)=al*p1*p3+um*vo(1,3)*vo(3,1)
      v(2,3,1,3)=al*vo(1,2)*p3+um*vo(1,3)*vo(3,2)
      v(3,3,1,3)=al*vo(1,3)*p3+
     &  um*(2.d0*vo(1,3)*p3+vo(3,1)*p1+vo(1,2)*vo(3,2))
!
      v(1,1,2,3)=a2u*vo(2,1)*vo(3,1)+
     &  um*(vo(2,3)*p3+vo(3,2)*p2)
      v(2,1,2,3)=al*vo(3,1)*p2+um*vo(2,1)*vo(3,2)
      v(3,1,2,3)=al*vo(2,3)*vo(3,1)+um*vo(2,1)*p3
      v(1,2,2,3)=al*vo(2,1)*vo(3,2)+um*vo(3,1)*p2
      v(2,2,2,3)=al*vo(3,2)*p2+
     &  um*(vo(2,3)*p3+2.d0*vo(3,2)*p2+vo(2,1)*vo(3,1))
      v(3,2,2,3)=al*vo(2,3)*vo(3,2)+um*p2*p3
      v(1,3,2,3)=al*vo(2,1)*p3+um*vo(2,3)*vo(3,1)
      v(2,3,2,3)=al*p2*p3+um*vo(2,3)*vo(3,2)
      v(3,3,2,3)=al*vo(2,3)*p3+
     &  um*(2.d0*vo(2,3)*p3+vo(3,2)*p2+vo(2,1)*vo(3,1))
!
      v(1,1,3,3)=a2u*vo(3,1)**2+um*(p3*p3+vo(3,2)**2)
      v(2,1,3,3)=au*vo(3,2)*vo(3,1)
      v(3,1,3,3)=au*vo(3,1)*p3
      v(1,2,3,3)=v(2,1,3,3)
      v(2,2,3,3)=a2u*vo(3,2)**2+um*(p3*p3+vo(3,1)**2)
      v(3,2,3,3)=au*vo(3,2)*p3
      v(1,3,3,3)=v(3,1,3,3)
      v(2,3,3,3)=v(3,2,3,3)
      v(3,3,3,3)=a2u*p3*p3+um*(vo(3,1)**2+vo(3,2)**2)
!
      return
      end
