!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine calcvol(n1,n2,n3,n4,cotet,volume)
!
!     calculates the volume of a tetrahedral element with nodes
!     n1, n2, n3 and n4
!
      implicit none
!
      integer n1,n2,n3,n4
!
      real*8 cotet(3,*),s(3),volume
!
      s(1)=(cotet(2,n2)-cotet(2,n1))*(cotet(3,n3)-cotet(3,n1))-
     &     (cotet(3,n2)-cotet(3,n1))*(cotet(2,n3)-cotet(2,n1))
      s(2)=(cotet(3,n2)-cotet(3,n1))*(cotet(1,n3)-cotet(1,n1))-
     &     (cotet(1,n2)-cotet(1,n1))*(cotet(3,n3)-cotet(3,n1))
      s(3)=(cotet(1,n2)-cotet(1,n1))*(cotet(2,n3)-cotet(2,n1))-
     &     (cotet(2,n2)-cotet(2,n1))*(cotet(1,n3)-cotet(1,n1))
      volume=((cotet(1,n4)-cotet(1,n1))*s(1)+
     &        (cotet(2,n4)-cotet(2,n1))*s(2)+
     &        (cotet(3,n4)-cotet(3,n1))*s(3))/6.d0
!
      return
      end
