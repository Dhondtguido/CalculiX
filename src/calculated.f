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
      subroutine calculated(nktet,d,dmin,ipoed,iedg,cotet)
!
!     determine the length of all edges in the actual mesh 
!
      implicit none
!
      integer i,nktet,ipoed(*),iedg(3,*),index,n1,n2
!
      real*8 d(*),dmin,cotet(3,*)
!
!     determine the size of all edges
!
      dmin=1.d30
!
      loop: do i=1,nktet
         index=ipoed(i)
         do
            if(index.eq.0) cycle loop
!
            n1=iedg(1,index)
            n2=iedg(2,index)
!
            d(index)=dsqrt((cotet(1,n1)-cotet(1,n2))**2+
     &                     (cotet(2,n1)-cotet(2,n2))**2+
     &                     (cotet(3,n1)-cotet(3,n2))**2)
!
            if(d(index).lt.dmin) dmin=d(index)
!
            index=iedg(3,index)
         enddo
      enddo loop
!
      return
      end

