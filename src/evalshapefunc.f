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
      subroutine evalshapefunc(xil,etl,xl2,nopes,p)
!     
      implicit none
!     
      integer i,j,nopes,iflag
!     
      real*8  xl2(3,8),xs2(3,2),shp2(7,8),p(3),xsj2(3),
     &  xil,etl
!     
      iflag=1
      do j=1,3
         p(j)=0.d0
      enddo 
      if(nopes.eq.8)then
         call shape8q(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.4)then
         call shape4q(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.6)then
         call shape6tri(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      else
         call shape3tri(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      endif 
      do i=1,nopes
         do j=1,3
            p(j)=p(j)+xl2(j,i)*shp2(4,i)
         enddo
      enddo 
!     
      return
      end
