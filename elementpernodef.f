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
      subroutine elementpernodef(iponoel,inoel,lakonf,ipkonf,konf,nef)
!
      implicit none
!
      character*8 lakonf(*)
!
      integer iponoel(*),inoel(2,*),ipkonf(*),konf(*),i,j,nef,
     &  inoelfree,nope,indexe,node
!
!
!
!     determining the elements belonging to the nodes of
!     the elements
!
      inoelfree=1
      do i=1,nef
         if(ipkonf(i).lt.0) cycle
         if(lakonf(i)(1:1).ne.'F') cycle
         if(lakonf(i)(4:4).eq.'8') then
            nope=8
         elseif(lakonf(i)(4:4).eq.'4') then
            nope=4
         elseif(lakonf(i)(4:4).eq.'6') then
            nope=6
         else
            cycle
         endif
         indexe=ipkonf(i)
         do j=1,nope
            node=konf(indexe+j)
            inoel(1,inoelfree)=i
            inoel(2,inoelfree)=iponoel(node)
            iponoel(node)=inoelfree
            inoelfree=inoelfree+1
         enddo
      enddo
!
      return
      end
