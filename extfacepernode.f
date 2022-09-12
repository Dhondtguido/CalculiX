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
      subroutine extfacepernode(iponoelfa,inoelfa,lakonfa,ipkonfa,
     &       konfa,nsurfs,inoelsize)
!
      implicit none
!
      character*8 lakonfa(*)
!
      integer iponoelfa(*),inoelfa(3,*),ipkonfa(*),konfa(*),i,j,nsurfs,
     &  inoelfree,nope,indexe,node,inoelsize
!
!
!
!     lists which external faces correspond to a given node i
!     iponoelfa(i) points to an entry j in field inoelfa where:
!     inoelfa(1,j): face number as catalogued in fields konfa, lakonfa
!     inoelfa(2,j): local node number in the topology description
!     inoelfa(3,j): pointer to the next face to which i belongs, or, if
!                   none is left: zero
!     
      inoelfree=1
      do i=1,nsurfs
         if(ipkonfa(i).lt.0) cycle
         if(lakonfa(i)(2:2).eq.'8') then
            nope=8
         elseif(lakonfa(i)(2:2).eq.'4') then
            nope=4
         elseif(lakonfa(i)(2:2).eq.'6') then
            nope=6
         elseif(lakonfa(i)(2:2).eq.'3') then
            nope=3
         endif
         indexe=ipkonfa(i)
         do j=1,nope
            node=konfa(indexe+j)
            inoelfa(1,inoelfree)=i
            inoelfa(2,inoelfree)=j
            inoelfa(3,inoelfree)=iponoelfa(node)
            iponoelfa(node)=inoelfree
            inoelfree=inoelfree+1
         enddo
      enddo
!
!     size of field inoelfa
!
      inoelsize=inoelfree-1
!
      return
      end
