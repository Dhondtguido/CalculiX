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
      subroutine createfint(ne,ipkon,lakon,kon,nactdof,mi,
     &   fn0,fint)
!
      implicit none
!
      character*8 lakon(*)
! 
      integer i,j,k,ne,indexe,ipkon(*),nope,node,kon(*),mi(*),
     &  nactdof(0:mi(2),*)
!      
      real*8 fn0(0:mi(2),*),fint(*)
!
!
!
!     calculation of the internal force vector for all
!     active degrees of freedom
!     
      do i=1,ne
! 
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
c        Bernhardi start
         if(lakon(i)(1:5).eq.'C3D8I') then
            nope=11
         elseif(lakon(i)(4:5).eq.'20') then
c        Bernhardi end
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         endif
         
         do j=1,nope
            node=kon(indexe+j)
            do k=1,3
               if(nactdof(k,node).gt.0) then
                  fint(nactdof(k,node))=fint(nactdof(k,node))
     &                 +fn0(k,indexe+j)
               endif
            enddo
         enddo
      enddo
!      
      return
      end
