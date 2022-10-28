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
      subroutine smalldist(co,distmin,lakon,ipkon,kon,ne)
!
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),i,j,ne,n,j1,j2,indexe,neigh2(2,1),
     &  neigh4(2,6),neigh6(2,9),neigh8(2,12),neigh10(2,12),
     &  neigh15(2,18),neigh20(2,24)
!
      real*8 dist,distmin,co(3,*)
!
!
!
      neigh2=reshape((/1,2/),(/2,1/))
      neigh4=reshape((/1,2,2,3,3,1,1,4,2,4,3,4/),(/2,6/))
      neigh6=reshape((/1,2,2,3,3,1,4,5,5,6,6,4,1,4,2,5,3,6/),(/2,9/))
      neigh8=reshape((/1,2,2,3,3,4,4,1,5,6,6,7,7,8,8,5,
     &                 1,5,2,6,3,7,4,8/),(/2,12/))
c      neigh10=reshape((/1,5,5,2,2,6,6,3,3,7,7,1,
c     &                  1,8,8,4,2,9,9,4,3,10,10,4/),(/2,12/))
c      neigh15=reshape((/1,7,7,2,2,8,8,3,3,9,9,1,
c     &                  4,10,10,5,5,11,11,6,6,12,12,4,
c     &                  1,13,13,4,2,14,14,5,3,15,15,6/),(/2,18/))
c      neigh20=
c     &     reshape((/1,9,9,2,2,10,10,3,3,11,11,4,4,12,12,1,
c     &               5,13,13,6,6,14,14,7,7,15,15,8,8,16,16,5,
c     &               1,17,17,5,2,18,18,6,3,19,19,7,4,20,20,8/),(/2,24/))

!     determining the smallest distance between nodes 
!     
      distmin=1.d30
!
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).eq.'F') cycle
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'2') then
            if((lakon(i)(7:7).eq.'A').or.
     &           (lakon(i)(7:7).eq.'L').or.
     &           (lakon(i)(7:7).eq.'S').or.
     &           (lakon(i)(7:7).eq.'E')) then
               n=4
            else
               n=12
            endif
            do j=1,n
               j1=kon(indexe+neigh8(1,j))
               j2=kon(indexe+neigh8(2,j))
               dist=(co(1,j1)-co(1,j2))**2+
     &              (co(2,j1)-co(2,j2))**2+
     &              (co(3,j1)-co(3,j2))**2
               distmin=min(dist/4.d0,distmin)
            enddo
         elseif(lakon(i)(1:8).eq.'ESPRNGA1') then
            n=1
            do j=1,n
               j1=kon(indexe+neigh2(1,j))
               j2=kon(indexe+neigh2(2,j))
               dist=(co(1,j1)-co(1,j2))**2+
     &              (co(2,j1)-co(2,j2))**2+
     &              (co(3,j1)-co(3,j2))**2
               distmin=min(dist,distmin)
            enddo
         elseif(lakon(i)(4:4).eq.'8') then
            if((lakon(i)(7:7).eq.'A').or.
     &           (lakon(i)(7:7).eq.'L').or.
     &           (lakon(i)(7:7).eq.'S').or.
     &           (lakon(i)(7:7).eq.'E')) then
               n=4
            else
               n=12
            endif
            do j=1,n
               j1=kon(indexe+neigh8(1,j))
               j2=kon(indexe+neigh8(2,j))
               dist=(co(1,j1)-co(1,j2))**2+
     &              (co(2,j1)-co(2,j2))**2+
     &              (co(3,j1)-co(3,j2))**2
               distmin=min(dist,distmin)
            enddo
         elseif(lakon(i)(4:4).eq.'4') then
            n=6
            do j=1,n
               j1=kon(indexe+neigh4(1,j))
               j2=kon(indexe+neigh4(2,j))
               dist=(co(1,j1)-co(1,j2))**2+
     &              (co(2,j1)-co(2,j2))**2+
     &              (co(3,j1)-co(3,j2))**2
               distmin=min(dist,distmin)
            enddo
         elseif(lakon(i)(4:5).eq.'10') then
            n=6
            do j=1,n
               j1=kon(indexe+neigh4(1,j))
               j2=kon(indexe+neigh4(2,j))
               dist=(co(1,j1)-co(1,j2))**2+
     &              (co(2,j1)-co(2,j2))**2+
     &              (co(3,j1)-co(3,j2))**2
               distmin=min(dist/4.d0,distmin)
            enddo
         elseif(lakon(i)(4:4).eq.'6') then
            if((lakon(i)(7:7).eq.'A').or.
     &           (lakon(i)(7:7).eq.'L').or.
     &           (lakon(i)(7:7).eq.'S').or.
     &           (lakon(i)(7:7).eq.'E')) then
               n=3
            else
               n=9
            endif
            do j=1,n
               j1=kon(indexe+neigh6(1,j))
               j2=kon(indexe+neigh6(2,j))
               dist=(co(1,j1)-co(1,j2))**2+
     &              (co(2,j1)-co(2,j2))**2+
     &              (co(3,j1)-co(3,j2))**2
               distmin=min(dist,distmin)
            enddo
         elseif(lakon(i)(4:5).eq.'15') then
            if((lakon(i)(7:7).eq.'A').or.
     &           (lakon(i)(7:7).eq.'L').or.
     &           (lakon(i)(7:7).eq.'S').or.
     &           (lakon(i)(7:7).eq.'E')) then
               n=3
            else
               n=9
            endif
            do j=1,n
               j1=kon(indexe+neigh6(1,j))
               j2=kon(indexe+neigh6(2,j))
               dist=(co(1,j1)-co(1,j2))**2+
     &              (co(2,j1)-co(2,j2))**2+
     &              (co(3,j1)-co(3,j2))**2
               distmin=min(dist/4.d0,distmin)
            enddo
         else
            cycle
         endif
      enddo
!
      distmin=dsqrt(distmin)*1.d-6
      distmin=max(distmin,1.d-8)
!
      write(*,'(1x,a14,e14.7)') 'Perturbation: ',distmin
      write(*,*)
!     
      return
      end
