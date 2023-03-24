!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine cattri(ne,lakon,ipkon,kon,kontri,ntri,mastelnr)
!     
!     catalogueing the tetrahedral elements of the mesh
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer i,j,ne,ipkon(*),kon(*),indexe,ntri,kontri(3,*),mastelnr(*)
!     
!     catalogue the crack elements
!     
      ntri=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(lakon(i).ne.'C3D6  L ') cycle
        ntri=ntri+1
        indexe=ipkon(i)
        do j=1,3
          kontri(j,ntri)=kon(indexe+j+6)
        enddo
        mastelnr(ntri)=i
c        write(*,*) 'cattri ',ntri,(kontri(j,ntri),j=1,3)
      enddo
!     
      return
      end

