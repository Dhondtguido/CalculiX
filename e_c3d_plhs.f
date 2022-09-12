!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine e_c3d_plhs(lakonl,sm,nelem,ipvar,var)
!     
!     computation of the pressure element matrix for the element with
!     the topology in konl
!     
      implicit none
!     
      character*8 lakonl
!     
      integer nelem,i,j,k,nope,mint3d,ipvar(*),index
!     
      real*8 shp(4,20),sm(8,8),weight,var(*)
!     
      if(lakonl(4:4).eq.'8') then
        nope=8
      elseif(lakonl(4:4).eq.'4') then
        nope=4
      elseif(lakonl(4:4).eq.'6') then
        nope=6
      endif
!     
      if(lakonl(4:5).eq.'8R') then
        mint3d=1
      elseif(lakonl(4:4).eq.'8') then
        mint3d=8
      elseif(lakonl(4:4).eq.'4') then
        mint3d=1
      elseif(lakonl(4:5).eq.'6 ') then
        mint3d=2
      endif
!     
!     initialisation of sm
!     
      do i=1,nope
        do j=1,nope
          sm(i,j)=0.d0
        enddo
      enddo
!     
!     computation of the matrix: loop over the Gauss points
!     
      index=ipvar(nelem)
      do k=1,mint3d
!     
!     copying the shape functions, their derivatives and the
!     Jacobian determinant from field var
!     
        do j=1,nope
          do i=1,4
            index=index+1
            shp(i,j)=var(index)
          enddo
        enddo
!     
!     Jacobian determinant
!     
        index=index+1
        weight=var(index)
        index=index+1
        do j=1,nope
!     
          do i=1,j
!     
!     lhs pressure matrix (only for semi-implicit incompressible
!     calculations)
!     
            sm(i,j)=sm(i,j)
     &           +(shp(1,i)*shp(1,j)+
     &           shp(2,i)*shp(2,j)+
     &           shp(3,i)*shp(3,j))*weight
          enddo
        enddo
      enddo
!     
      return
      end
