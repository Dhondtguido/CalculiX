!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine e_c3d_vlhs(lakonl,sm,nelem,ipvar,var)
!
!     computation of the velocity element matrix for the element with
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
      if(lakonl(4:4).ne.'4') then
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
c            do i=1,4
c              index=index+1
              index=index+4
c              shp(i,j)=var(index)
              shp(4,j)=var(index)
c            enddo
          enddo
          index=index+1
          weight=var(index)
!     
          index=index+1
!     
          do j=1,nope
            do i=1,j
!     
!     lhs velocity matrix
!     
              sm(i,j)=sm(i,j)
     &             +shp(4,i)*shp(4,j)*weight
            enddo
          enddo
        enddo
      else
!
!        C3D4: analytical solution (agrees with a 4 point integration
!              scheme
!
        index=ipvar(nelem)+17
        weight=var(index)
        do j=1,4
          do i=1,j-1
            sm(i,j)=0.05d0*weight
          enddo
          sm(j,j)=0.1d0*weight
        enddo
      endif
!
      return
      end

