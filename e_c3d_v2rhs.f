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
      subroutine e_c3d_v2rhs(konl,lakonl,bb,nelem,v,dtimef,mi,
     &     ipvar,var,nk)
!     
!     computation of the velocity element matrix and rhs for the element with
!     element with the topology in konl: step 3 (correction **)
!     
!     bb: rhs 
!     
      implicit none
!     
      character*8 lakonl
!     
      integer konl(8),nelem,i,j,k,nope,mint3d,mi(*),ipvar(*),index,nk
!     
      real*8 shp(4,8),bb(3,8),xsjmod,vl(0:mi(2),8),v(nk,0:mi(2)),
     &     var(*),dtimef,xsj,ddpress(3)
!     
      if(lakonl(4:4).eq.'4') then
        nope=4
        mint3d=1
      elseif(lakonl(4:4).eq.'6') then
        nope=6
        mint3d=2
      elseif(lakonl(4:5).eq.'8R') then
        nope=8
        mint3d=1
      elseif(lakonl(4:4).eq.'8') then
        nope=8
        mint3d=8
      endif
!     
!     initialisation of the rhs
!     
      do i=1,nope
        do j=1,3
          bb(j,i)=0.d0
        enddo
      enddo
!     
!     change in pressure
!     
      do i=1,nope
        vl(4,i)=v(konl(i),4)
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
        index=index+1
        xsj=var(index)
!     
        xsjmod=dtimef*xsj
!     
        index=index+1
!     
!     only for the semi-implicit procedure: calculate ddpress;
!     
        do j=1,3
          ddpress(j)=0.d0
        enddo
        do i=1,nope
          do j=1,3
            ddpress(j)=ddpress(j)+shp(j,i)*vl(4,i)
          enddo
        enddo
!     
!     3. contribution to the second part of the
!     momentum equation
!     
        do j=1,nope
          bb(1,j)=bb(1,j)-xsjmod*shp(4,j)*ddpress(1)
          bb(2,j)=bb(2,j)-xsjmod*shp(4,j)*ddpress(2)
          bb(3,j)=bb(3,j)-xsjmod*shp(4,j)*ddpress(3)
        enddo
!     
      enddo
!
c      write(*,*) 'e_c3d_v2rhs'
c      do k=1,8
c        write(*,*) nelem,k,(bb(j,k),j=1,3)
c      enddo
!     
      return
      end
