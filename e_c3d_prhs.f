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
      subroutine e_c3d_prhs(nk,konl,lakonl,ff,nelem,v,dtimef,mi,theta1,
     &     ipvar,var)
!     
!     computation of the pressure element matrix and rhs for the element with
!     element with the topology in konl: step 2
!     
!     sm: lhs matrix
!     ff: rhs 
!     
      implicit none
!     
      character*8 lakonl
!     
      integer konl(8),nk,nelem,index,i,j,k,ipvar(*),
     &     mi(*),nope,mint3d
!     
      real*8 shp(4,8),ff(8),theta1,xsjmod,vl(0:mi(2),8),v(nk,0:mi(2)),
     &     xsj,delrhovel(3),aux1(3),var(*),dtimef
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
        ff(i)=0.d0
      enddo
!     
!     first part of change of momentum (V*)
!
      do i=1,nope
        do j=1,3
          vl(j,i)=v(konl(i),j)
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
        index=index+1
        xsj=var(index)
        index=index+1
!     
        xsjmod=dtimef*xsj
!     
!     change of rho*V^*
!     
        do j=1,3
          delrhovel(j)=0.d0
        enddo
!        
        do i=1,nope
          do j=1,3
            delrhovel(j)=delrhovel(j)+shp(4,i)*vl(j,i)
          enddo
        enddo
!     
        do j=1,3
          aux1(j)=xsjmod*theta1*delrhovel(j)
        enddo
!     
!     determination of rhs
!     
        do j=1,nope
!     
          ff(j)=ff(j)+
     &         shp(1,j)*aux1(1)+shp(2,j)*aux1(2)+shp(3,j)*aux1(3)
        enddo
!     
      enddo
c      write(*,*) 'e_c3d_prhs'
c      do k=1,8
c        write(*,*) nelem,k,ff(k)
c      enddo
!     
      return
      end
