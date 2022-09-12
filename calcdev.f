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
      subroutine calcdev(vold,vcon,v,nk,iturbulent,mi,vconmax,vmax,
     &     iexplicit,nka,nkb)
!     
!     calculates the change in solution
!     
      implicit none
!     
      integer iturbulent,mi(*),nk,i,j,iexplicit,nka,nkb
!     
      real*8 v(nk,0:mi(2)),vold(0:mi(2),*),vcon(nk,0:mi(2)),vmax(0:6),
     &     vconmax(0:6)
!     
!     first subiteration: calculate the size of the conservative
!     fields
!     
      do j=0,6
        vconmax(j)=0.d0
      enddo
!     
      if(iexplicit.eq.1) then
!     
!     for incompressible fluids the density is stored
!     in vcon(4,*), the change in density in v(*,4)
!     
        do i=nka,nkb
          do j=0,4
            vconmax(j)=vconmax(j)+vcon(i,j)**2
          enddo
        enddo
      else
        do i=nka,nkb
          do j=0,3
            vconmax(j)=vconmax(j)+vcon(i,j)**2
          enddo
!     
!     for incompressible fluids the pressure is stored
!     in vold(4,*), the change in pressure in v(*,4)
!     
          vconmax(4)=vconmax(4)+vold(4,i)**2
        enddo
      endif
      if(iturbulent.ne.0) then
        do i=nka,nkb
          do j=5,6
            vconmax(j)=vconmax(j)+vcon(i,j)**2
          enddo
        enddo
      endif
!     
!     all subiterations: calculate the size of the change of
!     the conservative variables
!     
      do j=0,6
        vmax(j)=0.d0
      enddo
!     
      if(iexplicit.eq.1) then
!     
!     for incompressible fluids the density is stored
!     in vcon(*,4), the change in density in v(*,4)
!     
        do i=nka,nkb
          do j=0,4
            vmax(j)=vmax(j)+v(i,j)**2
          enddo
        enddo
      else
        do i=nka,nkb
          do j=0,3
            vmax(j)=vmax(j)+v(i,j)**2
          enddo
!     
!     for incompressible fluids the pressure is stored
!     in vold(4,*), the change in pressure in v(*,4)
!     
          vmax(4)=vmax(4)+v(i,4)**2
        enddo
      endif
      if(iturbulent.ne.0) then
        do i=nka,nkb
          do j=5,6
            vmax(j)=vmax(j)+v(i,j)**2
          enddo
        enddo
      endif
!     
      return
      end
