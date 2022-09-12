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
      subroutine updatecon(vold,vcon,v,nk,ithermal,iturbulent,
     &     mi,compressible,nka,nkb)
!     
!     updating the conservative variables
!     
      implicit none
!     
      integer iturbulent,mi(*),nk,ithermal(*),i,j,compressible,nka,nkb
!     
      real*8 v(nk,0:mi(2)),vold(0:mi(2),*),vcon(nk,0:mi(2))
!     
!     volumetric energy density
!     
      if(ithermal(1).gt.1) then
        do i=nka,nkb
          vcon(i,0)=vcon(i,0)+v(i,0)
        enddo
      endif
!     
!     volumetric momentum density
!     pressure (liquid) or density (gas)
!     
      do i=nka,nkb
!     
        do j=1,3
          vcon(i,j)=vcon(i,j)+v(i,j)
        enddo
!     
        if(compressible.eq.1) then
!     
!     explicit compressible: v contains the change in
!     density, vcon the density
!     
          vcon(i,4)=vcon(i,4)+v(i,4)
        else
          vold(4,i)=vold(4,i)+v(i,4)
        endif
      enddo
!     
!     volumetric turbulent density
!     
      if(iturbulent.ne.0) then
        do i=nka,nkb
          if(vcon(i,5)+v(i,5).gt.1.d-10) then
            vcon(i,5)=vcon(i,5)+v(i,5)
          else
            v(i,5)=0.d0
          endif
          if(vcon(i,6)+v(i,6).gt.0.d0) then
            vcon(i,6)=vcon(i,6)+v(i,6)
          else
            v(i,6)=0.d0
          endif
        enddo
      endif
!     
      return
      end
      
