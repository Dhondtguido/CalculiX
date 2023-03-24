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
c     Bernhardi start
      subroutine genmodes(i,kon,ipkon,lakon,ne,nk,nk_,co)
!
!     generate nodes for incompatible modes
!
      implicit none
!
      character*8 lakon(*)
!
      real*8 co(3,*),coords(3)
!
      integer i,kon(*),ipkon(*),ne,nope,nopeexp,
     &  nk,nk_,j,indexe,k
!
      indexe=ipkon(i)
!
!     check for elements which may have been deactivated
!
      if(indexe.lt.-1) then
         indexe=-2-indexe
      endif
!
      if(lakon(i)(1:5).eq.'C3D8I')then
         nope=8
         nopeexp=3
      else
         write(*,*) "*ERROR in genmodes: wrong element type, element=",
     &               lakon(i)
         call exit(201)
      endif
!
!     generating additional nodes for the incompatible element. 
!      
!     determining the mean value of the coordinates of the element
!      
      do k=1,3
         coords(k)=0.d0
         do j=1,nope
            coords(k)=coords(k)+co(k,kon(indexe+j))
         enddo
         coords(k)=coords(k)/8.d0
      enddo
!
      do j=1,nopeexp
         nk=nk+1
           if(nk.gt.nk_) then
              write(*,*) '*ERROR in genmodes: increase nk_'
              call exit(201)
           endif
         kon(indexe+nope+j)=nk
         do k=1,3
            co(k,nk)=coords(k)
         enddo
      enddo
!
      return
      end
c     Bernhardi end

