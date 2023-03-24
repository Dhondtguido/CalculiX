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
      subroutine gasmechbc(vold,nload,sideload,
     &     nelemload,xload,mi)
!     
      implicit none
!     
      character*20 sideload(*) 
!     
      integer i,nload,node,nelemload(2,*),mi(*)
!     
      real*8 vold(0:mi(2),*),xload(2,*)
!
!     updating the boudary conditions in a mechanical
!     calculation coming from a previous thermal calculation
!     
!     updating the pressure boundary conditions
!     
      do i=1,nload
         if(sideload(i)(3:4).eq.'NP') then
            node=nelemload(2,i)
            xload(1,i)=vold(2,node)
         endif
      enddo
!      
      return
      end








