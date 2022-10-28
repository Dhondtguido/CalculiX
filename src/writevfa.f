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
      subroutine writevfa(vfa,nface,nactdohinv,ielfa)
!
!     writing facial values
!
      implicit none
!
      integer nface,nactdohinv(*),ielfa(4,*),i,iel
!
      real*8 vfa(0:7,*)
!
      do i=1,nface
         if(ielfa(2,i).ge.0) cycle
         iel=ielfa(1,i)
         iel=nactdohinv(iel)
         write(*,*) 'writevfa ',iel,ielfa(4,i),
     &      vfa(0,i),vfa(1,i),vfa(2,i),vfa(3,i)
      enddo
!     
      return
      end
