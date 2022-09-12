!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine resultsp(nk,nactdoh,v,sol,mi)
!
!     calculates the pressure correction (STEP 2) in the nodes
!
      implicit none
!
      integer nk,nactdoh(*),i,mi(*)
!
      real*8 sol(*),v(nk,0:mi(2))
!
!     extracting the pressure correction from the solution
!
      do i=1,nk
        if(nactdoh(i).gt.0) then
          v(i,4)=sol(nactdoh(i))
        else
          v(i,4)=0.d0
        endif
      enddo
!
      return
      end
