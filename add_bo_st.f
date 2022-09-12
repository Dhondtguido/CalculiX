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
      subroutine add_bo_st(au,jq,irow,i,j,value)
!
!     stores the boundary stiffness coefficient (i,j) with value "value"
!     in the stiffness matrix stored in spare matrix format
!
      implicit none
!
      integer jq(*),irow(*),i,j,ipointer,id
!
      real*8 au(*),value
!
!
!
      call nident(irow(jq(j)),i,jq(j+1)-jq(j),id)
!
      ipointer=jq(j)+id-1
!
      if(irow(ipointer).ne.i) then
c         write(*,*) i,j,ipointer,irow(ipointer)
         write(*,*) '*ERROR in add_bo_st: coefficient should be 0'
         call exit(201)
      else
         au(ipointer)=au(ipointer)+value
      endif
!
      return
      end













