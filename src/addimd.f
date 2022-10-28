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
      subroutine addimd(imd,nmd,node)
!
!     adds entity "node" to field imd. imd contains the
!     entities selected by the user in which results are to be
!     calculated in a modal dynamics calculation
!
      implicit none
!
      integer imd(*),nmd,node,id,l
!
      call nident(imd,node,nmd,id)
      do
         if(id.gt.0) then
            if(imd(id).eq.node)exit
         endif
         nmd=nmd+1
         do l=nmd,id+2,-1
            imd(l)=imd(l-1)
         enddo
         imd(id+1)=node
         exit
      enddo
!
      return
      end
