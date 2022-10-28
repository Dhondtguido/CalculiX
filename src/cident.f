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
!
!     identifies the position id of px in an ordered array
!     x of integers; 
!
!     id is such that x(id).le.px and x(id+1).gt.px
!
      subroutine cident(x,px,n,id)
      implicit none
      character*81 x,px
      integer n,id,n2,m
      dimension x(n)
      id=0
      if(n.eq.0) return
      n2=n+1
      do
         m=(n2+id)/2
         if(px.ge.x(m)) then
            id=m
         else
            n2=m
         endif
         if((n2-id).eq.1) return
      enddo
      end
