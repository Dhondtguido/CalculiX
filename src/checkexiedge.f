!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998 Guido Dhondt
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
      subroutine checkexiedge(n1,n2,ipoed,iedg,node)
!
!     check whether the edge consisting of nodes n1 and n2 (n1<n2)
!     still exists      
!
      implicit none
!     
      integer n1,n2,ipoed(*),iedg(3,*),node,index
!     
!     
!
      index=ipoed(n1)
      do
        if(index.eq.0) then
          node=0
c          write(*,*) '*INFO in checkexiedge'
c          write(*,*) '      the edge to be split and consisting of'
c          write(*,*) '      the nodes ',n1,' and ',n2,' does not exist'
c          write(*,*) '      any more; the next edge is addressed'
          exit
        endif
        if(iedg(2,index).eq.n2) exit
        index=iedg(3,index)
      enddo
!     
      return
      end
