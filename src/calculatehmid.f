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
      subroutine calculatehmid(nktet_,h,ipoed,iedg,iedgmid)
!     
!     calculating the desired size of h in all midnodes
!     
      implicit none
!     
      integer nktet_,ipoed(*),iedg(3,*),iedgmid(*),node1,node2,
     &     nodem,i,index
!     
      real*8 h(*)
!     
!     the desired edge length at the midnodes is the mean of the
!     desired length at its neighbors
!     
      do i=1,nktet_
        index=ipoed(i)
!     
        do
          if(index.eq.0) exit
!
          node1=iedg(1,index)
          node2=iedg(2,index)
          nodem=iedgmid(index)
!
          h(nodem)=(h(node1)+h(node2))/2.d0
!     
          index=iedg(3,index)
        enddo
      enddo
!     
      return
      end
