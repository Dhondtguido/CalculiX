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
      subroutine catedges_crackprop(ipoed,iedg,ntri,ieled,kontri,
     &     nedg,ier)
!
!     catalogueing the edges of the tetrahedral elements of the mesh
!
      implicit none
!
      integer ipoed(*),iedg(3,*),ntri,ie(2,3),ieled(2,*),nodee(2),
     &     n,node,index,kontri(3,*),nedg,ifreeed,i,j,ier
!
!     nodes belonging to the six edges
!
      data ie /1,2,2,3,3,1/
!
!     initialization of iedg
!
      do i=1,3*ntri
         iedg(3,i)=i+1
      enddo
      iedg(3,3*ntri)=0
!
!     setting up the edge databank
!
      ifreeed=1
!
      do j=1,ntri
!
!        loop over all edges per triangle
!
         loop:do i=1,3
            nodee(1)=kontri(ie(1,i),j)
            nodee(2)=kontri(ie(2,i),j)
            n=2
            call insertsorti(nodee,n)
!
!           check whether edge is already catalogued
!
            node=nodee(1)
            index=ipoed(node)
!
            do
               if(index.eq.0) exit
               if(iedg(2,index).eq.nodee(2)) then
                  ieled(2,index)=j
                  cycle loop
               endif
               index=iedg(3,index)
            enddo
!
            index=ifreeed
            ifreeed=iedg(3,ifreeed)
            if(ifreeed.eq.0) then
               write(*,*) '*ERROR in catedges_crackprop: increase'
               write(*,*) '       the dimension of iedg'
               ier=1
c               call exit(201)
            endif
            iedg(1,index)=nodee(1)
            iedg(2,index)=nodee(2)
            iedg(3,index)=ipoed(node)
            ipoed(node)=index
            ieled(1,index)=j
         enddo loop
      enddo
!
      nedg=ifreeed-1
!
      return
      end

