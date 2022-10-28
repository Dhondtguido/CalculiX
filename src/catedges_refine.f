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
      subroutine catedges_refine(netet_,iedg,kontet,ipoed,ifreeed,
     &  iedtet,ipoeled,ieled,ifreele)
!
!     catalogueing the edges of the tetrahedral elements of the mesh
!
      implicit none
!
      integer i,j,netet_,iedg(3,*),kontet(4,*),nodee(2),ie(2,6),n,
     &  node,index,ipoed(*),ifreeed,iedtet(6,*),
     &  ieled(2,*),iedge,ifreele,ipoeled(*)
!
!     nodes belonging to the six edges
!
      data ie /1,2,2,3,1,3,1,4,2,4,3,4/
!
!     initialization of ipoed and iedg
!
      do i=1,6*netet_
         iedg(3,i)=i+1
      enddo
      iedg(3,6*netet_)=0
!
!     setting up the edge databank
!
      do j=1,netet_
         if(kontet(1,j).eq.0) cycle
!
!        loop over all edges per tet
!
         do i=1,6
            nodee(1)=kontet(ie(1,i),j)
            nodee(2)=kontet(ie(2,i),j)
            n=2
c            kflag=1
            call insertsorti(nodee,n)
c            call isortii(nodee,idum,n,kflag)
!
!           check whether edge is already catalogued
!
            node=nodee(1)
            index=ipoed(node)
!
            do
               if(index.eq.0) exit
               if(iedg(2,index).eq.nodee(2)) exit
               index=iedg(3,index)
            enddo
!
            if(index.eq.0) then
               index=ifreeed
               ifreeed=iedg(3,ifreeed)
               if(ifreeed.eq.0) then
                  write(*,*) '*ERROR in catedges: increase'
                  write(*,*) '       the dimension of iedg'
                  call exit(201)
               endif
               iedg(1,index)=nodee(1)
               iedg(2,index)=nodee(2)
               iedg(3,index)=ipoed(node)
               ipoed(node)=index
            endif
            iedtet(i,j)=index
         enddo
      enddo
!
!     initialization of ieled
!
      do i=1,6*netet_
         ieled(2,i)=i+1
      enddo
      ieled(2,6*netet_)=0
!
!     generating the element per edge relationship
!
      do j=1,netet_
         if(kontet(1,j).eq.0) cycle
         do i=1,6
            iedge=iedtet(i,j)
            index=ifreele
            ieled(1,index)=j
            ifreele=ieled(2,index)
            if(ifreele.eq.0) then
               write(*,*) '*ERROR in catedges: increase the'
               write(*,*) '       dimension of ieln'
               call exit(201)
            endif
            ieled(2,index)=ipoeled(iedge)
            ipoeled(iedge)=index
         enddo
      enddo
!
      return
      end

