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
      subroutine determineextern(ifac,itetfa,iedg,ipoed,
     &     iexternedg,iexternfa,iexternnode,nktet_,ipofa)
!
!     determines which faces, edges and nodes are external
!     
!     a face is external if it belongs to only one element
!     an edge is external if it belongs to an external face    
!     a node is external if it belongs to an external face 
!
      implicit none
!
      integer i,nktet_,index,indexe,ipofa(*),itetfa(2,*),ifac(4,*),
     &  iexternfa(*),iexternedg(*),iexternnode(*),n1,n2,n3,ipoed(*),
     &  iedg(3,*)
!
      loop: do i=1,nktet_
         index=ipofa(i)
         do
            if(index.eq.0) cycle loop
            if(itetfa(2,index).ne.0) then
               index=ifac(4,index)
               cycle
            endif
!
!           index is an external face
!
            iexternfa(index)=1
!
!           checking the nodes belonging to the face
!
            n1=ifac(1,index)
            n2=ifac(2,index)
            n3=ifac(3,index)
!
            iexternnode(n1)=1
            iexternnode(n2)=1
            iexternnode(n3)=1
!
!           checking the edges belonging to the face
!           edge n1-n2
!
            if(n1.lt.n2) then
               indexe=ipoed(n1)
               do
                  if(indexe.eq.0) exit
                  if(iedg(2,indexe).eq.n2) then
                     iexternedg(indexe)=indexe
c                     write(*,*) 'determineextern',iedg(1,indexe),
c     &                   iedg(2,indexe),iexternedg(indexe)
                     exit
                  endif
                  indexe=iedg(3,indexe)
               enddo
            else
               indexe=ipoed(n2)
               do
                  if(indexe.eq.0) exit
                  if(iedg(2,indexe).eq.n1) then
                     iexternedg(indexe)=indexe
c                     write(*,*) 'determineextern',iedg(1,indexe),
c     &                   iedg(2,indexe),iexternedg(indexe)
                     exit
                  endif
                  indexe=iedg(3,indexe)
               enddo
            endif
!
!           checking the edges belonging to the face
!           edge n2-n3
!
            if(n2.lt.n3) then
               indexe=ipoed(n2)
               do
                  if(indexe.eq.0) exit
                  if(iedg(2,indexe).eq.n3) then
                     iexternedg(indexe)=indexe
c                     write(*,*) 'determineextern',iedg(1,indexe),
c     &                   iedg(2,indexe),iexternedg(indexe)
                     exit
                  endif
                  indexe=iedg(3,indexe)
               enddo
            else
               indexe=ipoed(n3)
               do
                  if(indexe.eq.0) exit
                  if(iedg(2,indexe).eq.n2) then
                     iexternedg(indexe)=indexe
c                     write(*,*) 'determineextern',iedg(1,indexe),
c     &                   iedg(2,indexe),iexternedg(indexe)
                     exit
                  endif
                  indexe=iedg(3,indexe)
               enddo
            endif
!
!           checking the edges belonging to the face
!           edge n3-n1
!
            if(n3.lt.n1) then
               indexe=ipoed(n3)
               do
                  if(indexe.eq.0) exit
                  if(iedg(2,indexe).eq.n1) then
                     iexternedg(indexe)=indexe
c                     write(*,*) 'determineextern',iedg(1,indexe),
c     &                   iedg(2,indexe),iexternedg(indexe)
                     exit
                  endif
                  indexe=iedg(3,indexe)
               enddo
            else
               indexe=ipoed(n1)
               do
                  if(indexe.eq.0) exit
                  if(iedg(2,indexe).eq.n3) then
                     iexternedg(indexe)=indexe
c                     write(*,*) 'determineextern',iedg(1,indexe),
c     &                   iedg(2,indexe),iexternedg(indexe)
                     exit
                  endif
                  indexe=iedg(3,indexe)
               enddo
            endif
!
            index=ifac(4,index)
            cycle
         enddo
      enddo loop
!
      return
      end
