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
      subroutine removetet_refine(kontet,ifatet,ifreetet,ifac,itetfa,
     &  ifreefa,ipofa,ielement,ipoeln,ieln,ifreeln,
     &  ipoeled,ieled,ifreele,iedtet,ipoed,iedg,ifreeed,
     &  iexternfa,iexternedg)
!
!     removes a tet and updates the data base
!     equivalent to subroutine removetet, apart from the extra
!     update of the ieln data base
!
      implicit none
!
      integer nodes(4),kontet(4,*),ifatet(4,*),ifac(4,*),
     &  itetfa(2,*),ipofa(*),ipoeln(*),ieln(2,*),ifreeln,
     &  ifreetet,ifreefa,iedge,
     &  ipoeled(*),ieled(2,*),iedtet(6,*),ifreele,ipoed(*),
     &  iedg(3,*),ifreeed,iexternfa(*),iexternedg(*),
     &  index,i,iface,node,indexold,ielement
!
!
!
      do i=1,4
         nodes(i)=kontet(i,ielement)
      enddo
!
!     for each face belonging to the element:
!     - remove the element from database itetfa
!     - if the face does not belong to any other element,
!       remove the face from database ifac
!
      do i=1,4
         iface=abs(ifatet(i,ielement))
!
         if(itetfa(1,iface).eq.ielement) then
            itetfa(1,iface)=itetfa(2,iface)
            itetfa(2,iface)=0
         else
            itetfa(2,iface)=0
         endif
!
!        remove the face from database ifac if it does
!        not belong to any element any more (except if
!        external face)
!
         if((itetfa(1,iface).eq.0).and.(iexternfa(iface).le.0)) then
            node=ifac(1,iface)
            index=ipofa(node)
            if(index.eq.iface) then
               ipofa(node)=ifac(4,index)
               ifac(4,index)=ifreefa
               ifreefa=index
            else
               do
                  indexold=index
                  index=ifac(4,index)
                  if(index.eq.0) then
                     write(*,*) '*ERROR in removetet_refine: face to be'
                     write(*,*) '       deleted is not catalogued'
                     write(*,*) '       in field ifac'
                     call exit(201)
                  endif
                  if(index.eq.iface) then
                     ifac(4,indexold)=ifac(4,index)
                     ifac(4,index)=ifreefa
                     ifreefa=index
                     exit
                  endif
               enddo
            endif
         endif
      enddo
!
!     update the element per node information
!
      do i=1,4
         index=ipoeln(nodes(i))
         indexold=0
         do
            if(ieln(1,index).eq.ielement) exit
            indexold=index
            index=ieln(2,index)
            if(index.eq.0) then
               write(*,*) '*ERROR in removetet_refine: error in'
               write(*,*) '       ieln database'
               call exit(201)
            endif
         enddo
         if(indexold.eq.0) then
            ipoeln(nodes(i))=ieln(2,index)
         else
            ieln(2,indexold)=ieln(2,index)
         endif
         ieln(2,index)=ifreeln
         ifreeln=index
      enddo
!
!     for each edge belonging to the element:
!     - remove the element from database ieled
!     - if the edge does not belong to any other element,
!       remove the edge from database iedg (except if
!       external edge)
!
      do i=1,6
         iedge=iedtet(i,ielement)
         index=ipoeled(iedge)
         indexold=0
         do
            if(ieled(1,index).eq.ielement) exit
            indexold=index
            index=ieled(2,index)
            if(index.eq.0) then
               write(*,*) '*ERROR in removetet_refine: error in'
               write(*,*) '       ieled database'
               call exit(201)
            endif
         enddo
         if(indexold.eq.0) then
            ipoeled(iedge)=ieled(2,index)
         else
            ieled(2,indexold)=ieled(2,index)
         endif
         ieled(2,index)=ifreele
         ifreele=index
!
         if((ipoeled(iedge).eq.0).and.(iexternedg(iedge).eq.0)) then
            node=iedg(1,iedge)
            index=ipoed(node)
            if(index.eq.iedge) then
               ipoed(node)=iedg(3,index)
               iedg(3,index)=ifreeed
               ifreeed=index
            else
               do
                  indexold=index
                  index=iedg(3,index)
                  if(index.eq.0) then
                     write(*,*) '*ERROR in removetet_refine'
                     write(*,*) '       database iedg corrupted'
                     call exit(201)
                  endif
                  if(index.eq.iedge) then
                     iedg(3,indexold)=iedg(3,index)
                     iedg(3,index)=ifreeed
                     ifreeed=index
                     exit
                  endif
               enddo
            endif
c            iedg(1,iedge)=0
c            iedg(2,iedge)=0
         endif
      enddo
!
!     remove the element
!
      kontet(1,ielement)=0
c      kontet(2,ielement)=0
c      kontet(3,ielement)=0
      kontet(4,ielement)=ifreetet
      ifreetet=ielement
!
      return
      end
            
