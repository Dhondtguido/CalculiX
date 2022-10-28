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
      subroutine searchmidneigh(inn,iponn,nktet_,iexternedg,
     &     ipoed,iedg,ipoeled,ieled,ifreenn,iedgmid,iedtet)
!     
!     look for the midnodes' neighbors
!     
      implicit none
!     
      integer nktet_,inn(2,*),iponn(*),iexternedg(*),
     &     ipoed(*),iedg(3,*),i,j,k,iedge,ielem,ieled(2,*),
     &     indexe,ifreenn,ipoeled(*),ifreennstart,neigh,iedgmid(*),
     &     iedtet(6,*),node
!     
!     determining neighboring elements
!     
      do i=1,nktet_
!
!       loop over all edges for which i is the lowest vertex node
!
        iedge=ipoed(i)
        do
          if(iedge.eq.0) exit
          if(iexternedg(iedge).eq.0) then
            node=iedgmid(iedge)
            ifreennstart=ifreenn
!
!           loop over all elements to which edge "iedge" belongs
!
            indexe=ipoeled(iedge)
            do
              if(indexe.eq.0) exit
              ielem=ieled(1,indexe)
              loop: do j=1,6
                neigh=iedgmid(iedtet(j,ielem))
                if(neigh.ne.node) then
!     
!                 check whether the middle node has already been
!                 catalogued
!     
                  do k=ifreennstart,ifreenn-1
                    if(inn(1,k).eq.neigh) cycle loop
                  enddo
                  inn(1,ifreenn)=neigh
                  inn(2,ifreenn)=iponn(node)
                  iponn(node)=ifreenn
                  ifreenn=ifreenn+1
                endif
              enddo loop
              indexe=ieled(2,indexe)
            enddo
          endif
          iedge=iedg(3,iedge)
        enddo
      enddo
!     
      return
      end
