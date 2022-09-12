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
      subroutine midexternaledges(iexternedg,nexternedg,iedgext,
     &   ifreeed,ieled,ipoeled,iedg,iedtet,kontetor)
!
!     stores the nodes belonging to the external edges of the
!     unrefined mesh
!
      implicit none
!
      integer nexternedg,iexternedg(*),iedgext(3,*),i,j,ifreeed,
     &  iel,ieled(2,*),ipoeled(*),iedg(3,*),iedtet(6,*),kontetor(6,*)
!
!
!
      nexternedg=0
!
!     loop over all edges. For the unrefined mesh the edges are stored
!     in iedg in a consecutive order. 
!     An edge i is external if iexternedg(i)!=0; at exit the value
!     of iexternedg(i) is the number of the external edge in field
!     iedgext
!
      do i=1,ifreeed-1
         if(iexternedg(i).ne.0) then
            nexternedg=nexternedg+1
            iexternedg(i)=nexternedg
!
!           end nodes of the edge
!
            iedgext(1,nexternedg)=iedg(1,i)
            iedgext(3,nexternedg)=iedg(2,i)
!
!           any element containing the edge          
!
            iel=ieled(1,ipoeled(i))
!
!           recovering the middle node of the edge
!           if the edge is linear, this node is zero
!
            do j=1,6
               if(iedtet(j,iel).eq.i) then
                  iedgext(2,nexternedg)=kontetor(j,iel)
                  exit
               endif
            enddo
         endif
      enddo
!
      return
      end
