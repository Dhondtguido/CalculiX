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
      subroutine edgedivide(nnewnodes,nktet_,ipoed,iexternedg,iedg,
     &     d,h,n,r,iext,jfix)
!     
!     determines the edges to be divided and the number of sub-intervals
!     for each such edge based on the h-field
!     
      implicit none
!     
      integer nnewnodes,i,nktet_,index,ipoed(*),iexternedg(*),iedg(3,*),
     &     n1,n2,n(*),iext,jfix(*)
!     
      real*8 d(*),h(*),h1,h2,r(*)
!
!
!     
      nnewnodes=0
!     
      loop: do i=1,nktet_
        index=ipoed(i)
        do
          if(index.eq.0) cycle loop
!     
!     if iext=0 only internal edges are treated
!     if iext=1 only external edges are treated
!     
          if(((iext.eq.0).and.(iexternedg(index).ne.0)).or.
     &         ((iext.eq.1).and.(iexternedg(index).eq.0))) then
            index=iedg(3,index)
            cycle
          endif
!     
          n1=iedg(1,index)
          n2=iedg(2,index)
!     
!     h1 and h2 are the desired edge length at nodes
!     n1 and n2, respectively
!     
          h1=h(n1)
          h2=h(n2)
!     
          n(index)=int(2*d(index)/(h1+h2)-1)
!     
!     limiting the new nodes to maximum 1 per edge
!     an edge adjacent to domains-not-to-be-refined is not
!     split            
!     
          if((jfix(n1).eq.1).or.(jfix(n2).eq.1)) then
            n(index)=0
          elseif(n(index).gt.1) then
            n(index)=1
          elseif(n(index).lt.0) then
            n(index)=0
          endif
!     
          r(index)=(h2-h1)/(n(index)+2)
!     
          nnewnodes=nnewnodes+n(index)
!     
          index=iedg(3,index)
        enddo
      enddo loop
!     
      return
      end
