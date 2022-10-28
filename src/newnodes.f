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
      subroutine newnodes(nktet_,ipoed,n,iedg,h,d,r,
     &     conewnodes,cotet,ibasenewnodes,ipoeled,ieled,doubleglob,
     &     integerglob,nnewnodes,iedgnewnodes,hnewnodes,n1newnodes,
     &     n2newnodes)
!     
!     determines the edges to be divided and the number of sub-intervals
!     for each such edge based on the h-field
!     
      implicit none
!     
      integer indexn,i,nktet_,ipoed(*),n(*),iedg(3,*),
     &     n1,n2,k,ibasenewnodes(*),ipoeled(*),index,j,ieled(2,*),
     &     nnewnodes,n1newnodes(*),n2newnodes(*),
     &     integerglob(*),nktri,netet,ne,nkon,nfaces,nfield,
     &     nselect,iselect(1),imastset,istartset(1),iendset(1),
     &     ialset(1),nterms,konl(20),iedgnewnodes(*),nelem,loopa
!     
      real*8 alphasum,h(*),h1,d(*),r(*),scale,hscale,rscale,
     &     conewnodes(3,*),cotet(3,*),dist,
     &     doubleglob(*),coords(3),ratio(20),hnewnodes(*)
!     
      indexn=0
!     
      loop: do i=1,nktet_
      index=ipoed(i)
      do
        if(index.eq.0) cycle loop
!     
!     no subdivision or external edge: cycle
!     
        if(n(index).eq.0) then
          index=iedg(3,index)
          cycle
        endif
!     
!     scale subdivision lengths: 
!     
        alphasum=0.d0
        n1=iedg(1,index)
        n2=iedg(2,index)
        h1=h(n1)
!     
!     the subintervals are numbered 0,...,n(index)
!     the length of the subintervals is progressively
!     increased by r(index) (may be negative)
!     
        do k=0,n(index)
          alphasum=alphasum+h1+(k+1)*r(index)
        enddo
        scale=d(index)/alphasum
!     
!     loop over all new positions: numbered from
!     1,...,n(index)
!     
        alphasum=0.d0
        hscale=scale*h1
        rscale=scale*r(index)
        do j=1,n(index)
          alphasum=alphasum+hscale+j*rscale
          do k=1,3
            conewnodes(k,indexn+j)=cotet(k,n1)+alphasum*
     &           (cotet(k,n2)-cotet(k,n1))/d(index)
          enddo
!     
!     determine an element to which the new node belongs
!     
          ibasenewnodes(indexn+j)=ieled(1,ipoeled(index))
          iedgnewnodes(indexn+j)=index
          n1newnodes(indexn+j)=n1
          n2newnodes(indexn+j)=n2
        enddo
        indexn=indexn+n(index)
        index=iedg(3,index)
      enddo
      enddo loop
!     
!     initializing fields
!     
      nktri=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=1
      nselect=1
      iselect(1)=1
      imastset=0
!     
      do i=1,nnewnodes
        coords(1)=conewnodes(1,i)
        coords(2)=conewnodes(2,i)
        coords(3)=conewnodes(3,i)
        loopa=8
        call basis(doubleglob(1),doubleglob(netet+1),
     &       doubleglob(2*netet+1),
     &       doubleglob(3*netet+1),doubleglob(4*netet+1),
     &       doubleglob(5*netet+1),integerglob(6),
     &       integerglob(netet+6),
     &       integerglob(2*netet+6),doubleglob(6*netet+1),
     &       integerglob(3*netet+6),nktri,netet,
     &       doubleglob(4*nfaces+6*netet+1),nfield,
     &       doubleglob(nktri+4*nfaces+6*netet+1),
     &       integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &       integerglob(2*ne+7*netet+6),
     &       integerglob(nkon+2*ne+7*netet+6),
     &       coords(1),coords(2),coords(3),hnewnodes(i),
     &       ratio,iselect,nselect,
     &       istartset,iendset,ialset,imastset,
     &       integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,
     &       loopa,dist)
!     
      enddo
!     
      return
      end
