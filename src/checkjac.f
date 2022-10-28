!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine checkjac(cotet,node,pnew,kontet,c1,jflag,
     &     iedtet,iedgmid,ipoeled,ieled,iedge,ibadnodes,nbadnodes,
     &     iwrite)
!
!     check the Jacobian in all integration points of all elements
!     belonging to midnode "node" if the coordinates of "node" are
!     changed from cotet(1..3,node) to pnew(1..3); if the Jacobian
!     is strictly positive, the coordinates are changed, if not, they
!     are left unchanged
!     "node" belongs to edge "iedge"      
!     
      implicit none
!
      integer node,kontet(4,*),i,ielem,indexe,iflag,jflag,iwrite,
     &     j,n,iedtet(6,*),iedgmid(*),ipoeled(*),ieled(2,*),
     &     iedge,k,m,ibadnodes(*),nbadnodes,id,listed,neighedge
!
      real*8 cotet(3,*),pnew(3),pold(3),xsj,c1,c2,xi,et,ze,xl(3,10),
     &     shp(4,10)
!
      include "gauss.f"
!
      iflag=2
!
!     backup of the old coordinates
!
      do i=1,3
        pold(i)=cotet(i,node)
      enddo
!
!     storing the new coordinates in cotet
!
      c2=c1
      do i=1,3
        cotet(i,node)=c2*pnew(i)+(1.d0-c2)*pold(i)
      enddo
!
!     trial loops; if bad elements the projection is reduced and the
!     elements are rechecked
!
      n=3
      loop1: do j=1,n
      indexe=ipoeled(iedge)
!
!       loop over all elements belonging to edge "iedge"
!
        loop2: do
          if(indexe.eq.0) exit loop1
          ielem=ieled(1,indexe)
!
!         vertex nodes
!
          do k=1,4
            do m=1,3
              xl(m,k)=cotet(m,kontet(k,ielem))
            enddo
          enddo
!
!         middle nodes
!
          do k=1,6
            do m=1,3
              xl(m,k+4)=cotet(m,iedgmid(iedtet(k,ielem)))
            enddo
          enddo
!
          do k=1,4
            xi=gauss3d5(1,k)
            et=gauss3d5(2,k)
            ze=gauss3d5(3,k)
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            if(xsj.le.0.d0) then
c              write(*,*) 'checkjac bad element ',ielem,' node ',node
              exit loop2
            endif
          enddo
          indexe=ieled(2,indexe)
        enddo loop2
!
!       last loop: revert to non-projected coordinates
!
        if(j.eq.n) then
          do i=1,3
            cotet(i,node)=pold(i)
          enddo
          exit
        endif
!
!       not last loop: reduce the projection and retry
!
        c2=c2/2.d0
        do i=1,3
          cotet(i,node)=c2*pnew(i)+(1.d0-c2)*pold(i)
        enddo
      enddo loop1
!
!     if node was not completely projected: list all edges of
!     neighboring elements (corresponds 1-to-1 to midnodes)
!     -> application of fminsirefine for all neighboring edges/midnodes
!        which are not on the free surface      
!
      if(j.ne.1) then
        indexe=ipoeled(iedge)
        do
          if(indexe.eq.0) exit
!
!         adjacent element
!
          ielem=ieled(1,indexe)
          do k=1,6
            neighedge=iedtet(k,ielem)
!
!           edge of adjacent element
!
            listed=0
            call nident(ibadnodes,neighedge,nbadnodes,id)
            if(id.gt.0) then
              if(ibadnodes(id).eq.neighedge) then
                listed=1
              endif
            endif
            if(listed.eq.0) then
              nbadnodes=nbadnodes+1
              do m=nbadnodes,id+2,-1
                ibadnodes(m)=ibadnodes(m-1)
              enddo
              ibadnodes(id+1)=neighedge
            endif
          enddo
          indexe=ieled(2,indexe)
        enddo
      endif
!
!     if jflag=1 and j>1: last chance for a complete projection
!     did not work
!
      if((jflag.eq.1).and.(j.ne.1)) then
        write(*,*) '*WARNING in checkjac: projection of midnode ',node
        write(*,*) '         had to be reduced to keep the adjacent'
        write(*,*) '         elements regular'
        write(*,*)
        write(40,*) node
        iwrite=1
      endif
!      
      return
      end
