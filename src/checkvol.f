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
      subroutine checkvol(cotet,node,pnew,ipoeln,ieln,kontet,c1,jflag,
     &     ibadnodes,nbadnodes,iwrite)
!
!     check the volume of all elements belonging to "node" if the
!     coordinates of "node" are changed from cotet(1..3,node) to
!     pnew(1..3); if the volume is strictly positive, the coordinates are
!     changed, if not they are left unchanged
!
!     used in projectvertexnodes.f
!
      implicit none
!
      integer node,ipoeln(*),ieln(2,*),kontet(4,*),i,ielem,indexe,jflag,
     &     j,n,ibadnodes(*),nbadnodes,id,listed,m,iwrite,neigh
!
      real*8 cotet(3,*),pnew(3),pold(3),volume,c1,c2
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
      n=3
      loop: do j=1,n
        indexe=ipoeln(node)
        do
          if(indexe.eq.0) exit loop
          ielem=ieln(1,indexe)
          call calcvol(kontet(1,ielem),kontet(2,ielem),kontet(3,ielem),
     &         kontet(4,ielem),cotet,volume)
          if(volume.le.0.d0) then
            exit
          endif
          indexe=ieln(2,indexe)
        enddo
        if(j.eq.n) then
          do i=1,3
            cotet(i,node)=pold(i)
          enddo
          exit
        endif
        c2=c2/2.d0
        do i=1,3
          cotet(i,node)=c2*pnew(i)+(1.d0-c2)*pold(i)
        enddo
      enddo loop
!
!     if node was not completely projected: list as bad node
!     -> application of fminsirefine for all neighboring vertex
!        nodes which are not on the free surface      
!
      if(j.ne.1) then
        indexe=ipoeln(node)
        do
          if(indexe.eq.0) exit
          ielem=ieln(1,indexe)
          do i=1,4
            neigh=kontet(i,ielem)
            listed=0
            call nident(ibadnodes,neigh,nbadnodes,id)
            if(id.gt.0) then
              if(ibadnodes(id).eq.neigh) then
                listed=1
              endif
            endif
            if(listed.eq.0) then
              nbadnodes=nbadnodes+1
              do m=nbadnodes,id+2,-1
                ibadnodes(m)=ibadnodes(m-1)
              enddo
              ibadnodes(id+1)=neigh
            endif
          enddo
          indexe=ieln(2,indexe)
        enddo
      endif
!
!     if jflag=1 and j>1: last chance for a complete projection
!     did not work
!
      if((jflag.eq.1).and.(j.ne.1)) then
        write(*,*) '*WARNING in checkvol: projection of vertex node ',
     &       node
        write(*,*) '         had to be reduced to keep the adjacent'
        write(*,*) '         elements regular'
        write(40,*) node
        iwrite=1
      endif
!      
      return
      end
