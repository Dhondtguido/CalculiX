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
      subroutine genmidnodes(nktet_,ipoed,iedgmid,iexternedg,
     &     iedgext,cotet,nktet,iedg,jfix,ipoeled,ieled,kontet,
     &     iedtet,iwrite)
!     
!     generates midnodes in the middle between the neighboring vertex
!     nodes
!     
      implicit none
!     
      integer nktet_,ipoed(*),index,i,k,iedgmid(*),iexternedg(*),
     &     iedgext(3,*),node1,node2,iedg(3,*),nktet,jfix(*),ichange,
     &     ielem,iflag,indexe,iwrite,m,ipoeled(*),ieled(2,*),
     &     kontet(4,*),iedtet(6,*)
!     
      real*8 cotet(3,*),xi,et,ze,shp(4,10),xsj,xl(3,10)
!
      include "gauss.f"
!
      iflag=2
!     
!     generate middle nodes 
!     
      do i=1,nktet_
        index=ipoed(i)
!     
        do
          if(index.eq.0) exit
!
          node1=iedg(1,index)
          node2=iedg(2,index)
!     
          if((jfix(node1).eq.1).and.(jfix(node2).eq.1).and.
     &         (iexternedg(index).gt.0)) then
            iedgmid(index)=iedgext(2,iexternedg(index))
          else
            nktet=nktet+1
            iedgmid(index)=nktet
            do k=1,3
              cotet(k,nktet)=
     &             (cotet(k,node1)+cotet(k,node2))/2.d0
            enddo
          endif
!     
          index=iedg(3,index)
        enddo
      enddo
!     
!     check quality of elements (may be bad due to the midnodes on
!     the fixed edges; these are not necessarily in the middle between
!     the vertex nodes)
!
      do
        ichange=0
        loop: do i=1,nktet_
          index=ipoed(i)
!     
          do
            if(index.eq.0) exit
!     
            node1=iedg(1,index)
            node2=iedg(2,index)
!     
            if((jfix(node1).eq.1).and.(jfix(node2).eq.1).and.
     &           (iexternedg(index).gt.0)) then
!     
!     check whether all elements in the shell of the edge
!     have positive Jacobians
!     
              indexe=ipoeled(index)
              do
                if(indexe.eq.0) exit
                ielem=ieled(1,indexe)
!     
!     vertex nodes
!     
                do k=1,4
                  do m=1,3
                    xl(m,k)=cotet(m,kontet(k,ielem))
                  enddo
                enddo
!     
!     middle nodes
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
                    do m=1,3
                      cotet(m,iedgmid(index))=
     &                     (cotet(m,node1)+cotet(m,node2))/2.d0
                    enddo
                    write(*,*) '*WARNING in genmidnodes: '
                    write(*,*) '         fixed midnode ',iedgmid(index)
                    write(*,*)
     &                   '         had to be moved into the middle'
                    write(*,*)
     &                   '         of its neighboring vertex nodes'
                    write(*,*) '         to keep the adjacent'
                    write(*,*) '         elements regular'
                    write(*,*)
                    write(40,*) iedgmid(index)
                    iwrite=1
                    ichange=1
                    cycle loop
                  endif
                enddo
                indexe=ieled(2,indexe)
              enddo
            endif
!     
            index=iedg(3,index)
          enddo
        enddo loop
        if(ichange.eq.0) exit
      enddo
!     
      return
      end
