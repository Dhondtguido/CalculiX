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
      subroutine genmpc(inodestet,nnodestet,co,doubleglob,integerglob,
     &     ipompc,nodempc,coefmpc,nmpc,nmpc_,labmpc,mpcfree,ikmpc,
     &     ilmpc)
!     
!     generating MPC's connecting the nodes in the old tet mesh in
!     which SPC's or point forces were defined with the nodes in the
!     refined tet mesh
!     
      implicit none
!     
      character*20 labmpc(*)
!     
      integer inodestet(*),nnodestet,integerglob(*),nktet,netet,ne,nkon,
     &     nfaces,nfield,nselect,imastset,iselect(1),nterms,
     &     nelem,ialset(1),mpcfreeold,
     &     iendset(1),istartset(1),konl(20),loopa,
     &     node,i,j,k,m,ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &     ikmpc(*),ilmpc(*),idof,id
!     
      real*8 co(3,*),doubleglob(*),coords(3),ratio(20),value,
     &     coefmpc(*),dist
!     
      nktet=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=0
      nselect=0
      imastset=0
      loopa=8
!     
      loop1:do i=1,nnodestet
        node =inodestet(i)
!     
        do j=1,3
          coords(j)=co(j,node)
        enddo
!     
        call basis(doubleglob(1),doubleglob(netet+1),
     &       doubleglob(2*netet+1),
     &       doubleglob(3*netet+1),doubleglob(4*netet+1),
     &       doubleglob(5*netet+1),integerglob(6),integerglob(netet+6),
     &       integerglob(2*netet+6),doubleglob(6*netet+1),
     &       integerglob(3*netet+6),nktet,netet,
     &       doubleglob(4*nfaces+6*netet+1),nfield,
     &       doubleglob(13*nktet+4*nfaces+6*netet+1),
     &       integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &       integerglob(2*ne+7*netet+6),
     &       integerglob(nkon+2*ne+7*netet+6),
     &       coords(1),coords(2),coords(3),value,ratio,iselect,nselect,
     &       istartset,iendset,ialset,imastset,
     &       integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,loopa,
     &       dist)
!     
!     if the old node number was kept in the new mesh it means that
!     this node was not moved and no MPC's are needed
!     
        do j=1,nterms
          if(konl(j).eq.node) then
c            write(*,*) '*INFO in genmpc: no MPC generated for node',node
            cycle loop1
          endif
        enddo
!     
!     check whether SPC was applied in the node; if not, the
!     old node is used as dependent node
!     
        do j=1,3
          idof=8*(node-1)+j
          call nident(ikmpc,idof,nmpc,id)
!     
!     temporarily more than one MPC with the same idof is allowed
!     
c          if(id.gt.0) then
c            if(ikmpc(id).eq.idof) then
c              write(*,*) '*ERROR in genmpc: a MPC was defined'
c              write(*,*) '       in the dependent node'
c              call exit(201)
c            endif
c          endif
!     
          nmpc=nmpc+1
!     
          if(nmpc.gt.nmpc_) then
            write(*,*) '*ERROR reading *EQUATION: increase nmpc_'
            return
          endif
!     
!     MPC is characterized by the label RM (refine mesh)
!     
          labmpc(nmpc)='RM                  '
          ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
          do m=nmpc,id+2,-1
            ikmpc(m)=ikmpc(m-1)
            ilmpc(m)=ilmpc(m-1)
          enddo
          ikmpc(id+1)=idof
          ilmpc(id+1)=nmpc
!     
          nodempc(1,mpcfree)=node
          nodempc(2,mpcfree)=j
          coefmpc(mpcfree)=1.d0
          mpcfree=nodempc(3,mpcfree)           
!     
          do k=1,nterms
            if(dabs(ratio(k)).gt.1.d-10) then
              nodempc(1,mpcfree)=konl(k)
              nodempc(2,mpcfree)=j
              coefmpc(mpcfree)=-ratio(k)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
!     
              if(mpcfree.eq.0) then
                write(*,*) 
     &               '*ERROR reading *EQUATION: increase memmpc_'
                return
              endif
            endif
          enddo
          nodempc(3,mpcfreeold)=0
        enddo
        cycle
      enddo loop1
!     
      return
      end
