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
      subroutine gen3dsurf(iponoel,inoel,iponoelmax,kon,ipkon,
     &     lakon,ne,iponor,knor,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &     mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,nk_,
     &     co,nmethod,iperturb,nset,set,istartset,iendset,ialset,
     &     ikboun,ilboun,nboun,nboun_,nodeboun,ndirboun,xboun,iamboun,
     &     typeboun,mi,vold)
!     
!     connects nodes of 1-D and 2-D elements, belonging to nodal
!     surfaces, to the nodes of their expanded counterparts
!     
      implicit none
!     
      logical fixed
!     
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
      character*81 set(*)
!     
      integer iponoel(*),inoel(3,*),iponoelmax,kon(*),ipkon(*),ne,
     &     iponor(2,*),knor(*),ipompc(*),nodempc(3,*),nmpc,nmpc_,
     &     ikmpc(*),ilmpc(*),rig(*),ntrans,inotr(2,*),i,node,mpcfree,
     &     indexx,ielem,j,indexe,indexk,idir,nk,nk_,iamplitude,
     &     newnode,idof,id,mpcfreenew,k,nam,nmethod,iperturb(*),
     &     iendset(*),ialset(*),nset,ipos,l,idummy,ikboun(*),ilboun(*),
     &     nboun,nboun_,nodeboun(*),ndirboun(*),iamboun(*),mi(*),
     &     istartset(*)
!     
      real*8 coefmpc(*),trab(7,*),co(3,*),val,xboun(*),vold(0:mi(2),*)
!     
      label='                    '
      fixed=.false.
!     
      do i=1,nset
        ipos=index(set(i),' ')
!     
!     only nodal surfaces are involved
!     
        if(set(i)(ipos-1:ipos-1).ne.'S') cycle
        loop: do l=istartset(i),iendset(i)
          node=ialset(l)
          if(node.gt.iponoelmax) cycle
          indexx=iponoel(node)
!
!     next loop is needed, since an element may have
!     been deactivated
!
          do
            if(indexx.eq.0) cycle loop
            ielem=inoel(1,indexx)
            indexe=ipkon(ielem)
            if(indexe.ge.0) exit
            indexx=inoel(3,indexx)
          enddo
c          if(indexx.eq.0) cycle
c          ielem=inoel(1,indexx)
          j=inoel(2,indexx)
c          indexe=ipkon(ielem)
          indexk=iponor(2,indexe+j)
!     
          if(rig(node).eq.0) then
!     
!     2d element shell element: generate MPC's
!     
            if(lakon(ielem)(7:7).eq.'L') then
              newnode=knor(indexk+1)
              do idir=1,3
                idof=8*(newnode-1)+idir
                call nident(ikmpc,idof,nmpc,id)
                if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase nmpc_'
                    call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                    ikmpc(j)=ikmpc(j-1)
                    ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(1,mpcfree)=knor(indexk+3)
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-2.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
                endif
              enddo
            elseif(lakon(ielem)(7:7).eq.'B') then
!     
!     1d beam element: generate MPC's
!     
              newnode=knor(indexk+1)
              do idir=1,3
                idof=8*(newnode-1)+idir
                call nident(ikmpc,idof,nmpc,id)
                if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase nmpc_'
                    call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                    ikmpc(j)=ikmpc(j-1)
                    ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase memmpc_'
                    call exit(201)
                  endif
                  do k=2,4
                    nodempc(1,mpcfree)=knor(indexk+k)
                    nodempc(2,mpcfree)=idir
                    coefmpc(mpcfree)=1.d0
                    mpcfree=nodempc(3,mpcfree)
                    if(mpcfree.eq.0) then
                      write(*,*) 
     &                     '*ERROR in gen3dboun: increase memmpc_'
                      call exit(201)
                    endif
                  enddo
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-4.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
                endif
!     
              enddo
            else
!     
!     2d plane stress, plane strain or axisymmetric
!     element: dependent node is replaced by new node in the middle
!     
!     keeping the old node and generating an additional MPC leads
!     to problems since the old node is not restricted in the 
!     z-direction, only the new node in the middle is. If the old
!     node is used subsequently in a contact spring element all
!     its DOFs are used and the unrestricted z-DOF leads to a 
!     singular equation system
!     
c     write(*,*) ialset(l),' replaced by ',knor(indexk+2)
c     co(1,knor(indexk+2))=co(1,ialset(l))
c     co(2,knor(indexk+2))=co(2,ialset(l))
c     co(3,knor(indexk+2))=co(3,ialset(l))
c     ialset(l)=knor(indexk+2)
!     
!     2d plane stress, plane strain or axisymmetric
!     element: MPC in all but z-direction
!     
              newnode=knor(indexk+2)
              do idir=1,2
                idof=8*(newnode-1)+idir
                call nident(ikmpc,idof,nmpc,id)
                if(((id.le.0).or.(ikmpc(id).ne.idof)).and.
     &               (idir.ne.3)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                    write(*,*) 
     &                   '*ERROR in gen3dmpc: increase nmpc_'
                    call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                    ikmpc(j)=ikmpc(j-1)
                    ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dmpc: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-1.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dmpc: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
                endif
              enddo
!     
!     fixing the original node in z-direction
!     
              val=0.d0
              k=3
              if(nam.gt.0) iamplitude=0
              type='M'
              call bounadd(node,k,k,val,nodeboun,
     &             ndirboun,xboun,nboun,nboun_,iamboun,
     &             iamplitude,nam,ipompc,nodempc,coefmpc,
     &             nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &             ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,
     &             labmpc,type,typeboun,nmethod,iperturb,
     &             fixed,vold,idummy,mi,label)
!     
            endif
          endif
        enddo loop
      enddo
!     
      return
      end
      
      
