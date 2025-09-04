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
      subroutine gen3dmpc(ipompc,nodempc,coefmpc,nmpc,nmpc_,mpcfree,
     &     ikmpc,ilmpc,labmpc,iponoel,inoel,iponoelmax,kon,ipkon,lakon,
     &     ne,iponor,xnor,knor,rig)
!     
!     connects nodes of 1-D and 2-D elements, for which MPC's were
!     defined, to the nodes of their expanded counterparts
!     
      implicit none
!     
      logical dependent
!     
      character*8 lakon(*)
      character*20 labmpc(*)
!     
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,ikmpc(*),
     &     ilmpc(*),iponoel(*),inoel(3,*),iponoelmax,kon(*),ipkon(*),
     &     ne,iponor(2,*),knor(*),rig(*),i,index1,node,index2,ielem,
     &     indexe,j,indexk,newnode,idir,idof,id,mpcfreenew,k,idofold,
     &     idofnew,idold,idnew
!     
      real*8 coefmpc(*),xnor(*)
!     
      do i=1,nmpc
        index1=ipompc(i)
        dependent=.true.
        loop: do
          node=nodempc(1,index1)
          if(node.le.iponoelmax) then
            if(rig(node).ne.0) then
!
!             knot, therefore 1D/2D elements are attached
!             check for MPC types which are not allowed in combination
!             with 1D/2D elements.
!
              if((labmpc(i)(1:8).eq.'STRAIGHT').or.
     &             (labmpc(i)(1:5).eq.'PLANE').or.
     &             (labmpc(i)(1:4).eq.'BEAM').or.
     &             (labmpc(i)(1:8).eq.'MEANROT ').or.
     &             (labmpc(i)(1:6).eq.'CYCLIC')) then
                write(*,*) '*ERROR in gen3dmpc: nodes belonging to'
                write(*,*) '       1D/2D elements must not be'
                write(*,*) '       subject to a ',labmpc(i)(1:20),
     &               ' MPC.'
                call exit(201)
              endif
!              
              if(nodempc(2,index1).gt.3) then
                if(rig(node).lt.0) then
                  write(*,*) '*ERROR in gen3dmpc: in node ',node
                  write(*,*) '  a rotational DOF is constrained'
                  write(*,*) '  by a SPC; however, the elements'
                  write(*,*) '  to which this node belongs do not'
                  write(*,*) '  have rotational DOFs'
                  call exit(201)
                endif
                nodempc(1,index1)=rig(node)
                nodempc(2,index1)=nodempc(2,index1)-3
!     
!     adapting ikmpc and ilmpc (only for the dependent
!     term)
!     
                if(dependent) then
                  idofold=8*(node-1)+nodempc(2,index1)+3
                  call nident(ikmpc,idofold,nmpc,idold)
                  idofnew=8*(rig(node)-1)+nodempc(2,index1)
                  call nident(ikmpc,idofnew,nmpc,idnew)
                  if(idold.le.idnew) then
                    do j=idold,idnew-1
                      ikmpc(j)=ikmpc(j+1)
                      ilmpc(j)=ilmpc(j+1)
                    enddo
                  else
                    do j=idold,idnew+1,-1
                      ikmpc(j)=ikmpc(j-1)
                      ilmpc(j)=ilmpc(j-1)
                    enddo
                  endif
                  ikmpc(idnew)=idofnew
                  ilmpc(idnew)=i
                endif
!     
              endif
            else
              index2=iponoel(node)
!     
!     next loop is needed, since an element may have
!     been deactivated
!     
              do
!     
!     check for nodes not belonging to 1d or 2d elements
!     
                if(index2.eq.0) then
                  index1=nodempc(3,index1)
                  dependent=.false.
                  if(index1.eq.0) exit loop
                  cycle loop
                endif
!
!             node belongs to 1D/2D elements
!             check for MPC types which are not allowed in combination
!             with 1D/2D elements.
!
                if((labmpc(i)(1:8).eq.'STRAIGHT').or.
     &               (labmpc(i)(1:5).eq.'PLANE').or.
     &               (labmpc(i)(1:4).eq.'BEAM').or.
     &               (labmpc(i)(1:8).eq.'MEANROT ')) then
                  write(*,*) '*ERROR in gen3dmpc: nodes belonging to'
                  write(*,*) '       1D/2D elements must not be'
                  write(*,*) '       subject to a ',labmpc(i)(1:20),
     &                 ' MPC.'
                  call exit(201)
                endif
!
              ielem=inoel(1,index2)
                indexe=ipkon(ielem)
                if(indexe.ge.0) exit
                index2=inoel(3,index2)
              enddo
              j=inoel(2,index2)
              indexk=iponor(2,indexe+j)
!     
!     2d element shell element
!     
              if(lakon(ielem)(7:7).eq.'L') then
                if(labmpc(i)(1:5).eq.'RIGID') then
                  write(*,*) '*ERROR in gen3dmpc: nodes belonging to'
                  write(*,*) '       shell elements must not be'
                  write(*,*) '       subject to a ',labmpc(i)(1:20),
     &                 ' MPC.'
                  call exit(201)
                endif
                newnode=knor(indexk+1)
                idir=nodempc(2,index1)
                idof=8*(newnode-1)+idir
                call nident(ikmpc,idof,nmpc,id)
                if((id.le.0).or.(ikmpc(id).ne.idof)) then
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
                  nodempc(1,mpcfree)=knor(indexk+3)
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
                  coefmpc(mpcfree)=-2.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dmpc: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
                endif
!     
!               u(n_2)=u(n) (proposed by victor kemp; causes problems in
!               shellf.inp)    
!     
c                newnode=knor(indexk+2)
c                idof=8*(newnode-1)+idir
c                call nident(ikmpc,idof,nmpc,id)
c                if((id.le.0).or.(ikmpc(id).ne.idof)) then
c                  nmpc=nmpc+1
c                  if(nmpc.gt.nmpc_) then
c                    write(*,*) 
c     &                   '*ERROR in gen3dboun: increase nmpc_'
c                    call exit(201)
c                  endif
c                  labmpc(nmpc)='                    '
c                  ipompc(nmpc)=mpcfree
c                  do j=nmpc,id+2,-1
c                    ikmpc(j)=ikmpc(j-1)
c                    ilmpc(j)=ilmpc(j-1)
c                  enddo
c                  ikmpc(id+1)=idof
c                  ilmpc(id+1)=nmpc
c                  nodempc(1,mpcfree)=newnode
c                  nodempc(2,mpcfree)=idir
c                  coefmpc(mpcfree)=1.d0
c                  mpcfree=nodempc(3,mpcfree)
c                  if(mpcfree.eq.0) then
c                    write(*,*) 
c     &                   '*ERROR in gen3dboun: increase memmpc_'
c                    call exit(201)
c                  endif
c                  nodempc(1,mpcfree)=node
c                  nodempc(2,mpcfree)=idir
c                  coefmpc(mpcfree)=-1.d0
c                  mpcfreenew=nodempc(3,mpcfree)
c                  if(mpcfreenew.eq.0) then
c                    write(*,*) 
c     &                   '*ERROR in gen3dboun: increase memmpc_'
c                    call exit(201)
c                  endif
c                  nodempc(3,mpcfree)=0
c                  mpcfree=mpcfreenew
c                endif
!     
!     fixing the temperature degrees of freedom
!     
                if(idir.eq.0) then
!     
!     t(n_3)=t(n)
!     
                  newnode=knor(indexk+3)
                  idof=8*(newnode-1)+idir
                  call nident(ikmpc,idof,nmpc,id)
                  if((id.le.0).or.(ikmpc(id).ne.idof)) then
                    nmpc=nmpc+1
                    if(nmpc.gt.nmpc_) then
                      write(*,*) 
     &                     '*ERROR in gen3dboun: increase nmpc_'
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
     &                     '*ERROR in gen3dboun: increase memmpc_'
                      call exit(201)
                    endif
                    nodempc(1,mpcfree)=node
                    nodempc(2,mpcfree)=idir
                    coefmpc(mpcfree)=-1.d0
                    mpcfreenew=nodempc(3,mpcfree)
                    if(mpcfreenew.eq.0) then
                      write(*,*) 
     &                     '*ERROR in gen3dboun: increase memmpc_'
                      call exit(201)
                    endif
                    nodempc(3,mpcfree)=0
                    mpcfree=mpcfreenew
                  endif
                endif
              elseif(lakon(ielem)(7:7).eq.'B') then
!     
!     1d beam element
!     
                if(labmpc(i)(1:5).eq.'RIGID') then
                  write(*,*) '*ERROR in gen3dmpc: nodes belonging to'
                  write(*,*) '       beam elements must not be'
                  write(*,*) '       subject to a ',labmpc(i)(1:20),
     &                 ' MPC.'
                  call exit(201)
                endif
                newnode=knor(indexk+1)
                idir=nodempc(2,index1)
                idof=8*(newnode-1)+idir
                call nident(ikmpc,idof,nmpc,id)
                if((id.le.0).or.(ikmpc(id).ne.idof)) then
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
                  do k=2,4
                    nodempc(1,mpcfree)=knor(indexk+k)
                    nodempc(2,mpcfree)=idir
                    coefmpc(mpcfree)=1.d0
                    mpcfree=nodempc(3,mpcfree)
                    if(mpcfree.eq.0) then
                      write(*,*) 
     &                     '*ERROR in gen3dmpc: increase memmpc_'
                      call exit(201)
                    endif
                  enddo
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-4.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dmpc: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
                endif
!     
!     fixing the temperature degrees of freedom
!     
                if(idir.eq.0) then
                  do k=2,4
!     
!     t(n_k)=t(n), k=2,4
!     
                    newnode=knor(indexk+k)
                    idof=8*(newnode-1)+idir
                    call nident(ikmpc,idof,nmpc,id)
                    if((id.le.0).or.(ikmpc(id).ne.idof)) then
                      nmpc=nmpc+1
                      if(nmpc.gt.nmpc_) then
                        write(*,*) 
     &                       '*ERROR in gen3dboun: increase nmpc_'
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
     &                       '*ERROR in gen3dboun: increase memmpc_'
                        call exit(201)
                      endif
                      nodempc(1,mpcfree)=node
                      nodempc(2,mpcfree)=idir
                      coefmpc(mpcfree)=-1.d0
                      mpcfreenew=nodempc(3,mpcfree)
                      if(mpcfreenew.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dboun: increase memmpc_'
                        call exit(201)
                      endif
                      nodempc(3,mpcfree)=0
                      mpcfree=mpcfreenew
                    endif
                  enddo
                endif
              else
!     
!     2d plane stress, plane strain or axisymmetric
!     element
!     
                newnode=knor(indexk+2)
                idir=nodempc(2,index1)
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
              endif
            endif
          endif
          index1=nodempc(3,index1)
          dependent=.false.
          if(index1.eq.0) exit
        enddo loop
      enddo
!     
      return
      end
