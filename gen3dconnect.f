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
      subroutine gen3dconnect(kon,ipkon,lakon,ne,iponoel,inoel,
     &  iponoelmax,rig,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,
     &  nmpc_,mpcfree,ikmpc,ilmpc,labmpc,vold,ikboun,ilboun,nboun,
     &  nboun_,nodeboun,ndirboun,xboun,iamboun,typeboun,ithermal,
     &  mi,trab,ntrans,nmethod,nk,nk_,nam,inotr,iperturb,co)
!
!     connects expanded 1-D and 2-D elements with genuine 3D elements
!     or spring elements or mass elements
!
      implicit none
!
      logical fixed
!
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
!
      integer kon(*),ipkon(*),ne,iponoel(*),inoel(3,*),iponoelmax,
     &  rig(*),iponor(2,*),knor(*),ipompc(*),nodempc(3,*),nmpc,nmpc_,
     &  mpcfree,ikmpc(*),ilmpc(*),i,indexes,nope,l,node,index2,ielem,
     &  indexe,j,indexk,newnode,idir,idof,id,mpcfreenew,k,inotr(2,*),
     &  ikboun(*),ilboun(*),nboun,nboun_,nodeboun(*),ndirboun(*),
     &  iamboun(*),ithermal(*),mi(*),ntrans,nmethod,nk,nk_,nam,
     &  iperturb(*),idummy,iamplitude
!
      real*8 xnor(*),coefmpc(*),vold(0:mi(2),*),val,xboun(*),
     &  trab(7,*),co(3,*)
!
      label='                    '
      fixed=.false.
!
!     generating MPC's to connect shells and beams with solid
!     elements or spring elements or mass elements
!
      do i=1,ne
         indexes=ipkon(i)
         if(indexes.lt.0) cycle
!
!        looking for solid elements or spring elements or
!        mass elements only
!
         if((lakon(i)(7:7).ne.' ').and.(lakon(i)(1:1).ne.'E')) cycle
!
!        determining the number of nodes belonging to the element
!
         if(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(1:1).eq.'E') then
            read(lakon(i)(8:8),'(i1)') nope
            nope=nope+1
         elseif(lakon(i)(1:4).eq.'MASS') then
            nope=1
         else
            cycle
         endif
!
         loop: do l=1,nope
            node=kon(indexes+l)
            if(node.le.iponoelmax) then
               if(rig(node).eq.0) then
                 index2=iponoel(node)
!
!     next loop is needed, since an element may have
!     been deactivated
!
                  do
                    if(index2.eq.0) cycle loop
                    ielem=inoel(1,index2)
                    indexe=ipkon(ielem)
                    if(indexe.ge.0) exit
                    index2=inoel(3,index2)
                  enddo
                  j=inoel(2,index2)
                  indexk=iponor(2,indexe+j)
!
!                 2d shell element: the exterior expanded nodes
!                 are connected to the non-expanded node
!
                  if(lakon(ielem)(7:7).eq.'L') then
                     newnode=knor(indexk+1)
                     do idir=0,3
                        idof=8*(newnode-1)+idir
                        call nident(ikmpc,idof,nmpc,id)
                        if((id.le.0).or.(ikmpc(id).ne.idof)) then
                           nmpc=nmpc+1
                           if(nmpc.gt.nmpc_) then
                              write(*,*) 
     &                          '*ERROR in gen3dconnect: increase nmpc_'
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
     &                        '*ERROR in gen3dconnect: increase memmpc_'
                              call exit(201)
                           endif
                           nodempc(1,mpcfree)=knor(indexk+3)
                           nodempc(2,mpcfree)=idir
                           coefmpc(mpcfree)=1.d0
                           mpcfree=nodempc(3,mpcfree)
                           if(mpcfree.eq.0) then
                              write(*,*) 
     &                        '*ERROR in gen3dconnect: increase memmpc_'
                              call exit(201)
                           endif
                           nodempc(1,mpcfree)=node
                           nodempc(2,mpcfree)=idir
                           coefmpc(mpcfree)=-2.d0
                           mpcfreenew=nodempc(3,mpcfree)
                           if(mpcfreenew.eq.0) then
                              write(*,*) 
     &                        '*ERROR in gen3dconnect: increase memmpc_'
                              call exit(201)
                           endif
                           nodempc(3,mpcfree)=0
                           mpcfree=mpcfreenew
                        endif
                     enddo
                  elseif(lakon(ielem)(7:7).eq.'B') then
!
!                    1d beam element: corner nodes are connected to
!                    the not-expanded node
!
                     newnode=knor(indexk+1)
                     do idir=0,3
                        idof=8*(newnode-1)+idir
                        call nident(ikmpc,idof,nmpc,id)
                        if((id.le.0).or.(ikmpc(id).ne.idof)) then
                           nmpc=nmpc+1
                           if(nmpc.gt.nmpc_) then
                              write(*,*) 
     &                          '*ERROR in gen3dconnect: increase nmpc_'
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
     &                        '*ERROR in gen3dconnect: increase memmpc_'
                              call exit(201)
                           endif
                           do k=2,4
                              nodempc(1,mpcfree)=knor(indexk+k)
                              nodempc(2,mpcfree)=idir
                              coefmpc(mpcfree)=1.d0
                              mpcfree=nodempc(3,mpcfree)
                              if(mpcfree.eq.0) then
                                 write(*,*) 
     &                        '*ERROR in gen3dconnect: increase memmpc_'
                                 call exit(201)
                              endif
                           enddo
                           nodempc(1,mpcfree)=node
                           nodempc(2,mpcfree)=idir
                           coefmpc(mpcfree)=-4.d0
                           mpcfreenew=nodempc(3,mpcfree)
                           if(mpcfreenew.eq.0) then
                              write(*,*) 
     &                        '*ERROR in gen3dconnect: increase memmpc_'
                              call exit(201)
                           endif
                           nodempc(3,mpcfree)=0
                           mpcfree=mpcfreenew
                        endif
                     enddo
                  else
!     
!                    2d plane stress, plane strain or axisymmetric
!                    element: the expanded middle node (this is the
!                    "governing" node for these elements) is connected to
!                    the non-expanded node
!
                     newnode=knor(indexk+2)
                     do idir=0,2
                        idof=8*(newnode-1)+idir
                        call nident(ikmpc,idof,nmpc,id)
                        if((id.le.0).or.(ikmpc(id).ne.idof)) then
                           nmpc=nmpc+1
                           if(nmpc.gt.nmpc_) then
                              write(*,*) 
     &                          '*ERROR in gen3dconnect: increase nmpc_'
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
     &                        '*ERROR in gen3dconnect: increase memmpc_'
                              call exit(201)
                           endif
                           nodempc(1,mpcfree)=node
                           nodempc(2,mpcfree)=idir
                           coefmpc(mpcfree)=-1.d0
                           mpcfreenew=nodempc(3,mpcfree)
                           if(mpcfreenew.eq.0) then
                              write(*,*) 
     &                        '*ERROR in gen3dconnect: increase memmpc_'
                              call exit(201)
                           endif
                           nodempc(3,mpcfree)=0
                           mpcfree=mpcfreenew
                        endif
                     enddo
!
!                    fixing the original node in z-direction
!
!                    since the original node belongs to a 3-D element
!                    (or a spring...) all dofs in this node are
!                    active and the dof in the z-direction has to be
!                    fixed (can lead to a singular matrix, e.g. is all
!                    elements are plane strain and there is one spring1
!                    element).
!
                     if(ithermal(2).ne.2) then
                        val=0.d0
                        k=3
                        if(nam.gt.0) iamplitude=0
                        type='M'
                        call bounadd(node,k,k,val,nodeboun,
     &                       ndirboun,xboun,nboun,nboun_,iamboun,
     &                       iamplitude,nam,ipompc,nodempc,coefmpc,
     &                       nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &                       ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,
     &                       labmpc,type,typeboun,nmethod,iperturb,
     &                       fixed,vold,idummy,mi,label)
                     endif
!     
                  endif
               endif
            endif
         enddo loop
      enddo
!
      return
      end


